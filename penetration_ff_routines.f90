!************************************************
! this module contains routines for using a short-range force
! field consisting entirely of buckingham repusive terms to 
! treat charge penetration regions.  This force field is continuously
! pieced together with the decomposed long range force field
! currently, there is one available option for how to merge these force fields.
!
! the first option merges them at an energy threshold (usually zero).
! while this switching distance should be read in, this module checks everything for 
! consistency by using the routine rtsafe to find roots which give this switching distance
! for each atom pair
!
! the second option is currently disabled, as there is no way to make this differentiable
! with the current functional form. 
! old option switches them at the energy minimum, and the amoeba routine is used to find
! the pairwise minima and thus the distances for switching
!****************************************************


module penetration_ff_routines
  use routines

contains

  !*******************************************
  ! this subroutine initializes the use of the penetration force field
  !*******************************************
  Subroutine initialize_penetration_ff_routines
    use global_variables
    integer :: i_atom, j_atom,iter
    real*8,parameter:: r_pen_tol = 1d-8,r_pen_min = 2d0 , r_pen_max = 8d0,energy_thresh=1d-6
    real*8 :: lr_energy, sr_energy , junk, small=.1d0, d_thresh = 1d-5 , e_thresh = 1d-6 , f_thresh = 1d-3, test_distance, p_energy,p_energy_delta, delta=1d-4, p_deriv, l_deriv
    real*8,dimension(2)::y
    real*8,dimension(2,1)::p


    !******************* disable switch at minimum option *******************
    ! rather than delete the working code we already wrote for this option for the 
    ! unlikely event that is it used in the future, here we disable the option
    Select Case(switch_at_minimum)
    Case("no")
100    continue
    case default
       stop "parameter 'switch_at_minimum' must be set to 'no' in the current code"
    End Select
    !************************************************************************


    !************************* old stuff **********************************
    Select Case(switch_at_minimum)
    Case("yes")

       ! find switching distance for all solute-framework atom types
       ! since we only use penetration force field for solute-framework interactions,
       ! set solute-solute switching distance to zero

       atype_penetration_ff(:,:,2)=0d0

       do i_atom = 1, n_atom_type
          do j_atom = 1, n_atom_type
             ! if solute-framework interaction
             if ( ( (atype_solute(i_atom) .eq. 0 ) .and. (atype_solute(j_atom) .eq. 1 ) ) .or. ( (atype_solute(i_atom) .eq. 1 ) .and. (atype_solute(j_atom) .eq. 0 ) ) ) then
                p(1,1)=2d0
                p(2,1)=5d0
                call two_atom_penetration_energy(p(1,1),y(1),junk,i_atom,j_atom)
                call two_atom_penetration_energy(p(2,1),y(2),junk,i_atom,j_atom)
                ! get distance at minimum
                call amoeba(p,y,i_atom,j_atom,r_pen_tol,iter)
                atype_penetration_ff(i_atom,j_atom,2) = p(1,1)
                ! store shift
                sr_energy = atype_penetration_ff(i_atom,j_atom,1) * exp ( - atype_lj_parameter(i_atom,j_atom,2) * atype_penetration_ff(i_atom,j_atom,2) )
                atype_penetration_ff(i_atom,j_atom,3) = sr_energy - y(1)
             endif
          enddo
       enddo

    End Select
    !**************************************************************************



    ! We are a using a piecewise force field with continuous derivatives, all the necessary parameters including the shift distance have been read in. Here, just check for consistency
    do i_atom = 1, n_atom_type
       do j_atom = 1, n_atom_type
          ! if we are not using a piecewise potential for this atom pair, then switch distance will have been set to zero
          if ( atype_penetration_ff(i_atom,j_atom,2) > small ) then

             test_distance = rtsafe( r_pen_min , r_pen_max , r_pen_tol , i_atom , j_atom ) 
             if ( abs( atype_penetration_ff(i_atom,j_atom,2) - test_distance ) > d_thresh ) then
                write(*,*) "piecewise force field between atom types ", atype_name(i_atom) , atype_name(j_atom), " has inconsistent switch distance."
                write(*,*) "input distance was ",atype_penetration_ff(i_atom,j_atom,2), " calculated distance is ",test_distance
!!$                   stop
             endif


             ! now check continuity of force field
             p_energy = atype_penetration_ff(i_atom,j_atom,1) * exp(-atype_penetration_ff(i_atom,j_atom,4) * atype_lj_parameter(i_atom,j_atom,2)* atype_penetration_ff(i_atom,j_atom,2) ) - atype_penetration_ff(i_atom,j_atom,3)
             p_energy_delta = atype_penetration_ff(i_atom,j_atom,1) * exp(-atype_penetration_ff(i_atom,j_atom,4) * atype_lj_parameter(i_atom,j_atom,2)* ( atype_penetration_ff(i_atom,j_atom,2) - delta ) ) - atype_penetration_ff(i_atom,j_atom,3)
             p_deriv = ( p_energy - p_energy_delta ) /delta
             ! energy continuity
             if ( abs( p_energy - penetration_force_field_threshold ) > e_thresh ) then
                write(*,*) "two body penetration force field is not continuous for atom types " , atype_name(i_atom) , atype_name(j_atom)
                write(*,*) "penetration ff energy at shift distance is ",  p_energy
!!$                   stop
             endif
             call two_atom_penetration_energy( atype_penetration_ff(i_atom, j_atom, 2) , lr_energy , l_deriv , i_atom , j_atom ) 
             ! derivative continuity
             if ( abs( p_deriv - l_deriv ) > f_thresh ) then
                write(*,*) "two body penetration force field derivative is not continuous for atom types " , atype_name(i_atom) , atype_name(j_atom)
                write(*,*) "derivative from the left (penetration ff) is ", p_deriv, "derivative from teh right (decomposed ff ) is " , l_deriv
!!$                   stop
             endif
          endif
       enddo
    enddo



  end  subroutine initialize_penetration_ff_routines



!****************************************
! this subroutine returns the two body energy between atom types
! as meant to be called by rtsafe.  We are determining the distance at which
! this two body term equals penetration_force_field_threshold
! since we are finding roots with rtsafe, subtract this term
!****************************************
Subroutine two_atom_penetration_energy(dist,energy,denergy,atom1,atom2)
  use global_variables
  real*8,intent(in) :: dist
  integer,intent(in) :: atom1,atom2
  real*8,intent(out) :: energy,denergy


  real*8 :: bkghm_energy, bkghm_denergy,r_ij,r_ij2,r_ij6,r_ij8,r_ij10,r_ij12,term6,term8,term10,term12,fac6,fac8,fac10,fac12
  real*8,dimension(3) :: r_vec,dterm6,dterm8,dterm10,dterm12

  ! energy
  bkghm_energy = atype_lj_parameter(atom1,atom2,1)*exp(-atype_lj_parameter(atom1,atom2,2)*dist)

  r_ij = dist
  r_ij2 = r_ij**2
  r_ij6 = r_ij2**3
  r_ij8=r_ij6*r_ij2
  r_ij10=r_ij8*r_ij2
  term6= - C6_C10_damp(atom1,atom2,r_ij,6) * atype_lj_parameter(atom1,atom2,3)/r_ij6
  term8= - C6_C10_damp(atom1,atom2,r_ij,8) * atype_lj_parameter(atom1,atom2,4)/r_ij8
  term10= - C6_C10_damp(atom1,atom2,r_ij,10) * atype_lj_parameter(atom1,atom2,5)/r_ij10

  energy = bkghm_energy + term6 + term8 + term10

  Select Case(C12_dispersion)
     Case("yes")
   r_ij12=r_ij10*r_ij2  
   term12= - C6_C10_damp(atom1,atom2,r_ij,12) * atype_lj_parameter(atom1,atom2,6)/r_ij12     
   energy = energy + term12
  End Select


  Select Case(switch_at_minimum)
     Case("no")
  ! subtract threshold energy 
  energy = energy - penetration_force_field_threshold

  ! derivative
  bkghm_denergy = -atype_lj_parameter(atom1,atom2,2) * atype_lj_parameter(atom1,atom2,1)*exp(-atype_lj_parameter(atom1,atom2,2)*dist)
  r_vec=0d0
  r_vec(1)=dist
  fac6= atype_lj_parameter(atom1,atom2,3)/r_ij6
  fac8= atype_lj_parameter(atom1,atom2,4)/r_ij8
  fac10= atype_lj_parameter(atom1,atom2,5)/r_ij10
  dterm6 =  C6_C10_damp(atom1,atom2,r_ij,6) * fac6 * (6d0 *r_vec/r_ij2) - dC6_C10_damp(atom1,atom2,r_vec,r_ij,6) * fac6
  dterm8 =  C6_C10_damp(atom1,atom2,r_ij,8) * fac8 * (8d0 *r_vec/r_ij2) - dC6_C10_damp(atom1,atom2,r_vec,r_ij,8) * fac8
  dterm10 = C6_C10_damp(atom1,atom2,r_ij,10) * fac10 * (10d0 *r_vec/r_ij2) - dC6_C10_damp(atom1,atom2,r_vec,r_ij,10) * fac10

  denergy = bkghm_denergy + dterm6(1) + dterm8(1) + dterm10(1)

  Select Case(C12_dispersion)
     Case("yes")
   fac12= atype_lj_parameter(atom1,atom2,6)/r_ij12    
   dterm12 = C6_C10_damp(atom1,atom2,r_ij,12) * fac12 * (12d0 *r_vec/r_ij2) - dC6_C10_damp(atom1,atom2,r_vec,r_ij,12) * fac12
   denergy = denergy + dterm12(1)
  End Select

  Case("yes")
     ! derivative will not be used here, just set to zero
     denergy=0d0
  End Select

end subroutine two_atom_penetration_energy





!************************************ Numerical Recipes Code ***********************************


!*********************************************
! this is a bracketed newton-raphson root finding function
! from numerical recipes
!*********************************************
FUNCTION rtsafe(x1,x2,xacc,atom1,atom2)
IMPLICIT NONE
REAL*8, INTENT(IN) :: x1,x2,xacc
integer,intent(in) :: atom1,atom2
REAL*8 :: rtsafe
INTEGER, PARAMETER :: MAXIT=100
INTEGER :: j
REAL*8 :: df,dx,dxold,f,fh,fl,temp,xh,xl
call two_atom_penetration_energy(x1,fl,df,atom1,atom2)
call two_atom_penetration_energy(x2,fh,df,atom1,atom2)
if ((fl > 0.0 .and. fh > 0.0) .or. (fl < 0.0 .and. fh < 0.0)) then
 stop 'root must be bracketed in rtsafe'
endif
if (fl == 0.0) then
 rtsafe=x1
 RETURN
else if (fh == 0.0) then
 rtsafe=x2
 RETURN
else if (fl < 0.0) then
 xl=x1
 xh=x2
else
 xh=x1
 xl=x2
end if
rtsafe=0.5d0*(x1+x2)
dxold=abs(x2-x1)
dx=dxold
call two_atom_penetration_energy(rtsafe,f,df,atom1,atom2)
do j=1,MAXIT
 if (((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f) > 0.0 .or. &
      abs(2.0d0*f) > abs(dxold*df) ) then
    dxold=dx
    dx=0.5d0*(xh-xl)
    rtsafe=xl+dx
    if (xl == rtsafe) RETURN
 else
    dxold=dx
    dx=f/df
    temp=rtsafe
    rtsafe=rtsafe-dx
    if (temp == rtsafe) RETURN
 end if
 if (abs(dx) < xacc) RETURN
 call two_atom_penetration_energy(rtsafe,f,df,atom1,atom2)
 if (f < 0.0) then
    xl=rtsafe
 else
    xh=rtsafe
 end if
end do

stop 'rtsafe: exceeded maximum iterations'

END FUNCTION rtsafe


!************************************
! this finds minima, from numerical recipes, slightly modified
!************************************
	SUBROUTINE amoeba(p,y,atom1,atom2,ftol,iter)
	IMPLICIT NONE
	INTEGER, INTENT(OUT) :: iter
	REAL*8, INTENT(IN) :: ftol
        integer,intent(in) :: atom1,atom2
	REAL*8, DIMENSION(:), INTENT(INOUT) :: y
	REAL*8, DIMENSION(:,:), INTENT(INOUT) :: p
	INTEGER, PARAMETER :: ITMAX=100
	REAL*8, PARAMETER :: TINY=1.0e-10
	INTEGER :: ihi,ndim
	REAL*8, DIMENSION(size(p,2)) :: psum
	INTEGER :: i,ilo,inhi
	REAL*8 :: rtol,ysave,ytry,ytmp,df
	ndim=size(p,2)
        if ( (size(p,2) .ne. size(p,1)-1 ) .or. (size(p,2) .ne. size(y)-1) ) then
           stop "error in input array sizes in amoeba"
        endif


        iter=0
	psum(:)=sum(p(:,:),dim=1)
	do
       
           ilo=iminloc(y(:))
		ihi=imaxloc(y(:))
		ytmp=y(ihi)
		y(ihi)=y(ilo)
		inhi=imaxloc(y(:))
		y(ihi)=ytmp
		rtol=2.0d0*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
		if (rtol < ftol) then
			call swap(y(1),y(ilo))
			call swap2(p(1,:),p(ilo,:))
			RETURN
		end if
		if (iter >= ITMAX) then
                   stop "ITMAX exceeded in amoeba"
                endif
		ytry=amotry(-1.0d0,ihi,p,psum,y,atom1,atom2)
		iter=iter+1
		if (ytry <= y(ilo)) then
			ytry=amotry(2.0d0,ihi,p,psum,y,atom1,atom2)
			iter=iter+1
		else if (ytry >= y(inhi)) then
			ysave=y(ihi)
			ytry=amotry(0.5d0,ihi,p,psum,y,atom1,atom2)
			iter=iter+1
			if (ytry >= ysave) then
				p(:,:)=0.5d0*(p(:,:)+spread(p(ilo,:),1,size(p,1)))
				do i=1,ndim+1
					if (i /= ilo) then
                                        call two_atom_penetration_energy(p(i,1),y(i),df,atom1,atom2)
                                        endif
				end do
				iter=iter+ndim
				psum(:)=sum(p(:,:),dim=1)
			end if
		end if
	end do
	END SUBROUTINE amoeba


	FUNCTION amotry(fac,ihi,p,psum,y,atom1,atom2)
	IMPLICIT NONE
	REAL*8, INTENT(IN) :: fac
        integer, intent(in) :: ihi,atom1,atom2
	REAL*8, DIMENSION(:,:), INTENT(INOUT) :: p
        real*8,dimension(:),intent(inout) :: psum,y
	REAL*8 :: amotry
	REAL*8 :: fac1,fac2,ytry,df
	REAL*8, DIMENSION(size(p,2)) :: ptry
        integer :: ndim

	ndim=size(p,2)
	fac1=(1.0d0-fac)/ndim
	fac2=fac1-fac
	ptry(:)=psum(:)*fac1-p(ihi,:)*fac2
        call two_atom_penetration_energy(ptry(1),ytry,df,atom1,atom2)
        	if (ytry < y(ihi)) then
		y(ihi)=ytry
		psum(:)=psum(:)-p(ihi,:)+ptry(:)
		p(ihi,:)=ptry(:)
	end if
	amotry=ytry
	END FUNCTION amotry


end module penetration_ff_routines
