MODULE f1dim_mod
  INTEGER   :: ncom
  real*8,dimension(:),pointer::pcom
  REAL*8, DIMENSION(:), POINTER :: xicom
CONTAINS
  !BL
  FUNCTION f1dim(x,chg,mass,box,tot_n_mole,n_mole,n_atom,dfti_desc,dfti_desc_inv,log_file)
    use MKL_DFTI
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: x
    integer,intent(in)::tot_n_mole,n_mole, log_file
    integer,dimension(:),intent(inout)::n_atom
    real*8, dimension(:,:),intent(inout) :: chg,mass
    real*8, dimension(3,3),intent(in)::box
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv
    REAL*8 :: f1dim
    integer::i,j
    INTERFACE
       function kaisq(p,chg,mass,box,tot_n_mole,n_mole,n_atom,dfti_desc,dfti_desc_inv,log_file)
         use MKL_DFTI
         use global_variables
         use routines
         use ms_evb
         real*8  :: kaisq
         real*8,dimension(:), intent(in) :: p
         integer,intent(in)::tot_n_mole,n_mole, log_file
         integer,dimension(:),intent(inout)::n_atom
         real*8, dimension(:,:),intent(inout) :: chg,mass
         real*8, dimension(3,3),intent(in)::box
         TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv
       end function kaisq
    END INTERFACE
    real*8,dimension(size(pcom))::xt
    xt=pcom+(x*xicom)
    f1dim=kaisq(xt,chg,mass,box,tot_n_mole,n_mole,n_atom,dfti_desc,dfti_desc_inv,log_file)

!!$    write(*,*) "in f1dim", f1dim
!!$    write(*,*) xt
!!$    write(*,*) "leaving f1dim"

  END FUNCTION f1dim
END MODULE f1dim_mod



module conj_grad

contains

  !******************************************
  ! this subroutine finds the nearest local minimum energy
  ! configuration of the system through a conjugate gradient
  ! optimization of the atomic coordinates
  !
  ! a call to the multi-state empirical valence bond module
  ! ( ms_evb_calculate_total_force_energy ) may change data structures
  ! chg, n_atom, mass, if a proton transfer occurs
  ! therefore these data structures are given the intent(inout) attribute
  !******************************************
  subroutine conjugate_gradient_minimize(box,tot_n_mole,n_mole,n_atom,n_atom_drude,xyz,xyz_drude,mass,potential,chg,dfti_desc,dfti_desc_inv,log_file)
    use MKL_DFTI
    use ms_evb
    use routines
    use global_variables
    integer,intent(in)::n_mole,tot_n_mole,log_file
    integer,dimension(:),intent(inout)::n_atom,n_atom_drude
    real*8,dimension(:,:,:),intent(inout)::xyz,xyz_drude
    real*8,dimension(:,:),intent(inout)::mass,chg
    real*8,dimension(:,:),intent(in)::box
    real*8,intent(in)::potential
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv

    real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: force
    real*8, dimension(MAX_N_MOLE,3) :: r_com
    character(MAX_ANAME) :: atomname
    character(MAX_MNAME) :: moleculename
    real*8  :: energy
    real*8, dimension(:), allocatable :: p
    real*8 :: ftol=1d-4, fret
    integer :: i_mole,i_atom, i,ndim, index, iteration, count


    ! set up p, dimension is 3 * N
    ndim = 3 * total_atoms
    allocate(p(ndim))

    ! initialize one dimensional array of coordinates
    index=1 
    do i_mole=1,n_mole
       do i_atom=1,n_atom(i_mole)
          do i=1,3
             p(index) = xyz(i_mole, i_atom, i )
             index = index + 1
          enddo
       enddo
    enddo

    !******************************
    ! this is numerical recipes conjugate gradient minimizer, see numerical recipes for details
    !******************************

    write(*,*) "starting conjugate gradient optimization..."
    call frprmn(p,ftol,iteration,fret,chg,mass,box,tot_n_mole,n_mole,n_atom,dfti_desc,dfti_desc_inv,log_file)

    ! map data structures back
    index=1 
    do i_mole=1,n_mole
       do i_atom=1,n_atom(i_mole)
          do i=1,3
             xyz(i_mole,i_atom,i) = p(index)
             index = index + 1
          enddo
       enddo
    enddo
  ! get center of mass for these coordinates
  do i_mole=1,n_mole
     ! after updating atomic coordinates, calculate new center of mass of molecules
     r_com(i_mole,:) = pos_com( xyz, i_mole, n_atom, mass )
     ! translate molecules back into the box
     call shift_move(n_atom(i_mole),xyz(i_mole,:,:),r_com(i_mole,:),box)
  end do


    ! get final energy
    call ms_evb_calculate_total_force_energy( force, energy, tot_n_mole, n_mole, n_atom, xyz, r_com, chg, mass, box, dfti_desc,dfti_desc_inv,log_file )

    write(*,*) "conjugate gradient minimization finished after"
    write(*,*) iteration, " steps"
    write(*,*) "optimized coordinates"
    count=0
    do i_mole=1,n_mole
       moleculename = molecule_type_name(molecule_index(i_mole))
       do i_atom = 1, n_atom(i_mole)
          count = count + 1
          atomname = atype_name(atom_index(i_mole,i_atom))
          write( *, '(I5,2A5,I5,3F14.6)' ) i_mole, moleculename, atomname, count, xyz(i_mole,i_atom,:)
       enddo
    enddo

    write(*,*) ""
    write(*,*) "optimized energy", energy
    write(*,*) ""
    write(*,*) "forces on atoms"
    do i_mole=1,n_mole
       moleculename = molecule_type_name(molecule_index(i_mole))
       do i_atom = 1, n_atom(i_mole)
          count = count + 1
          atomname = atype_name(atom_index(i_mole,i_atom))
          write( *, '(I5,2A5,I5,3F14.6)' ) i_mole, moleculename, atomname, count, force(i_mole,i_atom,:)
       enddo
    enddo

    stop

  end subroutine conjugate_gradient_minimize




  !*********************************************
  ! frprmn is the main conjugate gradient fitting subroutine
  ! used in the fitting program
  ! 
  ! the rest of the subroutines in this file are called by frprmn
  ! 
  ! all of these subroutines are adopted from numerical recipes
  !*********************************************

  SUBROUTINE frprmn(p,ftol,iter,fret,chg,mass,box,tot_n_mole,n_mole,n_atom,dfti_desc,dfti_desc_inv,log_file)
    use MKL_DFTI
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: iter
    REAL*8, INTENT(IN) :: ftol
    REAL*8, INTENT(OUT) :: fret
    real*8,dimension(:), INTENT(INOUT) :: p
    integer,intent(in)::tot_n_mole,n_mole, log_file
    integer,dimension(:),intent(inout)::n_atom
    real*8, dimension(:,:),intent(inout) :: chg,mass
    real*8, dimension(3,3),intent(in)::box
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv
    INTERFACE
       function kaisq(p,chg,mass,box,tot_n_mole,n_mole,n_atom,dfti_desc,dfti_desc_inv,log_file)
         use MKL_DFTI
         use global_variables
         use routines
         use ms_evb
         real*8  :: kaisq
         real*8,dimension(:), intent(in) :: p
         integer,intent(in)::tot_n_mole,n_mole, log_file
         integer,dimension(:),intent(inout)::n_atom
         real*8, dimension(:,:),intent(inout) :: chg,mass
         real*8, dimension(3,3),intent(in)::box
         TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv
       end function kaisq
       function dkaisq(p,chg,mass,box,tot_n_mole,n_mole,n_atom,dfti_desc,dfti_desc_inv,log_file)
         use MKL_DFTI
         use global_variables
         use routines
         use ms_evb
         real*8,dimension(:), intent(in) :: p
         REAL*8, DIMENSION(size(p)) :: dkaisq
         integer,intent(in)::tot_n_mole,n_mole, log_file
         integer,dimension(:),intent(inout)::n_atom
         real*8, dimension(:,:),intent(inout) :: chg,mass
         real*8, dimension(3,3),intent(in)::box
         TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv
       end function dkaisq
    END INTERFACE
    INTEGER, PARAMETER :: ITMAX=10000
    REAL*8, PARAMETER :: EPS=1.0D-10
    INTEGER :: its
    REAL*8 :: dgg,fp,gam,gg
    REAL*8, DIMENSION(size(p)) :: g,h,xi
    ! test
    integer :: index, i_mole, i_atom
    ! test

    fp=kaisq(p,chg,mass,box,tot_n_mole,n_mole,n_atom,dfti_desc,dfti_desc_inv,log_file)
    xi=dkaisq(p,chg,mass,box,tot_n_mole,n_mole,n_atom,dfti_desc,dfti_desc_inv,log_file)

!!$        write(*,*) "kaisq",fp
!!$        write(*,*) "dkaisq"
    g=-xi
    h=g
    xi=h
    do its=1,ITMAX
       iter=its

       call linmin(p,xi,fret,chg,mass,box,tot_n_mole,n_mole,n_atom,dfti_desc,dfti_desc_inv,log_file)

       ! test
!!$       write(*,*) "iteration ", iter
!!$       index=1
!!$       do i_mole=1,n_mole
!!$          do i_atom=1,n_atom(i_mole)
!!$             write(*,*) i_mole, i_atom, p(index)
!!$             index=index+1
!!$          enddo
!!$       enddo
       !test

       if (2.0D0*abs(fret-fp) <= ftol*(abs(fret)+abs(fp)+EPS)) RETURN
       fp=fret
       xi=dkaisq(p,chg,mass,box,tot_n_mole,n_mole,n_atom,dfti_desc,dfti_desc_inv,log_file)
       gg=dot_product(g,g)
       !		dgg=dot_product(xi,xi)
       dgg=dot_product(xi+g,xi)
       if (gg == 0.0) RETURN
       gam=dgg/gg
       g=-xi
       h=g+gam*h
       xi=h
    end do
    call nrerror('frprmn: maximum iterations exceeded')
  END SUBROUTINE frprmn


  SUBROUTINE linmin(p,xi,fret,chg,mass,box,tot_n_mole,n_mole,n_atom,dfti_desc,dfti_desc_inv,log_file)
    USE f1dim_mod
    use MKL_DFTI
    IMPLICIT NONE
    REAL*8, INTENT(OUT) :: fret
    real*8,dimension(:),target,intent(inout)::p
    REAL*8, DIMENSION(:), TARGET, INTENT(INOUT) :: xi
    integer,intent(in)::tot_n_mole,n_mole, log_file
    integer,dimension(:),intent(inout)::n_atom
    real*8, dimension(:,:),intent(inout) :: chg,mass
    real*8, dimension(3,3),intent(in)::box
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv
    REAL*8, PARAMETER :: TOL=1.0D-4
    REAL*8 :: ax,bx,fa,fb,fx,xmin,xx
    integer::i
    ncom=size(xi)
    pcom=>p
    xicom=>xi
    ax=0.0
   
    !xx=1.0
    ! changed by Jesse 3/10/2015 
    xx=0.001d0


    call mnbrak(ax,xx,bx,fa,fx,fb,f1dim,chg,mass,box,tot_n_mole,n_mole,n_atom,dfti_desc,dfti_desc_inv,log_file)
    fret=brent(ax,xx,bx,f1dim,TOL,xmin,chg,mass,box,tot_n_mole,n_mole,n_atom,dfti_desc,dfti_desc_inv,log_file)
    xi=xmin*xi
    p=p+xi
  END SUBROUTINE linmin


  SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func,chg,mass,box,tot_n_mole,n_mole,n_atom,dfti_desc,dfti_desc_inv,log_file)
    use MKL_DFTI
    IMPLICIT NONE
    REAL*8, INTENT(INOUT) :: ax,bx
    REAL*8, INTENT(OUT) :: cx,fa,fb,fc
    integer,intent(in)::tot_n_mole,n_mole, log_file
    integer,dimension(:),intent(inout)::n_atom
    real*8, dimension(:,:),intent(inout) :: chg,mass
    real*8, dimension(3,3),intent(in)::box
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv
    INTERFACE
       FUNCTION func(x,chg,mass,box,tot_n_mole,n_mole,n_atom,dfti_desc,dfti_desc_inv,log_file)
         use MKL_DFTI
         IMPLICIT NONE
         integer,intent(in)::tot_n_mole,n_mole, log_file
         integer,dimension(:),intent(inout)::n_atom
         real*8, dimension(:,:),intent(inout) :: chg,mass
         real*8, dimension(3,3),intent(in)::box
         TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv
         REAL*8, INTENT(IN) :: x
         REAL*8 :: func
       END FUNCTION func
    END INTERFACE
    REAL*8, PARAMETER :: GOLD=1.618034D0,GLIMIT=100.0D0,TINY=1.0D0-20
    REAL*8 :: fu,q,r,u,ulim
    fa=func(ax,chg,mass,box,tot_n_mole,n_mole,n_atom,dfti_desc,dfti_desc_inv,log_file)
    fb=func(bx,chg,mass,box,tot_n_mole,n_mole,n_atom,dfti_desc,dfti_desc_inv,log_file)

!!$    write(*,*) "mnbrak"
!!$    write(*,*) fa, ax
!!$    write(*,*) fb, bx

    if (fb > fa) then
       call swap(ax,bx)
       call swap(fa,fb)
    end if
    cx=bx+GOLD*(bx-ax)

!!$    write(*,*) "cx", cx

    fc=func(cx,chg,mass,box,tot_n_mole,n_mole,n_atom,dfti_desc,dfti_desc_inv,log_file)

    do
       if (fb < fc) RETURN
       r=(bx-ax)*(fb-fc)
       q=(bx-cx)*(fb-fa)
       u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0D0*sign(max(abs(q-r),TINY),q-r))
       ulim=bx+GLIMIT*(cx-bx)
       if ((bx-u)*(u-cx) > 0.0) then
          fu=func(u,chg,mass,box,tot_n_mole,n_mole,n_atom,dfti_desc,dfti_desc_inv,log_file)
          if (fu < fc) then
             ax=bx
             fa=fb
             bx=u
             fb=fu
             RETURN
          else if (fu > fb) then
             cx=u
             fc=fu
             RETURN
          end if
          u=cx+GOLD*(cx-bx)
          fu=func(u,chg,mass,box,tot_n_mole,n_mole,n_atom,dfti_desc,dfti_desc_inv,log_file)
       else if ((cx-u)*(u-ulim) > 0.0) then
          fu=func(u,chg,mass,box,tot_n_mole,n_mole,n_atom,dfti_desc,dfti_desc_inv,log_file)
          if (fu < fc) then
             bx=cx
             cx=u
             u=cx+GOLD*(cx-bx)
             call shft(fb,fc,fu,func(u,chg,mass,box,tot_n_mole,n_mole,n_atom,dfti_desc,dfti_desc_inv,log_file))
          end if
       else if ((u-ulim)*(ulim-cx) >= 0.0) then
          u=ulim
          fu=func(u,chg,mass,box,tot_n_mole,n_mole,n_atom,dfti_desc,dfti_desc_inv,log_file)
       else
          u=cx+GOLD*(cx-bx)
          fu=func(u,chg,mass,box,tot_n_mole,n_mole,n_atom,dfti_desc,dfti_desc_inv,log_file)
       end if
       call shft(ax,bx,cx,u)
       call shft(fa,fb,fc,fu)
    end do
  CONTAINS
    !BL
    SUBROUTINE shft(a,b,c,d)
      REAL*8, INTENT(OUT) :: a
      REAL*8, INTENT(INOUT) :: b,c
      REAL*8, INTENT(IN) :: d
      a=b
      b=c
      c=d
    END SUBROUTINE shft
  END SUBROUTINE mnbrak


  FUNCTION brent(ax,bx,cx,func,tol,xmin,chg,mass,box,tot_n_mole,n_mole,n_atom,dfti_desc,dfti_desc_inv,log_file)
    use MKL_DFTI
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: ax,bx,cx,tol
    REAL*8, INTENT(OUT) :: xmin
    integer,intent(in)::tot_n_mole,n_mole, log_file
    integer,dimension(:),intent(inout)::n_atom
    real*8, dimension(:,:),intent(inout) :: chg,mass
    real*8, dimension(3,3),intent(in)::box
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv   
    REAL*8 :: brent
    INTERFACE
        FUNCTION func(x,chg,mass,box,tot_n_mole,n_mole,n_atom,dfti_desc,dfti_desc_inv,log_file)
         use MKL_DFTI
         IMPLICIT NONE
         integer,intent(in)::tot_n_mole,n_mole, log_file
         integer,dimension(:),intent(inout)::n_atom
         real*8, dimension(:,:),intent(inout) :: chg,mass
         real*8, dimension(3,3),intent(in)::box
         TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv
         REAL*8, INTENT(IN) :: x
         REAL*8 :: func
       END FUNCTION func  
    END INTERFACE
    INTEGER, PARAMETER :: ITMAX=100
    REAL*8, PARAMETER :: CGOLD=0.3819660D0,ZEPS=1.0D-3*epsilon(ax)
    INTEGER :: iter
    REAL*8 :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
    a=min(ax,cx)
    b=max(ax,cx)
    v=bx
    w=v
    x=v
    e=0.0
    fx=func(x,chg,mass,box,tot_n_mole,n_mole,n_atom,dfti_desc,dfti_desc_inv,log_file)
    fv=fx
    fw=fx
    do iter=1,ITMAX
       xm=0.5D0*(a+b)
       tol1=tol*abs(x)+ZEPS
       tol2=2.0D0*tol1
       if (abs(x-xm) <= (tol2-0.5D0*(b-a))) then
          xmin=x
          brent=fx
          RETURN
       end if
       if (abs(e) > tol1) then
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.0D0*(q-r)
          if (q > 0.0) p=-p
          q=abs(q)
          etemp=e
          e=d
          if (abs(p) >= abs(0.5D0*q*etemp) .or. &
               p <= q*(a-x) .or. p >= q*(b-x)) then
             e=merge(a-x,b-x, x >= xm )
             d=CGOLD*e
          else
             d=p/q
             u=x+d
             if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
          end if
       else
          e=merge(a-x,b-x, x >= xm )
          d=CGOLD*e
       end if
       u=merge(x+d,x+sign(tol1,d), abs(d) >= tol1 )
       fu=func(u,chg,mass,box,tot_n_mole,n_mole,n_atom,dfti_desc,dfti_desc_inv,log_file)
       if (fu <= fx) then
          if (u >= x) then
             a=x
          else
             b=x
          end if
          call shft(v,w,x,u)
          call shft(fv,fw,fx,fu)
       else
          if (u < x) then
             a=u
          else
             b=u
          end if
          if (fu <= fw .or. w == x) then
             v=w
             fv=fw
             w=u
             fw=fu
          else if (fu <= fv .or. v == x .or. v == w) then
             v=u
             fv=fu
          end if
       end if
    end do
    call nrerror('brent: exceed maximum iterations')
  CONTAINS
    !BL
    SUBROUTINE shft(a,b,c,d)
      REAL*8, INTENT(OUT) :: a
      REAL*8, INTENT(INOUT) :: b,c
      REAL*8, INTENT(IN) :: d
      a=b
      b=c
      c=d
    END SUBROUTINE shft
  END FUNCTION brent


  SUBROUTINE swap(a,b)
    REAL*8, INTENT(INOUT) :: a,b
    REAL*8 :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap

  SUBROUTINE nrerror(string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    write (*,*) 'nrerror: ',string
    STOP 'program terminated by nrerror'
  END SUBROUTINE nrerror

end module conj_grad




!******************************* kaisq routines ***************
! this is just a wrapper to energy evaluation
!**************************************************************
function kaisq(p,chg,mass,box,tot_n_mole,n_mole,n_atom,dfti_desc,dfti_desc_inv,log_file)
  use MKL_DFTI
  use global_variables
  use routines
  use ms_evb
  real*8  :: kaisq
  real*8,dimension(:), intent(in) :: p
  integer,intent(in)::tot_n_mole,n_mole, log_file
  integer,dimension(:),intent(inout)::n_atom
  real*8, dimension(:,:),intent(inout) :: chg,mass
  real*8, dimension(3,3),intent(in)::box
  TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv

  real*8,dimension(MAX_N_MOLE,MAX_N_ATOM,3)::xyz, force
  real*8,dimension(MAX_N_MOLE,3):: r_com
  integer :: i_mole, i_atom, i, index
  real*8 :: energy

  ! test
  real*8, dimension(MAX_N_MOLE,MAX_N_ATOM,3) :: dr, forcejunk, xyztest
  real*8  :: denergy, energy_new, energy_new_deriv
  ! test

  ! reconstruct xyz array from p
  index=1
  do i_mole=1,n_mole
     do i_atom=1,n_atom(i_mole)
        do i=1,3
           xyz(i_mole,i_atom,i) = p(index)
           index = index+1
        enddo
!!$           write(*,'(A5,3F14.6)') atype_name(atom_index(i_mole,i_atom)), xyz(i_mole,i_atom,:)
     enddo
  enddo

  ! get center of mass for these coordinates
  do i_mole=1,n_mole
     ! after updating atomic coordinates, calculate new center of mass of molecules
     r_com(i_mole,:) = pos_com( xyz, i_mole, n_atom, mass )
     ! translate molecules back into the box
     call shift_move(n_atom(i_mole),xyz(i_mole,:,:),r_com(i_mole,:),box)
  end do

  ! get energy
  call ms_evb_calculate_total_force_energy( force, energy, tot_n_mole, n_mole, n_atom, xyz, r_com, chg, mass, box, dfti_desc,dfti_desc_inv,log_file )

  kaisq = energy

!!$write(*,*) "energy", energy

  ! test
!!$  denergy=0d0
!!$  do i_mole=1,n_mole
!!$     do i_atom=1,n_atom(i_mole)
!!$        do i=1,3
!!$         call random_number( dr(i_mole,i_atom,i) )
!!$        enddo
!!$        dr(i_mole,i_atom,:) = dr(i_mole,i_atom,:) * 0.001d0
!!$        xyztest = xyz
!!$        xyztest(i_mole,i_atom,:) = xyz(i_mole,i_atom,:) + dr(i_mole,i_atom,:)
!!$        denergy = -dot_product(force(i_mole,i_atom,:), dr(i_mole,i_atom,:) )
!!$        energy_new_deriv = energy + denergy
!!$        call ms_evb_calculate_total_force_energy( forcejunk, energy_new, tot_n_mole, n_mole, n_atom, xyztest, r_com, chg, mass, box, dfti_desc,dfti_desc_inv,log_file )
!!$        write(*,*) i_mole, i_atom
!!$        write(*,*) energy, energy_new, energy_new_deriv
!!$     enddo
!!$  enddo
!!$ stop
  ! test

end function kaisq




!***************************************
! this is just a wrapper to force evaluation
!****************************************
function dkaisq(p,chg,mass,box,tot_n_mole,n_mole,n_atom,dfti_desc,dfti_desc_inv,log_file)
  use MKL_DFTI
  use global_variables
  use routines
  use ms_evb
  real*8,dimension(:), intent(in) :: p
  REAL*8, DIMENSION(size(p)) :: dkaisq
  integer,intent(in)::tot_n_mole,n_mole, log_file
  integer,dimension(:),intent(inout)::n_atom
  real*8, dimension(:,:),intent(inout) :: chg,mass
  real*8, dimension(3,3),intent(in)::box
  TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv

  real*8,dimension(MAX_N_MOLE,MAX_N_ATOM,3)::xyz, force
  real*8,dimension(MAX_N_MOLE,3):: r_com
  integer :: i_mole, i_atom, i, index
  real*8 :: energy

  ! reconstruct xyz array from p
  index=1
  do i_mole=1,n_mole
     do i_atom=1,n_atom(i_mole)
        do i=1,3
           xyz(i_mole,i_atom,i) = p(index)
           index = index+1
        enddo
     enddo
  enddo

  ! get center of mass for these coordinates
  do i_mole=1,n_mole
     ! after updating atomic coordinates, calculate new center of mass of molecules
     r_com(i_mole,:) = pos_com( xyz, i_mole, n_atom, mass )
     ! translate molecules back into the box
     call shift_move(n_atom(i_mole),xyz(i_mole,:,:),r_com(i_mole,:),box)
  end do

  ! get force
  call ms_evb_calculate_total_force_energy( force, energy, tot_n_mole, n_mole, n_atom, xyz, r_com, chg, mass, box, dfti_desc,dfti_desc_inv,log_file )

  ! now reorganize force to correspond to p
  index=1
  do i_mole=1,n_mole
     do i_atom=1,n_atom(i_mole)
        do i=1,3
           ! want gradient, not force (minus sign)
           dkaisq(index)=-force(i_mole,i_atom,i)
           index=index+1
        enddo
     enddo
  enddo


end function dkaisq




