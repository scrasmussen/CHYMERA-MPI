!! Using memory alloc'd by FFTW should be faster
!
#undef USE_FFTW_ALLOC


Module Lai2D_1
!=============================================================================
! Name        : Lai2D_1.F90
! Author(s)   : Craig Rasmussen
! Version     :
! Copyright   : University of Oregon, 2014
! Description : Solves the Poisson equation in two dimensions (th,r)
!
! Method :
!
!   Uses method of Lai, "A simple compact fourth-order Poisson solver on
!   "polar geometry", JCP, 182, 337-345, 2002.
!
!=============================================================================
use, intrinsic :: ISO_C_BINDING, only : C_PTR
use params, only : M, N, kind_r, kind_c
implicit none

type(C_PTR) :: pF, pFn, plan_forward, plan_backward
real(kind_r), pointer :: F(:,:)
complex(kind_c), pointer :: Fn(:,:)

Contains
!-----------------------------------------------------------------------------


Subroutine Initialize(U)
!=============================================================================
! Initialize the dependent variables
!
! Problem :
!
!   From Lai (2002) test function 1 :
!
!      u(x,y) = 3 e^{x+y} (x - x^2) (y - y^2) + 5
!
!   initially let's try: u(th,r) = r^2 cos(th) + 5
!
!   f(th,r) = 3 cos(th)
!
!=============================================================================
  use, intrinsic :: ISO_C_BINDING, only : C_SIZE_T, C_F_Pointer
  use Params, only : M, N, DR, DTH, R, TH, PI, ZERO
  use FFTW, only : fftw_alloc_real, fftw_plan_many_dft_r2c, &
                   fftw_plan_many_dft_c2r, FFTW_MEASURE
  implicit none
  real(kind_r), intent(OUT), pointer :: U(:,:)

  integer :: i, j
  integer :: rank, howmany, idist, odist, istride, ostride
!  integer(C_SIZE_T), parameter :: nModes = N + 2    ! padding fixes memory error?
  integer(C_SIZE_T), parameter :: nModes = N

!-----------------------------------------------------------------------------

  ! calculate plan for FFTW
  !

  print *, "nModes=", nModes, "N=", N, "M=", M
  print *, "PI=", PI

  pF  = fftw_alloc_real(nModes*(M+2))
  pFn = fftw_alloc_real((1+nModes/2)*(M+2))

#ifdef USE_FFTW_ALLOC
  call C_F_Pointer(pF,  F,  shape=[N,M+2])
  call C_F_Pointer(pFn, Fn, shape=[1+N/2,M+2])
#else
  allocate(F(N,M+2), Fn(1+N/2,M+2))
#endif

  rank    = 1     ! 1D transforms
  howmany = M+2   ! one transform for each interior radial grid point
  idist   = N     ! distance between first element of first array and first element of second array
  odist   = N
  istride = 1     ! distance between elements in the transform dimension
  ostride = 1

  plan_forward  = fftw_plan_many_dft_r2c(rank,[N],howmany,F ,[N],istride,idist,               &
                                                          Fn,[N],ostride,odist,FFTW_MEASURE)
  plan_backward = fftw_plan_many_dft_c2r(rank,[N],howmany,Fn,[N],istride,idist,               &
                                                          F ,[N],ostride,odist,FFTW_MEASURE)

  !... Allocate independent and dependent variables
  !    --------------------------------------------
  allocate(Th(0:N-1), R(0:M+1))
  allocate(U(0:N-1,M+2))

  !... Initialize arrays and set boundary conditions
  !    ---------------------------------------------

  do i = 0, M+1
     R(i) = (i - 0.5_kind_r)*DR
  end do

  do j = 0, N-1
     Th(j) = j*DTH
  end do

  do i = 2, M+1
     F(:,i) = 3*cos(Th)
  end do

  F(:,  1) = F(:,2)          ! inner boundary (for now)
  F(:,M+2) = 5 + cos(Th)     ! outer boundary

  ! TODO - check sin(Th), doesn't seem to work

  U = ZERO
  U(:,M+2) = F(:,M+2)

End Subroutine Initialize


Subroutine Finalize
!=============================================================================
! Free memory and FFTW plans
!=============================================================================
  use FFTW, only : fftw_destroy_plan, fftw_free
  implicit none

!-----------------------------------------------------------------------------

  call fftw_destroy_plan(plan_forward)
  call fftw_destroy_plan(plan_backward)

#ifdef USE_FFTW_ALLOC
  call fftw_free(pF)
  call fftw_free(pFn)
#endif

  nullify(F, Fn)

End Subroutine Finalize


Subroutine Transform
!=============================================================================
! Transform the mass density forcing function
!=============================================================================
  use FFTW, only : fftw_execute_dft_r2c, fftw_execute_dft_c2r
  implicit none

!-----------------------------------------------------------------------------

  print *
  print *, "Initial(F): shape=", shape(F)
  print *, F(:,1)
  print *, F(:,2)
  print *, F(:,3)
  print *, F(:,4)
  print *

  call fftw_execute_dft_r2c(plan_forward, F, Fn)
  Fn = Fn/N

  print *, "Forward(Fn): shape=", shape(Fn)
  print *, REAL(Fn(:,1))
  print *, Fn(:,2)
  print *, Fn(:,3)
  print *, Fn(:,4)
  print *

  call fftw_execute_dft_c2r(plan_backward, Fn, F)
  print *, "Backward(F):"
  print *, F(:,1)
  print *, F(:,2)
  print *, F(:,3)
  print *, F(:,4)
  print *

End Subroutine Transform


Subroutine Solve
!=============================================================================
! Solves the tridiagonal system
!=============================================================================
  use Params, only : M, N, DR, kind_r, kind_c
  implicit none

  integer, parameter :: NRHS = 1    ! the number of right-hand-sides
  integer, parameter :: LDB  = M    ! the leading dimension RHS array
  integer :: info                   ! returned status

  integer :: i, mode
  real(kind_r) :: a
  complex(kind_c) :: DL(M), D(M), DU(M), B(M,NRHS)

!-----------------------------------------------------------------------------

! TODO - make matrix coefficients correctly complex 
!

print *
print *, Fn(1,2:3)
print *
print *, Fn(2,2:3)
print *


  do mode = 0, N-1

     do i = 1, M
        a = 1.0d0/(2*i - 1)

        DL(i) = 1 - a
        DU(i) = 1 + a
        D (i) = -2*(1 + 2*(mode*a)**2)

        !    B(i,1) = Fn(mode+1,i+1)*DR*DR   !TODO-FIXME
        B(i,1) = Fn(mode+1,i+1)

        print *, "MATrix:", mode, i, a
        print *, "    DL:", DL(i)
        print *, "    DU:", DU(i)
        print *, "    D :", D (i)
        print *, "    B :", B (i,1)
        print *

        B(i,1) = B(i,1)*DR*DR

     end do

     !... boundary conditions
     !    -------------------

     D(1) = D(1) + ((-1)**mode) * DL(1)        ! lower boundary goes to diagonal  (TODO-MPI:only r==0)
     B(M,1) = B(M,1) - DU(M)*Fn(mode+1,M+2)    ! upper boundary goes to RHS

     print *, "    BC:", Fn(mode+1,M+2)
     print *, "    BC:", B(M,1)
     print *

     call zgtsv(M, NRHS, DL(2), D, DU, B, LDB, info)

     print *, "mode=", mode, "  info=", info
     print *, B(:,1)
     print *

  end do

End Subroutine Solve


End Module Lai2D_1
