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

type(C_PTR) :: pF, pFn, pG, pGn, plan_forward, plan_backward, plan_g
real(kind_r), pointer :: F(:,:), G(:)
complex(kind_c), pointer :: Fn(:,:), Gn(:)

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
                   fftw_plan_many_dft_c2r, fftw_plan_dft_r2c_1d, FFTW_MEASURE
  implicit none
  real(kind_r), intent(OUT), pointer :: U(:,:)

  integer :: i, j, mode
  integer :: rank, howmany, idist, odist, istride, ostride
  integer(C_SIZE_T), parameter :: nModes = N

!-----------------------------------------------------------------------------

  ! calculate plan for FFTW
  !

  print *, "nModes=", nModes
  print *, "PI=", PI

  mode = 2

  pF  = fftw_alloc_real(nModes*(M+2))
  pFn = fftw_alloc_real((1+nModes/2)*(M+2))
  pG  = fftw_alloc_real(nModes)
  pGn = fftw_alloc_real(1+nModes/2)

  call C_F_Pointer(pF,  F,  shape=[N,M+2])        ! includes halo region
  call C_F_Pointer(pFn, Fn, shape=[1+N/2,M+2])    ! includes halo region
  call C_F_Pointer(pG,  G,  shape=[N])
  call C_F_Pointer(pGn, Gn, shape=[1+N/2])

  rank    = 1     ! 1D transforms
  howmany = M     ! one transform for each interior radial grid point
  idist   = N     ! distance between first element of first array and first element of second array
  odist   = N
  istride = 1     ! distance between elements in the transform dimension
  ostride = 1

!TODO - DELETEME
!  plan_forward  = fftw_plan_dft_r2c_1d(N, F(:,mode), Fn, FFTW_MEASURE)
!  plan_backward = fftw_plan_dft_c2r_1d(N, Fn(:,mode), F(:,mode), FFTW_MEASURE)

  plan_forward  = fftw_plan_many_dft_r2c(rank,[N],howmany,F ,[N],istride,idist,               &
                                                          Fn,[N],ostride,odist,FFTW_MEASURE)
  plan_backward = fftw_plan_many_dft_c2r(rank,[N],howmany,Fn,[N],istride,idist,               &
                                                          F ,[N],ostride,odist,FFTW_MEASURE)
  plan_g        = fftw_plan_dft_r2c_1d(N, G, Gn, FFTW_MEASURE)

  !... Allocate independent and dependent variables
  !    --------------------------------------------
  allocate(Th(0:N-1), R(0:M+1))
  allocate(U(0:N-1,0:M+1), F(0:N-1,0:M+1))

  !... Initialize arrays and set boundary conditions
  !    ---------------------------------------------

  do i = 0, M+1
     R(i) = (i - 0.5_kind_r)*DR
  end do

  do j = 0, N-1
     Th(j) = j*DTH
  end do

  do i = 1, M+1
     F(:,i) = 3*cos(Th)
  end do

  ! TODO - check sin(Th), doesn't seem to work
  G(:) = 5 + cos(Th)     ! outer boundary

  U = ZERO
  U(:,M+1) = G(:)

End Subroutine Initialize


Subroutine Finalize()
!=============================================================================
! Free memory and FFTW plans
!=============================================================================
  use FFTW, only : fftw_destroy_plan, fftw_free
  implicit none

!-----------------------------------------------------------------------------

  call fftw_destroy_plan(plan_forward)
  call fftw_destroy_plan(plan_backward)
  call fftw_destroy_plan(plan_g)

  call fftw_free(pF)
  call fftw_free(pFn)
  call fftw_free(pG)
  call fftw_free(pGn)

End Subroutine Finalize


Subroutine Transform()
!=============================================================================
! Transform the mass density forcing function
!=============================================================================
  use FFTW, only : fftw_execute_dft_r2c, fftw_execute_dft_c2r
  implicit none
  integer :: mode

!-----------------------------------------------------------------------------

  mode = 2

  print *
  print *, "Initial:"
  print*, F(:,1)
  print*, F(:,2)
  print *

  call fftw_execute_dft_r2c(plan_forward, F, Fn)
  Fn = Fn/N

  print *, "Forward(Fn):"
  print *, Fn(:,1)
  print *, Fn(:,2)
  print *

  call fftw_execute_dft_r2c(plan_g, G, Gn)
  Gn = Gn/N

  print *, "Forward(Gn):"
  print *, Gn
  print *

  call fftw_execute_dft_c2r(plan_backward, Fn, F)
  print *, "Backward:"
  print *, F(:,1)
  print *, F(:,2)
  print *

stop

End Subroutine Transform


Subroutine Solve()
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

  mode = 2

! TODO - make matrix coefficients correctly complex 
!
  do i = 1, M
     a = 1.0d0 / (2*i - 1)

     DL(i) = 1 - a                        ! first value not used
     DU(i) = 1 + a                        ! last  value not used
     D (i) = -2*(1 + 2*(mode*a)**2)

!    B(i,1) = Fn(mode,i)*DR*DR   !TODO-FIXME
     B(i,1) = Fn(mode,i)*DR*DR

     print *, "MATrix:", i, a, DL(i), DU(i), D(i), B(i,1)

  end do

! TODO - need all modes because G has constant
  B(M,1) = B(M,1) - DU(M)*Gn(mode)

  call zgtsv(M, NRHS, DL(2), D, DU, B, LDB, info)

  print *, "info=", info
  print *, B(:,1)

End Subroutine Solve


End Module Lai2D_1
