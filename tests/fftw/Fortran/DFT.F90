!! Using memory alloc'd by FFTW should be faster
!
#define USE_FFTW_ALLOC


module DFT
!=============================================================================

use, intrinsic :: ISO_C_BINDING, only : C_PTR
use params, only : M, N, kind_r, kind_c
implicit none

type(C_PTR) :: pF, pFn, plan_forward, plan_backward
real(kind_r), pointer :: F(:,:)
complex(kind_c), pointer :: Fn(:,:)

contains
!-----------------------------------------------------------------------------


subroutine Initialize
!=============================================================================
! Initialize the dependent variables
!=============================================================================
  use, intrinsic :: ISO_C_BINDING, only : C_SIZE_T, C_F_Pointer
  use Params, only : M, N, DR, DTH, R, TH, PI, ZERO
  use FFTW, only : fftw_alloc_real, fftw_alloc_complex, FFTW_ESTIMATE, &
                   fftw_plan_many_dft_r2c, fftw_plan_many_dft_c2r
  implicit none

  integer :: i, j
  integer :: rank, howmany, idist, odist, istride, ostride
  integer(C_SIZE_T), parameter :: nModes = N

!-----------------------------------------------------------------------------

  ! calculate plan for FFTW
  !

  print *, "nModes=", nModes, "N=", N, "M=", M
  print *, "PI=", PI

#ifdef USE_FFTW_ALLOC
  pF  = fftw_alloc_real(nModes*(M+2))
  pFn = fftw_alloc_complex((1+nModes/2)*(M+2))
  call C_F_Pointer(pF,  F,  shape=[N,M+2])
  call C_F_Pointer(pFn, Fn, shape=[1+N/2,M+2])
#else
  allocate(F(N,M+2), Fn(1+N/2,M+2))
#endif

  rank    = 1     ! 1D transforms
  howmany = M+2   ! one transform for each interior radial grid point
  istride = 1     ! distance between elements in the transform dimension
  ostride = 1

  idist   = N     ! distance between first element of first array and first element of second array
  odist   = N/2+1
  plan_forward  = fftw_plan_many_dft_r2c(rank,[N],howmany,F ,[N],istride,idist,               &
                                                          Fn,[N],ostride,odist,FFTW_ESTIMATE)
  idist   = N/2+1
  odist   = N
  plan_backward = fftw_plan_many_dft_c2r(rank,[N],howmany,Fn,[N],istride,idist,               &
                                                          F ,[N],ostride,odist,FFTW_ESTIMATE)

  !... Allocate independent and dependent variables
  !    --------------------------------------------
  allocate(Th(0:N-1), R(0:M+1))

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
!    F(:,i) = 1 + cos(1*Th) + cos(2*Th) + cos(3*Th)
  end do

  F(:,  1) = F(:,2)          ! inner boundary (for now)
  F(:,M+2) = 5 + cos(Th)     ! outer boundary

  print *, "Initial(F): shape=", shape(F)


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

  nullify(F, Fn)

#ifdef USE_FFTW_ALLOC
  call fftw_free(pFn)
  call fftw_free(pF)
#endif

End Subroutine Finalize


Subroutine TransformFFTW
!=============================================================================
! Transform the mass density forcing function
!=============================================================================
  use FFTW, only : fftw_execute_dft_r2c, fftw_execute_dft_c2r
  implicit none

!-----------------------------------------------------------------------------

  Fn = (0,0)

  print *
  print *, "Initial(F): shape=", shape(F)
  print '(16(f7.3,1x))', F(:,1)
  print '(16(f7.3,1x))', F(:,2)
  print '(16(f7.3,1x))', F(:,3)
  print '(16(f7.3,1x))', F(:,4)
  print *

  call fftw_execute_dft_r2c(plan_forward, F, Fn)

  print *, "Forward(Fn): shape=", shape(Fn)
  print '(16(f7.3,1x))', Fn(:,1)
  print '(16(f7.3,1x))', Fn(:,2)
  print '(16(f7.3,1x))', Fn(:,3)
  print '(16(f7.3,1x))', Fn(:,4)
  print *

  call fftw_execute_dft_c2r(plan_backward, Fn, F)
  F = F/N                    ! normalize so that reverse transform is same as original input
  print *, "Backward(F):"
  print '(16(f7.3,1x))', F(:,1)
  print '(16(f7.3,1x))', F(:,2)
  print '(16(f7.3,1x))', F(:,3)
  print '(16(f7.3,1x))', F(:,4)
  print *

  print *
  print *, "Forward(REAL(Fn)): shape=", shape(Fn)
  print '(16(f7.3,1x))', REAL(Fn(:,1))
  print '(16(f7.3,1x))', REAL(Fn(:,2))
  print '(16(f7.3,1x))', REAL(Fn(:,3))
  print '(16(f7.3,1x))', REAL(Fn(:,4))
  print *

  print *
  print *, "Forward(IMAG(Fn)): shape=", shape(Fn)
  print '(16(f7.3,1x))', AIMAG(Fn(:,1))
  print '(16(f7.3,1x))', AIMAG(Fn(:,2))
  print '(16(f7.3,1x))', AIMAG(Fn(:,3))
  print '(16(f7.3,1x))', AIMAG(Fn(:,4))
  print *

STOP

End Subroutine TransformFFTW


Subroutine Transform
!=============================================================================
! Transform using slow direct DFT transform
!=============================================================================
  use params, only : Th
  implicit none

  integer :: i, j, k
  real(kind_r) :: re, im

!-----------------------------------------------------------------------------

  print *
  print *, "Initial(F): shape=", shape(F)
  print '(16(g10.3,1x))', F(:,1)
  print '(16(g10.3,1x))', F(:,2)
  print '(16(g10.3,1x))', F(:,3)
  print '(16(g10.3,1x))', F(:,4)
  print *

  Fn = (0,0)

  do i = 1, M+2
     do j = 1, N
        k = j - 1
        re = sum(F(:,i)*cos(k*Th))
        im = sum(F(:,i)*sin(k*Th))
        if (abs(re) < 1.0e-14) re = 0
        if (abs(im) < 1.0e-14) im = 0
        Fn(j,i) = cmplx(re,im)
     end do
  end do

#ifdef NO_LONGER
  k = 0
  prod = F(:,1)*cos(k*Th)
  s = sum(prod)
  if (abs(s) < 1.0e-14) s = 0

  Fn(1,1) = s

  print *
  print '(" sum: ", 16(g10.3,1x))', Fn(1,1)
  print *

  k = 1
  prod = F(:,1)*cos(k*Th)
  s = sum(prod)
  if (abs(s) < 1.0e-14) s = 0

  Fn(2,1) = s

  print *
  print '(" sum: ", 16(g10.3,1x))', Fn(2,1)
  print *
#endif

!  call fftw_execute_dft_r2c(plan_forward, F, Fn)
  Fn(:,1) = Fn(:,1)
  Fn(:,2) = Fn(:,2)
  Fn(:,3) = Fn(:,3)
  Fn(:,4) = Fn(:,4)
  Fn = Fn

  print *, "Forward(Fn): shape=", shape(Fn)
  print '(16(g10.3,1x))', Fn(:,1)
  print '(16(g10.3,1x))', Fn(:,2)
  print '(16(g10.3,1x))', Fn(:,3)
  print '(16(g10.3,1x))', Fn(:,4)
  print *

end subroutine Transform


end module DFT
