Module Params
!=============================================================================
! Name        : Params.f90
! Author(s)   : Craig Rasmussen
! Version     :
! Copyright   : University of Oregon, 2014
! Description : Parameters for solving Poisson's equation in two dimensions (phi,r)
!
!=============================================================================
  use, intrinsic :: ISO_C_BINDING, only : C_DOUBLE, C_DOUBLE_COMPLEX
  implicit none

  integer, parameter :: NDIMS  =   2
  integer, parameter :: M      =   2
  integer, parameter :: N      =   6
  integer, parameter :: H      =   1

  integer, parameter :: kind_r = C_DOUBLE
  integer, parameter :: kind_c = C_DOUBLE_COMPLEX

  real(kind_r), parameter :: ZERO = 0.0_kind_r

  real(kind_r), parameter :: PI = 2.0_kind_r*asin(1.0_kind_r)

  real(kind_r), parameter :: DR  = 2.0_kind_r/(2.0_kind_r*M + 1.0_kind_r)
  real(kind_r), parameter :: DTH = 2.0_kind_r*PI/N

  !... dependent variables
  !    -------------------
  real(kind_r), allocatable :: Th(:), R(:)

End Module Params
