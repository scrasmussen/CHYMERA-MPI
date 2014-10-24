Program Poisson2D
!=============================================================================
! Name        : poisson2D.F90
! Author(s)   : Craig Rasmussen
! Version     :
! Copyright   : University of Oregon, 2014
! Description : Solves the Poisson equation in two dimensions (phi,r)
!
! Method :
!
!   Uses method of Lai, "A simple compact fourth-order Poisson solver on
!   "polar geometry", JCP, 182, 337-345, 2002.
!
!=============================================================================
  use params,  only : M, N, H, R, TH, kind_r
  use Lai2D_1, only : Initialize, Finalize, Transform, Solve

  implicit none

  real(kind_r), pointer :: U(:,:)

  integer :: i, j

!-----------------------------------------------------------------------------

  !... Initialize arrays and set boundary conditions
  !    ---------------------------------------------

  call Initialize(U)

  print *, "theta:", Th
  print *, "r:", R

  call Transform

  call Finalize


End Program Poisson2D
