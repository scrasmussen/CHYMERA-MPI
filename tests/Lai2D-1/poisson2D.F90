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
  call Solve

  call Finalize

#ifdef NOT_YET
  call error_allocate(M, N, H)

  !call VTK_Output (M, N, H, T1, 0, 'T', aContext)
  !call VTK_Output (M, N, H, T2, 1, 'T', aContext)

  call VTK_Output (M, N, H, Err, 0, 'E', aContext)

  Print*, "I am Alive"

  !... Iterate solution
  !    ________________
  do i = 2, nsteps, 2
     call Iterate (M, N, H, T1, T2);
     call Parallel_HaloArray_Exchange (ha2, aTopology)
     call Iterate (M, N, H, T2, T1);
     call Parallel_HaloArray_Exchange (ha1, aTopology)
     call calc_error (m, n, h, T2)
     call VTK_Output (M, N, 0, Err, i  , 'E', aContext)
     !call VTK_Output (M, N, H, T2, i  , 'T', aContext)
     !call VTK_Output (M, N, H, T1, i+1, 'T', aContext)
  end do

  call Parallel_HaloArray_Deallocate ( ha1, T1 )
  call Parallel_HaloArray_Deallocate ( ha2, T2 )

  call Parallel_End (aContext)
#endif

End Program Poisson2D
