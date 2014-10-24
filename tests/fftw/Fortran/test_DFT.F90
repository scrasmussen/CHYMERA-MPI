program Test_DFT
  use params,  only : M, N, H, R, TH, PI, kind_r
  use DFT, only : Initialize, Finalize, Transform, TransformFFTW

  implicit none

!-----------------------------------------------------------------------------

  !... Initialize arrays and set boundary conditions
  !    ---------------------------------------------

  call Initialize

  print '(16(g10.3,1x))', "theta:", Th/(2*PI)
  print '(16(g10.3,1x))', "r:", R

  call Transform
  call TransformFFTW

  call Finalize


end program Test_DFT