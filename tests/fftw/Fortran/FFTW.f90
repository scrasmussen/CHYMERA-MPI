Module FFTW
!=============================================================================
! Name        : FFTW.f90
! Author(s)   : Craig Rasmussen
! Version     :
! Copyright   : University of Oregon, 2014
! Description : Module wrapper for the FFTW Fortran interfaces
!
! Method :
!
!   Provides a simple module wrapper of the FFTW Fortran language
!   interoperability interfaces provided by the FFTW distribution.
!   See http://www.fftw.org/fftw3.pdf for the C documentation and
!   the file, include/fftw3.f03 in the FFTW distribution.
!
!=============================================================================
Use, Intrinsic :: ISO_C_BINDING
Include "fftw3.f03"

End Module FFTW
