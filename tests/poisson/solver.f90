       program solver
! This program is an excised version of the poisson solver in CHYMERA
!       implicit none
! initialize everything.
 
#include "hydroparam.h"
#include "globals.h"
#include "units.h"


       COMMON /ITS/ITSTRT,ITSTOP,ITSTEP

! Call the setup routine, get density from fort.2 or fort.7.

       call setup(itstrt,itstop,idiag,isoadi,istor,itstep,isym,maxtrm)
       print *, "called setup"
! Open recording files

! Call the poisson solver
         redge=r(jmax1)
         print *, "set outer boundary"
         call bdygen(maxtrm,isym,redge)
         print *, "called subroutine bdygen()"
         call pot3(8,iprint)
         print *, "called subroutine pot3()"


! Output stuff like Fourier coefficients and potential.

!         potfile='potential.'//tim
!         open(unit=38,file=potfile,form='unformatted')
!         write(38) phi
!         write(38) time



       end
