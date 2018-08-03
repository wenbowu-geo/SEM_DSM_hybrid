subroutine partition_job()
use coupling_SEM_DSM_par
use constants

implicit none

#ifdef USE_MPI
  include "mpif.h"
#endif

#ifdef USE_MPI
   integer npack_remained
#endif

   if(npackage_SEM.lt.nproc) stop 'Error, npack should be larger than nproc'

#ifdef USE_MPI
   npack=int(npackage_SEM/nproc)
   npack_remained=mod(npackage_SEM,nproc)
   ipack_start=npack*myrank+1
   ipack_end=ipack_start+npack-1


   if(myrank.lt.npack_remained) then
      npack=npack+1
      ipack_start=ipack_start+myrank
      ipack_end=ipack_end+myrank+1
   else
      ipack_start=ipack_start+npack_remained
      ipack_end=ipack_end+npack_remained
   end if
#else 
   npack=npackage_SEM
   ipack_start=1
   ipack_end=npackage_SEM
#endif


   if(ipack_end-ipack_start+1.ne.npack) stop 'Error in counting packages'
   if(npack.gt.max_npackages) then
      stop 'the size of array is larger than max_array'
   end if



end subroutine partition_job
