subroutine fft(ipoint)
  use coupling_SEM_DSM_par
  use convolution_par, only:nfreq_Green
  use constants
  implicit none

  integer ipoint,icomp,it,ifreq

!other variables

  do it=1,npts_fine
     do icomp=1,ncomp
         fft_disp_bound_fine(it,icomp)=disp_bound_fine(it,icomp)
     end do
     do icomp=1,ncomp
         fft_traction_bound_fine(it,icomp)=traction_bound_fine(it,icomp)
     end do
  end do

  do icomp=1,ncomp
     if(ipoint.gt.npoints) stop 'error: ipoint>npoints'
     call four1(fft_disp_bound_fine(1,icomp),npts_fine,-1)
     do ifreq=1,nfreq_Green
        !deltat_refine- normalization factor
        fft_disp_bound_new(ifreq,icomp,ipoint)=fft_disp_bound_fine(ifreq,icomp)*deltat_refine  
     end do
!     print *,'fft myrank disp_bound',icomp,myrank
  end do
  do icomp=1,ncomp
!     print *,'fft myrank traction_bound',icomp,myrank
     call four1(fft_traction_bound_fine(1,icomp),npts_fine,-1)
     do ifreq=1,nfreq_Green
       !1/deltat_refine- normalization factor
       fft_traction_bound_new(ifreq,icomp,ipoint)=fft_traction_bound_fine(ifreq,icomp)*deltat_refine
     end do
  end do

end subroutine fft


