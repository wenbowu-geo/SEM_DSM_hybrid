subroutine write_fft()
   use convolution_par
   use constants
   implicit none

!other variables
   character(len=80) ::fft_file_disp,fft_file_disp_icomp
   integer ::istat,icomp,ifreq
   integer ::ios
   character(len=4),dimension(ncomp) ::comp
   data comp /'.bhx','.bhy','.bhz'/

   do istat=1,nstation
     fft_file_disp="./OUTPUT_FILES/"//trim(station_name(istat))
     do icomp=1,ncomp
        fft_file_disp_icomp=trim(fft_file_disp)//comp(icomp)
        open(unit=40,file=trim(fft_file_disp_icomp),action='write',form="formatted",status="unknown",iostat=ios)
        do ifreq=1,nfreq_Green
          write(40,*)final_disp(icomp,ifreq,nstation)
        end do
     end do
   end do

end subroutine write_fft
