subroutine do_convolution(ifreq)
use convolution_par
use coupling_SEM_DSM_par,only:nproc,myrank,npack,fft_disp_bound_new,fft_traction_bound_new,&
                        normal_vector,coord_bound,c,jacobian,NPOINTS,is_elastic,&
                        is_acoustic

use constants
implicit none

!input/output
integer ::ifreq

!other variables
integer ::ista,ipoint
integer ::icomp,icomp1,icomp2,icomp_force
integer ::idGreen
double precision,dimension(ncomp,ncomp,ncomp,ncomp) ::my_c_bound
double precision,dimension(ncomp,ncomp) ::rot1,rot2
complex(kind=8),dimension(ncomp) ::my_disp_bound,my_traction_bound
integer, dimension(ncomp,ncomp) ::icomp_sigma_table
complex(kind=8),dimension(ncomp,ncomp) ::disp_Green_sphe,disp_Green_carte
complex(kind=8),dimension(ncomp,ncomp,ncomp) ::epsilon_Green_sphe
complex(kind=8),dimension(ncomp,ncomp) ::traction_Green_carte
complex(kind=8) ::term1,term2

double precision, dimension(nfit_Green)::xarray_in,imag_yarray_in,real_yarray_in
double precision:: dist_thispoint,imag_y_out,real_y_out,error_esti
integer ::idist_in_table,ifit


!character(len=80) ::fft_file_check
!double precision,dimension(ncomp) ::normal_vector
!double precision,dimension(ncomp,ncomp,ncomp,ncomp) ::elastic_coef


integer ::ipoints1,ipoints2

data icomp_sigma_table  /1,4,5,4,2,6,5,6,3/


   ipoints1=1
   ipoints2=npoints
 do ipoint=ipoints1,ipoints2!1,npoints
   my_disp_bound(:)=fft_disp_bound_new(ifreq,:,ipoint)
   my_traction_bound(:)=fft_traction_bound_new(ifreq,:,ipoint)
!   if(myrank.eq.40.and.ipoint.eq.3) then
!     write(fft_file_check,"('OUTPUT_FILES/freq_',i5.5)")ifreq
!     open(unit=40,file=trim(fft_file_check),action='write',form="formatted",status="unknown")
!     write(40,*)my_disp_bound(1)
!     write(40,*)my_disp_bound(2)
!     write(40,*)my_disp_bound(3)
!     write(40,*)my_traction_bound(1)
!     write(40,*)my_traction_bound(2)
!     write(40,*)my_traction_bound(3)
!     close(40)
!   end if

   do ista=1,nstation

     idGreen=id_Green(ipoint,ista)
     dist_thispoint=distance(ipoint,ista)
     idist_in_table=idistance(ipoint,ista)
!     if(myrank.eq.34) print *,'idGreen_iproc41',idGreen,ista,normal_vector(:,ipoint)
     do ifit=1,nfit_Green
        xarray_in(ifit)=min_theta+((idist_in_table-nfit_Green/2+ifit) -1)*dtheta
     end do

     do icomp1=1,ncomp
       do icomp_force=1,ncomp
        if(is_elastic(ipoint)) then
!         disp_Green_sphe(icomp1,icomp_force)= &
!                        disp_Green_elas(icomp1,idGreen,icomp_force)*(1.d0-ratio_dist_remained) +&
!                        disp_Green_elas(icomp1,idGreen+1,icomp_force)*ratio_dist_remained
           do ifit=1,nfit_Green
!             xarray_in(ifit)=min_theta+(idist_in_table-nfit_Green/2+ifit)*dtheta
             imag_yarray_in(ifit)=aimag(disp_Green_elas(icomp1,idGreen-nfit_Green/2+ifit,icomp_force))
             real_yarray_in(ifit)=real(disp_Green_elas(icomp1,idGreen-nfit_Green/2+ifit,icomp_force))
           end do
         call updown(nfit_Green,xarray_in,imag_yarray_in,dist_thispoint,imag_y_out,error_esti)
         call updown(nfit_Green,xarray_in,real_yarray_in,dist_thispoint,real_y_out,error_esti)
         disp_Green_sphe(icomp1,icomp_force)=cmplx(real_y_out,imag_y_out) 
         if(ista==1.and.ipoint.eq.1.and.ifreq.eq.4096.and.myrank.eq.0) then
            print *,'check_updown',icomp1,icomp_force,xarray_in(1:4),imag_yarray_in(1:4),dist_thispoint,imag_y_out,&
             real_yarray_in(1:4),real_y_out,disp_Green_sphe(icomp1,icomp_force),disp_Green_elas(icomp1,idGreen,icomp_force)
         end if
!         disp_Green_sphe(icomp1,icomp_force)= &
!                        disp_Green_elas(icomp1,idGreen,icomp_force)
 !        disp_Green_sphe(icomp1,icomp_force)=cmplx(0.d0)
           do icomp2=1,ncomp
!             epsilon_Green_sphe(icomp1,icomp2,icomp_force)= &
!                    epsilon_Green(idGreen,icomp_sigma_table(icomp1,icomp2),icomp_force)*(1.d0-ratio_dist_remained) +&
!                    epsilon_Green(idGreen+1,icomp_sigma_table(icomp1,icomp2),icomp_force)*ratio_dist_remained
              do ifit=1,nfit_Green
                imag_yarray_in(ifit)=aimag(epsilon_Green(idGreen-nfit_Green/2+ifit,icomp_sigma_table(icomp1,icomp2),icomp_force))
                real_yarray_in(ifit)=real(epsilon_Green(idGreen-nfit_Green/2+ifit,icomp_sigma_table(icomp1,icomp2),icomp_force))
              end do
              call updown(nfit_Green,xarray_in,imag_yarray_in,dist_thispoint,imag_y_out,error_esti)
              call updown(nfit_Green,xarray_in,real_yarray_in,dist_thispoint,real_y_out,error_esti)
              epsilon_Green_sphe(icomp1,icomp2,icomp_force)=cmplx(real_y_out,imag_y_out)

!             epsilon_Green_sphe(icomp1,icomp2,icomp_force)= &
!                    epsilon_Green(idGreen,icomp_sigma_table(icomp1,icomp2),icomp_force)
         !test
        !if(idGreen<=ntheta)
        !epsilon_Green_sphe(icomp1,icomp2,icomp_force)=cmplx(0.d0)
           end do
        else
!         disp_Green_sphe(icomp1,icomp_force)=disp_Green_acous(icomp1,idGreen,icomp_force)*(1.d0-ratio_dist_remained) +&
!               disp_Green_acous(icomp1,idGreen+1,icomp_force)*ratio_dist_remained
!to be fixed
         disp_Green_sphe(icomp1,icomp_force)= cmplx(0.d0)
           do icomp2=1,ncomp
             if(icomp2.eq.icomp1) then
!                 epsilon_Green_sphe(icomp1,icomp2,icomp_force)= &
!                    pressure_Green(idGreen,icomp_force)*(1.d0-ratio_dist_remained) &
!                    +pressure_Green(idGreen+1,icomp_force)*ratio_dist_remained

!to be fixed
                  epsilon_Green_sphe(icomp1,icomp2,icomp_force)= cmplx(0.d0)
             else
                  epsilon_Green_sphe(icomp1,icomp2,icomp_force)= cmplx(0.d0)
             end if
           end do
        end if
       end do
     end do

     rot1(:,:)=rot_matrix_bound(:,:,ipoint,ista)
     rot2(:,:)=rot_matrix_station(:,:,ipoint,ista)
     disp_Green_carte=matmul(transpose(rot1),matmul(disp_Green_sphe,rot2))
!km to meter
!     disp_Green_carte=disp_Green_carte*km
!     traction_Green_carte(:,:)=cmplx(0.d0)
       traction_Green_carte(:,:)=cmplx(0.d0)
       call copy_matrix_c(c(1,ipoint),my_c_bound(1,1,1,1))

       call calcu_traction_Green_carte1(traction_Green_carte,rot1,rot2, &
                normal_vector(1,ipoint), &
                epsilon_Green_sphe(1,1,1),ncomp)

    do icomp=1,ncomp
       do icomp1=1,ncomp
          term1=disp_Green_carte(icomp1,icomp)*my_traction_bound(icomp1)
          term2=-my_disp_bound(icomp1)*traction_Green_carte(icomp1,icomp)
!          term2=cmplx(0.d0,0.d0)
!z component only depneds on GREEN functions for vertical single force
          disp(icomp,ifreq,ista)=disp(icomp,ifreq,ista)+(term1+term2)*jacobian(ipoint)
       end do
     end do

    end do
 end do 

end subroutine do_convolution



subroutine copy_matrix_c(c,my_c_bound)

use constants
implicit none

double precision, dimension(N_ELAS_COEF),intent(in) ::c
double precision,dimension(ncomp,ncomp,ncomp,ncomp),intent(out) ::my_c_bound
!other variables
integer ::i,j,k,l
integer ::index1,index2
!integer ::index_temp
integer, dimension(ncomp,ncomp) ::index_table
integer, dimension(6,6) ::table_c66
data index_table  /1,6,5,6,2,4,5,4,3/
data table_c66 /1,2,3,4,5,6,&
                2,7,8,9,10,11,&
                3,8,12,13,14,15,&
                4,9,13,16,17,18,&
                5,10,14,17,19,20,&
                6,11,15,18,20,21/
!       [1 5 6]        [11 12 13]
!       [5 2 4] =      [21 22 23]
!       [6 4 3]        [31 32 33]

do i=1,ncomp
  do j=1,ncomp
    index1=index_table(i,j)
    do k=1,ncomp
      do l=1,ncomp
        index2=index_table(k,l)
!        if(index1.gt.index2) then
!             index_temp=index1
!             index1=index2
!             index2=index_temp
!        end if
        my_c_bound(i,j,k,l)=c(table_c66(index1,index2))
      end do
    end do
  end do
end do


end subroutine





subroutine calcu_traction_Green_carte(traction_Green_carte,rot1,rot2, &
               c,normal_vector,epsilon_Green_sphe,ncomp)
   implicit none

!input/output
   integer ::ncomp
   complex(kind=8),dimension(ncomp,ncomp),intent(out) ::traction_Green_carte
   double precision,dimension(ncomp,ncomp),intent(in) ::rot1,rot2
   complex(kind=8),dimension(ncomp,ncomp,ncomp),intent(in) ::epsilon_Green_sphe
   double precision,dimension(ncomp,ncomp,ncomp,ncomp),intent(in) ::c
   double precision, dimension(ncomp),intent(in) ::normal_vector
!other variables
   integer ::i,j,k,l,n
   complex(kind=8),dimension(ncomp,ncomp) ::temp_sphe,temp_carte
   complex(kind=8),dimension(ncomp,ncomp,ncomp) ::epsilon_Green_carte
   complex(kind=8),dimension(ncomp,ncomp) ::traction_temp
   double precision ::devia_factor


   do n=1,ncomp
    temp_sphe(:,:)=epsilon_Green_sphe(:,:,n)
    temp_carte=matmul(transpose(rot1),matmul(temp_sphe,rot1))
    epsilon_Green_carte(:,:,n)=temp_carte(:,:)
   end do
   traction_temp(:,:)=cmplx(0.d0)
   do i=1,ncomp
    do n=1,ncomp
     do j=1,ncomp
      do k=1,ncomp
       do l=1,ncomp
         if(k.ne.l) then
            devia_factor=1.0
         else
            devia_factor=1.0
         end if
         traction_temp(i,n)=traction_temp(i,n)+c(i,j,k,l)*normal_vector(j)*&
                            epsilon_Green_carte(k,l,n)*devia_factor
       end do
      end do
     end do
    end do
   end do
 
   traction_Green_carte=matmul(traction_temp,rot2)

end subroutine

subroutine calcu_traction_Green_carte1(traction_Green_carte,rot1,rot2, &
               normal_vector,epsilon_Green_sphe,ncomp)
   implicit none

!input/output
   integer ::ncomp
   complex(kind=8),dimension(ncomp,ncomp),intent(out) ::traction_Green_carte
   double precision,dimension(ncomp,ncomp),intent(in) ::rot1,rot2
   complex(kind=8),dimension(ncomp,ncomp,ncomp),intent(in) ::epsilon_Green_sphe
   double precision, dimension(ncomp),intent(in) ::normal_vector
!other variables
   integer ::i,j,n
   complex(kind=8),dimension(ncomp,ncomp) ::temp_sphe,temp_carte
   complex(kind=8),dimension(ncomp,ncomp,ncomp) ::epsilon_Green_carte
   complex(kind=8),dimension(ncomp,ncomp) ::traction_temp


   do n=1,ncomp
    temp_sphe(:,:)=epsilon_Green_sphe(:,:,n)
    temp_carte=matmul(transpose(rot1),matmul(temp_sphe,rot1))
    epsilon_Green_carte(:,:,n)=temp_carte(:,:)
   end do
   traction_temp(:,:)=cmplx(0.d0)
   do i=1,ncomp
    do n=1,ncomp
      do j=1,ncomp
         traction_temp(i,n)=traction_temp(i,n)+normal_vector(j)*&
                            epsilon_Green_carte(i,j,n)
      end do
    end do
   end do

   traction_Green_carte=matmul(traction_temp,rot2)

end subroutine

