
module convolution_par

!   integer ::myrank,nproc
   complex(kind=8), dimension(:,:,:),allocatable ::disp,final_disp
   double precision, dimension(:),allocatable :: station_lat,station_lon
   character(len=80), dimension(:),allocatable ::station_name
   double precision, dimension(:,:), allocatable ::station_coord


   complex(kind=8), dimension(:,:,:),allocatable ::disp_Green_elas
   complex(kind=8), dimension(:,:,:),allocatable ::epsilon_Green
   complex(kind=8), dimension(:,:,:),allocatable ::disp_Green_acous
   complex(kind=8), dimension(:,:),allocatable ::pressure_Green

   integer ::nfreq_Green,ndep_acous,ndep_elas,ntheta,&
             nstat_green_elas,nstat_green_acous
   double precision ::min_dep,ddep,min_theta,dtheta


   double precision, dimension(:,:,:,:),allocatable ::rot_matrix_bound
   double precision, dimension(:,:,:,:),allocatable ::rot_matrix_station
!   double precision, dimension(:,:),allocatable ::gcarc
   integer, dimension(:,:),allocatable ::id_green
   double precision, dimension(:,:),allocatable ::igcarc_remained
       
   

   integer nstation,max_nstation
   double precision ::time_series_length,omega_imag
   real,dimension(:),allocatable ::work_time
   complex(kind=8), dimension(:),allocatable ::work_spc


end module
