integer, parameter ::max_nfit=131
integer, parameter ::min_nfit_SEM_resample = 8
integer, parameter ::max_npackages=4000
integer, parameter ::TINY=1.e-9
integer, parameter ::N_ELAS_COEF=21
integer, parameter ::ncomp=3
!Be careful to change max_array_size. It depends on your max-memory available.
integer, parameter ::max_array_size=700000000
integer, parameter ::ncomp_stress=6
double precision,parameter ::pi=3.1415926
double precision,parameter ::degtorad=pi/180.d0
double precision,parameter ::r_earth=6371.0*1.0e3
double precision,parameter ::radtodeg=180.d0/pi
double precision,parameter ::km=1000.d0
integer, parameter :: MAX_STRING_LEN = 512
logical, parameter :: IGNORE_JUNK = .true.,DONT_IGNORE_JUNK = .false.
integer, parameter :: IIN = 40,IOUT = 41


