integer, parameter ::max_nfit=131
integer, parameter ::min_nfit_SEM_resample = 8
!nfit_Green must be even and larger than 4. Larger value of nfit_Green gives better interpolation, but it would slow down the compuation.
!We recommend a even number between 4-16. If the distance interval of Green function table is small, nfit_Green=4 is fine. Otherwise,
!use larger nift_Green to get better interpolation.
integer, parameter ::nfit_Green=16
integer, parameter ::max_npackages=4000
integer, parameter ::TINY=1.e-9
integer, parameter ::N_ELAS_COEF=21
integer, parameter ::ncomp=3
!Be careful to change max_array_size. It depends on your max-memory available.
integer, parameter ::max_array_size=700000000
integer, parameter ::ncomp_stress=6
double precision,parameter ::pi=3.1415926535897932d0
double precision,parameter ::degtorad=pi/180.d0
double precision,parameter ::r_earth=6371.0*1.0e3
double precision,parameter ::radtodeg=180.d0/pi
double precision,parameter ::km=1000.d0
integer, parameter :: MAX_STRING_LEN = 512
logical, parameter :: IGNORE_JUNK = .true.,DONT_IGNORE_JUNK = .false.
integer, parameter :: IIN = 40,IOUT = 41


