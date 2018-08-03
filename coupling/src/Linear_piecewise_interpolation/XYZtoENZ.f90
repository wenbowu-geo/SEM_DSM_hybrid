subroutine XYZtoENZ(records,nfreq_Green,lat,lon)
    use constants
    implicit none
!imput/output
    integer ::nfreq_Green
    complex(kind=8),dimension(ncomp,nfreq_Green) ::records
    double precision ::lat,lon


!other variables
    double precision ::theta,phi
    complex(kind=8),dimension(ncomp) ::temp_copy
    double precision,dimension(ncomp) ::vertical,east,north
    integer ::ifreq

    theta=pi/2.d0-lat*degtorad
    phi=lon*degtorad
    vertical(1)=dsin(theta)*dcos(phi)
    vertical(2)=dsin(theta)*dsin(phi)
    vertical(3)=dcos(theta)
    north(1)=dcos(theta)*dcos(phi+pi)
    north(2)=dcos(theta)*dsin(phi+pi)
    north(3)=dsin(theta)
    east(1)=dcos(phi+pi/2.d0)
    east(2)=dsin(phi+pi/2.d0)
    east(3)=0.d0
   
    do ifreq=1,nfreq_Green
       temp_copy(:)=records(:,ifreq)
       records(1,ifreq)=sum(temp_copy*east)
       records(2,ifreq)=sum(temp_copy*north)
       records(3,ifreq)=sum(temp_copy*vertical)
    end do

end subroutine XYZtoENZ
