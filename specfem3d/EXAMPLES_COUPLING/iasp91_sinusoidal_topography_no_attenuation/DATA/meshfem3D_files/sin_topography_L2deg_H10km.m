clear
pi=3.141592653589793d0;
R_earth=6371000.0;
%The diameter of sinusoidal topography is 2 degrees (about 224 km).
L=R_earth*2*pi/360.0;
%The topography peak is 10000 meters.
H=10000;
degtorad=pi/180.0;

%lat_min
lat0=-3.0;
%lat_max
lat1=3.0;
%lon_min
lon0=-3.0;
%lon_max
lon1=3.0;
%lon space
dlat=0.01;
%lat space
dlon=0.01;
lat=lat0:dlat:lat1;
lon=lon0:dlon:lon1;
[lat_mesh,lon_mesh]=meshgrid(lat,lon);
nlat=length(lat);
nlon=length(lon);

for ilat=1:nlat
    for ilon=1:nlon
        [dist,az] = distance(lat(ilat),lon(ilon),0.0,0.0);
        dist_km=dist*R_earth*degtorad;
        if(dist_km<L/2.0)
            topo(ilat,ilon)=H/2.0*(cos(dist_km/(L/2.0)*pi)+1);
        else
            topo(ilat,ilon)=0.0;
        end
    end
end
 
fileID=fopen('latlon_surf_topo.dat','w');
fprintf(fileID,'%d \t %d \n',nlon,nlat);
fprintf(fileID,'%f \t %f \n',lon0,lat0);
fprintf(fileID,'%f \t %f \n',dlon,dlat);

for ilat=1:nlat
    for ilon=1:nlon
        fprintf(fileID,'%10.3f',topo(ilat,ilon));
    end
    fprintf(fileID,'\n');
end
fclose(fileID);
surf(lat_mesh,lon_mesh,topo);
%save latlon_surf_topo.dat topo -ASCII;
