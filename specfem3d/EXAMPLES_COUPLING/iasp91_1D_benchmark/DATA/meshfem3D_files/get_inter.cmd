echo -35000.0 | gawk '{for(i=1;i<111;i++) {for(j=1;j<111;j++) print $1;}}' >topo_Moho.dat
echo -20000.0 | gawk '{for(i=1;i<111;i++) {for(j=1;j<111;j++) print $1;}}' >topo_Midcrust.dat
echo  0.0 | gawk '{for(i=1;i<321;i++) {for(j=1;j<321;j++) print $1;}}' >topo_flat_surf.dat

echo  0.0 | gawk 'BEGIN{NX=100;NZ=100;orig_x=-0.2;orig_z=-0.2;dx=0.04;dz=0.04;} {print NX,NZ; print orig_x,orig_z; print dx,dz; for(iz=1;iz<NZ+1;iz++) {for(ix=1;ix<NX+1;ix++) print 0.0;}}' >real_bathymetry_topography
