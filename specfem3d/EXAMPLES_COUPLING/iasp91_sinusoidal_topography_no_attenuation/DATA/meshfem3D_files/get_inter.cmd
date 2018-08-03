echo -35000.0 | gawk '{for(i=1;i<111;i++) {for(j=1;j<111;j++) print $1;}}' >topo_Moho.dat
echo -20000.0 | gawk '{for(i=1;i<111;i++) {for(j=1;j<111;j++) print $1;}}' >topo_Midcrust.dat
echo  0.0 | gawk '{for(i=1;i<321;i++) {for(j=1;j<321;j++) print $1;}}' >topo_flat_surf.dat
