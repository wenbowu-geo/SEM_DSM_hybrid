echo -35000.0 | gawk '{for(i=1;i<111;i++) {for(j=1;j<111;j++) print $1;}}' >topo_Moho.dat
echo -20000.0 | gawk '{for(i=1;i<111;i++) {for(j=1;j<111;j++) print $1;}}' >topo_Midcrust.dat

#echo -3878.0 | gawk '{for(i=1;i<541;i++) {for(j=1;j<5411;j++) print $1;}}' >topo_botOC_notopo.dat
#echo  0.0 | gawk '{for(i=1;i<541;i++) {for(j=1;j<541;j++) print $1;}}' >topo_surf_notopo.dat
