#!/bin/bash
#
PWD=`pwd`

#vertical single force corresponding to vertical hybrid synthetics
echo Vertical single force
echo
echo
echo
cd fr/
sh run_DSM_fr.sh
echo
echo
echo
echo
echo
echo


#radial single force
echo Radial single force
echo
echo
echo
cd ../ftheta/
sh run_DSM_ftheta.sh
echo
echo
echo
echo
echo
echo


#transverse single force
echo Transverse single force
echo
echo
echo
cd ../fphi
sh run_DSM_fphi.sh
echo
echo
echo
echo
echo
echo
