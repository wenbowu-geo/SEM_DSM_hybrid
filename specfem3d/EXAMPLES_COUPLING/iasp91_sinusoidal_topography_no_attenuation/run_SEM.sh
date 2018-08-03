#!/bin/bash
#

echo "running iasp91_sinusolidal_topography example: `date`"
currentdir=`pwd`
mpiruncmd="mpirun"

# sets up directory structure in current example directoy
echo
echo "setting up example..."
echo
echo

mkdir -p OUTPUT_FILES

# cleans output files
if [ "$1" != "noclean" ]; then
  echo "cleaning OUTPUT_FILES/"
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo "!!!!!!!!!!!!IMPORTANT!!!!!!!!!!!!!!!!"
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo "(1) Remember to clean your OUTPUT_FILES/ before submitting a job."
  echo "(2) Make sure your previous job has been finished or manually killed before resubmiting the job."
  echo "    Because two jobs running could write information into the same file and get things messed up."
  echo 
  echo
  echo

  rm  -rf OUTPUT_FILES/*
fi

# Get the number of processors
NPROC=`grep ^NPROC DATA/Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`
#NNODE=`echo "$NPROC" | awk '{print int($1/16);}'`

mkdir -p OUTPUT_FILES/DATABASES_MPI

for iproc in `seq 0 $[$NPROC-1]`;
do
      zeropad_iproc=$(printf "%04d" $iproc)
      `mkdir -p OUTPUT_FILES/DATABASES_MPIiproc$zeropad_iproc`
done




################################Three distcontinuities##################################
#For the IASP91 model, our SEM box (thickness of 120km) contains three discontinuities 
#(free surface, upper crust interface at 15 km and Moho at 35 km). You need to list 
#the topography file of each interface in the file DATA/meshfem3D_files/interfaces.dat.
#The topography files guide xmeshfem3D to generate meshes consistent with the 1D DSM model.


############Moho and upper crust interface
#Here, we assume no topography on the upper crust interface and Moho, so you can easily 
#go to the directory DATA/meshfem3D_files/ and run get_Moho_MidCrust.cmd  to generate topo_Moho.dat 
#and topo_Midcrust.dat. Note that the depths of these two interfaces have to be consistent
#with the dsm model descibed in DATA/dsm_model_input.


############Free surface
#make the file ./DATA/meshfem3D_files/topo_surf.dat.
echo 
echo
echo
echo "Converting topogrpahy in two coordinate systems"
echo "Topography(longitude,latitude)-> topography(XI,ETA)"
echo "And save them the file ./DATA/meshfem3D_files/topo_surf.dat"

#../../bin/xcubedsphere_topo ith_interface topogrphy_as_function_latlon topography_basement min_allowed_topography max_allowed_topography
../../bin/xcubedsphere_topo 3 ./DATA/meshfem3D_files/latlon_surf_topo.dat 0.0 -20000.0 30000.0

#The first argument "3" indicates which interface we are working on. 
#There are three interfaces listed in the file ./DATA/meshfem3D_files/interfaces.dat.
#We are working on the third one, the free surface.

#The second parameter "./DATA/meshfem3D_files/topo_surf.dat" shows the name of topography 
#model in the longitude-latitude coordinate.

#The third parameter "0.0" is the basement topography of the free surface. 
#The code ../../bin/xcubedsphere_topo tapers the topography to the basement value near 
#the boundaries of SEM box.

#The last two parameters are the minnimum and maximum allowed topography. Any topography out of them will be truncated.
#In this example, the peak topography is 10000 meters, so we set max_allowed_topography=30000.0
######################################################################################



#make the file DATA/meshfem3D_files/real_bathymetry_topography
#This command below create the file DATA/meshfem3D_files/real_bathymetry_topography. 
#It is similar with what we do in last command.
#However, topography is not truncated and tapered here. xmeshfem3D will read this file, based on which 
#the elements are defined as fluid or solid media type.
../../bin/xcubedsphere_topo_FindOcean ./DATA/meshfem3D_files/latlon_surf_topo.dat  320 320 -0.1        -0.1      0.01d0    0.01d0

module load openmpi intel
echo 
echo
echo "running mesher"
mpirun -np $NPROC ../../bin/xmeshfem3D
echo
echo

echo
echo
echo "running generate_databases"
mpirun -np $NPROC ../../bin/xgenerate_databases
echo
echo

echo
echo
echo "running solver"
echo "It takes about 5min-15min..."
#Depends on the cluster you are using, you may need to replace sbatch by qsub.
#sbatch run_specfem3d.pbs
#echo
#echo

mpirun -np $NPROC ../../bin/xspecfem3D


echo
echo
echo


echo
echo
echo "Done! Next, move to the DSM step and run DSM for Green's Function."
echo `date`
