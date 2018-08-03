#!/bin/bash
#
PWD=`pwd`
dir_SEM=../../../../specfem3d/EXAMPLES_COUPLING/iasp91_1D_benchmark

#Number of porcessors used
NPROC=120
#Assmue each node composed of 16 CPU cores. If not, change it accordingly.
NNODE=`echo "$NPROC" | awk '{print int($1/16);}'`
distance_min=25.0
distance_max=35.0
distance_spacing=0.1
ndistance=`echo $distance_max $distance_min $distance_spacing | awk '{print int(($1-$2)/$3);}'`

echo "Running DSM to get Green's Function`date`"
currentdir=`pwd`
mpiruncmd="mpirun"

# sets up directory structure in current example directoy
echo
echo "setting up example..."
echo

# cleans output files
if [ "$1" != "noclean" ]; then
  echo "cleaning OUTPUT_FILES/"
  rm  -rf OUTPUT_GREENS/
fi


#create folders
mkdir -p OUTPUT_GREENS/
mkdir -p OUTPUT_GREENS/coef_cAnddcdr
mkdir -p OUTPUT_GREENS/disp_fluid
mkdir -p OUTPUT_GREENS/disp_solid
mkdir -p OUTPUT_GREENS/epsilon
mkdir -p OUTPUT_GREENS/potential
mkdir -p OUTPUT_GREENS/pressure
mkdir -p OUTPUT_GREENS/sigma
mkdir -p OUTPUT_GREENS/velo_solid


#copy the executable file
cp ../../../src/dsmti .

#copy the file dsm_input_file
echo
echo "copy the file dsm_input_file from the corresponding SEM folder"

#${dir_SEM}/DATA/dsm_model_input corresponds to vertical single force
cp ${dir_SEM}/DATA/dsm_model_input dsm_model_input
echo
echo


#copy the depth table of coupling GLL points in SEM simulations
cp ${dir_SEM}/OUTPUT_FILES/DATABASES_MPI/depth_table_elastic_innerbound depth_solid_list
cp ${dir_SEM}/OUTPUT_FILES/DATABASES_MPI/depth_table_acoustic_innerbound depth_fluid_list


ndepth_solid_media=`awk '{if(NR==1) print $2;}' depth_solid_list `
ndepth_fluid_media=`awk '{if(NR==1) print $2;}' depth_fluid_list `

echo $distance_min $distance_spacing $ndistance $ndepth_solid_media| awk ' {for(idep=1;idep<$4+1;idep++) { \
print $3; for(idist=1;idist<$3+1;idist++) print $1+(idist-1)*$2}}' >dist_solid_list

echo $distance_min $distance_spacing $ndistance $ndepth_fluid_media| awk ' {for(idep=1;idep<$4+1;idep++) { \
print $3; for(idist=1;idist<$3+1;idist++) print $1+(idist-1)*$2}}' >dist_fluid_list


#prepare PBS job  script
echo Preparing the PBS script for job submision
echo
echo

                echo submit_DSM.pbs

                echo "#!/bin/bash" >submit_DSM.pbs
                echo "# Sample PBS for parallel jobs" >> submit_DSM.pbs
                echo "#PBS -l nodes=${NNODE}:ppn=16,walltime=00:50:00" >> submit_DSM.pbs
	        #specify your job name
                echo "#PBS -N 1D_iasp91_DSMcomputing" >>submit_DSM.pbs
                #specify your email address for report
                echo "#PBS -M your_email_address" >>submit_DSM.pbs
                echo "#PBS -e log.err" >>submit_DSM.pbs
                echo "#PBS -o log.log" >> submit_DSM.pbs

                echo "# Load all the needed modules (e.g. openmpi, intel et al.) ">>submit_DSM.pbs
                echo "module load openmpi intel" >> submit_DSM.pbs
                echo "ulimit -s unlimited " >> submit_DSM.pbs
                echo "cd $PWD " >> submit_DSM.pbs
                echo "$mpiruncmd  ./dsmti <dsm_model_input" >> submit_DSM.pbs

#run PBS script

echo "submit PBS job"
qsub submit_DSM.pbs
#sbatch submit_DSM.pbs
echo
echo
echo "Waiting for job done! It takes about 20 mins."


echo `date`
