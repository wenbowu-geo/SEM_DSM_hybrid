#!/bin/bash
#
dir_SEM=../../../specfem3d/EXAMPLES_COUPLING/iasp91_1D_benchmark
#dir_SEM=../../../specfem3d/EXAMPLES_COUPLING/iasp91_sinusoidal_topography
dir_DSM=../../../DSM/example/iasp91_1D_benchmark

#parameters in DATA/Par_file
#How many points we take to taper the SEM displacements and tractions.
npoints_taper_SEM=100

#In this case, only vertical components are computed. 
#Thus you do not need to compute the Green's functions corresponding to f_theta and f_phi.

#compute vertical component? If .true., make sure you have computed the Green's functions associated with vertical single force f_r. 
vertical_component=.true.
#compute radial component? If .true., make sure you have computed the Green's functions associated with radial single force f_theta. 
radial_component=.false.
#compute tranverse component? If .true., make sure you have computed the Green's functions associated with tranverse single force f_phi. 
transverse_component=.false.


#Check the HPC information and provide a proper number of CPU processors below.
#Here, We assume each node composed of 16 CPUs. If you want to use the PBS script 
#generated here, NPROC must be larger than 16. Otherwise, the parameter NNODE (number of nodes)
#would be zero! OR, create your own PBS script.
NPROC=32

NNODE=`echo "$NPROC" | awk '{print int($1/16);}'`

echo "Coupling 3D SEM with 1D DSM `date`"
currentdir=`pwd`
mpiruncmd="mpirun"

# sets up directory structure in current example directoy
echo
echo "setting up example..."
echo

# cleans output files
if [ "$1" != "noclean" ]; then
  echo "cleaning seismograms/"
  rm  -rf seismograms/
  echo "cleaning Green's Functions"
#  rm  -rf DSM_input/
  echo "cleaning SEM outputs/"
  rm  -rf SEM_input/iproc*
  rm  -f SEM_input/package_list  
  rm  -f SEM_input/SEM_Par_Coupling


#create directories
mkdir -p OUTPUT_FILES
mkdir -p SEM_input
mkdir -p DSM_input OUTPUT_FILES seismograms

#vertical component
mkdir -p DSM_input/zcomp
mkdir -p DSM_input/zcomp/disp_solid
mkdir -p DSM_input/zcomp/sigma

fi


#copy SEM outputs
  echo
  echo "copy the SEM outputs into SEM_input/"
  echo "Large size of data, it takes time..."
  echo
  cp ${dir_SEM}/OUTPUT_FILES/DATABASES_MPI/package_list SEM_input/
  cp ${dir_SEM}/OUTPUT_FILES/DATABASES_MPI/SEM_Par_Coupling SEM_input/
  cp -r  ${dir_SEM}/OUTPUT_FILES/DATABASES_MPIiproc* SEM_input/
#Rename the directories DATABASES_MPIiproc*
  ls -d SEM_input/DATABASES_MPIiproc* | awk '{split($1,new_name,"DATABASES_MPI"); print "mv ",$1,"SEM_input/"new_name[2];}' |sh  


#copy DSM outputs
  echo
  echo "copy the DSM Green's Function into DSM_input/zcomp/"
  echo "Large size of data, it takes time..."
  cp  ${dir_DSM}/OUTPUT_GREENS/sigma/freq* DSM_input/zcomp/sigma
  cp  ${dir_DSM}/OUTPUT_GREENS/disp_solid/freq* DSM_input/zcomp/disp_solid 
  cp ${dir_DSM}/OUTPUT_GREENS/Green_Par_forConvolution DSM_input/Green_Par
  echo


#prepare the input file DATA/Par_file
  echo
  echo "Make the file DATA/Par_file"
  echo 
#Examples:
#npoints_taper_SEM =  100
#vertical_component= .true.
#radial_component= .false.
#transverse_component= .false.


  echo "npoints_taper_SEM = " $npoints_taper_SEM>DATA/Par_file
  echo "vertical_component= " $vertical_component >>DATA/Par_file
  echo "radial_component= " $radial_component >>DATA/Par_file
  echo "transverse_component= " $transverse_component >>DATA/Par_file




#Copy the executable file. Choose either linear interpolation or Lagrange interpolation. 
#Lagrange interpolation is more precise, but slower than linear interpolation.
echo "Copy the executable file"
#cp ../../src/Linear_piecewise_interpolation/coupling_integral .
cp ../../src/Lagrange_interpolation/coupling_integral .
echo


#prepare PBS job  script
#Your cluster may provide a different job submition tool. If so, change the below lines to fit your cluster.
echo "Preparing the PBS script for job submision"

                echo "submit_coupling.pbs"

                echo "#!/bin/bash" >submit_coupling.pbs
                echo "# Sample PBS for parallel jobs" >> submit_coupling.pbs
                #run this job using {NNODE} nodes and each node has 16 processors.
                #ppn represent processors per node.
                echo "#PBS -l nodes=${NNODE}:ppn=16,walltime=00:50:00" >> submit_coupling.pbs
                #specify your job name
                echo "#PBS -N 1D_iasp91_coupling" >>submit_coupling.pbs
                #specify your email address for report
                echo "#PBS -M your_email_address" >>submit_coupling.pbs
                echo "#PBS -e log.err" >>submit_coupling.pbs
                echo "#PBS -o log.log" >> submit_coupling.pbs
                echo "# Load all the needed modules (e.g. openmpi, intel et al.) ">>submit_coupling.pbs
                echo "module load openmpi intel" >> submit_coupling.pbs
                echo "ulimit -s unlimited " >> submit_coupling.pbs
                echo "cd $PWD " >> submit_coupling.pbs
                echo "$mpiruncmd  ./coupling_integral" >> submit_coupling.pbs


#submit the job
#sbatch submit_coupling.pbs
qsub submit_coupling.pbs
echo
echo
echo "Waiting for jobs Done! It takes a few minutes."
echo "When it is done, check your seismograms in the directory seismograms/."
echo "Keep in mind, only vertical components are valid."

