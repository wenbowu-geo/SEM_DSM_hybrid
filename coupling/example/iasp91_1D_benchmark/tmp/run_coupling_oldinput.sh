#!/bin/bash
#
dir_SEM=../../../specfem3d/EXAMPLES_COUPLING/1D_iasp91_benchmark
#dir_SEM=../../../specfem3d/EXAMPLES_COUPLING/iasp91_sinusoidal_topography
dir_DSM=../../../DSM/example/iasp91_1D_benchmark
npoints_taper_SEM=100
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
  rm  -f SEM_input/SEM_Par_file


mkdir -p OUTPUT_FILES
mkdir -p SEM_input
mkdir -p DSM_input OUTPUT_FILES seismograms
mkdir -p DSM_input/zcomp
mkdir -p DSM_input/zcomp/disp_solid
mkdir -p DSM_input/zcomp/sigma
mkdir -p seismograms


fi

#copy SEM outputs
  echo
  echo "copy the SEM outputs into SEM_input/"
  echo "Large size of data, it takes time..."
  echo
  cp ${dir_SEM}/OUTPUT_FILES/DATABASES_MPI/package_list SEM_input/
  cp ${dir_SEM}/OUTPUT_FILES/DATABASES_MPI/SEM_Par_Coupling SEM_input/
  cp ${dir_SEM}/DATA/Par_file SEM_input/SEM_Par_file
  cp ${dir_SEM}/DATA/meshfem3D_files/SEMtoTele_Par_file SEM_input/SEMtoTele_Par_file
  cp -r  ${dir_SEM}/OUTPUT_FILES/DATABASES_MPIiproc* SEM_input/
#Rename the folders DATABASES_MPIiproc*
   ls -d SEM_input/DATABASES_MPIiproc* | awk '{split($1,new_name,"DATABASES_MPI"); print "mv ",$1,"SEM_input/"new_name[2];}' |sh  


#copy DSM outputs
  echo
  echo "copy the DSM Green's Function into DSM_input/zcomp/"
  echo "Large size of data, it takes time..."
#  cp  ${dir_DSM}/OUTPUT_GREENS/sigma/freq* DSM_input/zcomp/sigma
#  cp  ${dir_DSM}/OUTPUT_GREENS/disp_solid/freq* DSM_input/zcomp/disp_solid 
#  cp ${dir_DSM}/OUTPUT_GREENS/green_par_forConvolution DSM_input/green_par
  echo


#prepare the input file DATA/Par_file
  echo
  echo "Make the file DATA/Par_file"
  echo 
  #counting the number of packages of SEM displacement and traction
  npackage_SEM=`wc -l <SEM_input/package_list`
#take out the parameters in SEM and DSM
  DT_SEM_orignal=`grep ^DT SEM_input/SEM_Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`
  NSTEP_SEM_orignal=`grep ^NSTEP SEM_input/SEM_Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`
  nstep_each_section_SEM=`grep ^NSTEP_BETWEEN_OUTPUTBOUND SEM_input/SEMtoTele_Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`
  decimate_factor=`grep ^DECIMATE_COUPLING SEM_input/SEMtoTele_Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`
  npoints_per_pack_SEM=`grep ^NPOINTS_PER_PACK SEM_input/SEMtoTele_Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`
  deltat_SEM=`echo $DT_SEM_orignal $decimate_factor |awk '{printf "%16.7f", $1*$2;}'`
  total_nstep_SEM=`echo $NSTEP_SEM_orignal $decimate_factor | awk '{if($1%$2==0) print int($1/$2); else print int($1/$2)+1;}'`
  nsection_SEM=`echo $total_nstep_SEM $nstep_each_section_SEM | awk '{if($1%$2==0) print int($1/$2); else print int($1/$2)+1;}'`

  time_length_DSM=`sed -n '5p'   DSM_input/green_par | tr -d ' '`
  nfreq_DSM=`sed -n '1p'   DSM_input/green_par | tr -d ' '`
  double_nfreq_DSM=$[2*$nfreq_DSM]
  deltat_DSM=`echo $time_length_DSM $double_nfreq_DSM | awk '{printf "%16.7f", $1/$2;}'`

#Examples:
#npackage_SEM = 107
#npoints_per_pack_SEM = 300
#nsection_SEM = 15
#nstep_each_section_SEM = 200
#total_nstep_SEM = 3000
#deltat_SEM = 0.05
#double_nfreq_DSM = 1024



  echo "npackage_SEM = " $npackage_SEM>DATA/Par_file
  echo "npoints_per_pack_SEM = " $npoints_per_pack_SEM>>DATA/Par_file
  echo "nsection_SEM = " $nsection_SEM>>DATA/Par_file
  echo "nstep_each_section_SEM = "$nstep_each_section_SEM>>DATA/Par_file
  echo "total_nstep_SEM = " $total_nstep_SEM>>DATA/Par_file
  echo "deltat_SEM = " $deltat_SEM>>DATA/Par_file
  echo "double_nfreq_DSM = " $double_nfreq_DSM>>DATA/Par_file
  echo "deltat_DSM = " $deltat_DSM>>DATA/Par_file
  echo "npoints_taper_SEM = " $npoints_taper_SEM>>DATA/Par_file



#Copy the executable file 
echo "Copy the executable file"
cp ../../src/Linear_piecewise_interpolation_oldinput/coupling_integral .
#cp ../../src/Lagrange_interpolation_oldinput/coupling_integral .
echo


#prepare PBS job  script
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
                echo "#PBS -M wenbow@princeton.edu" >>submit_coupling.pbs
                echo "#PBS -e log.err" >>submit_coupling.pbs
                echo "#PBS -o log.log" >> submit_coupling.pbs
                echo "# Load all the needed modules (e.g. openmpi, intel et al.) ">>submit_coupling.pbs
                echo "module load openmpi" >> submit_coupling.pbs
                echo "ulimit -s unlimited " >> submit_coupling.pbs
                echo "cd $PWD " >> submit_coupling.pbs
                echo "$mpiruncmd  ./coupling_integral" >> submit_coupling.pbs

sbatch submit_coupling.pbs


