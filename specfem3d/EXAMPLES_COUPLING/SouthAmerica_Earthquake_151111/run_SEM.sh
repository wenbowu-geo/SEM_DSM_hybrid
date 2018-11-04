#!/bin/bash
#

echo "running SoutAmerica_Earthquake_151111 example: `date`"
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
#Here, We assume each node composed of 16 CPUs. If you want to use the PBS script 
#generated here, NPROC must be larger than 16. Otherwise, the parameter NNODE (number of nodes)
#would be zero! OR, create your own PBS script.

NNODE=`echo "$NPROC" | awk '{print int($1/16);}'`

mkdir -p OUTPUT_FILES/DATABASES_MPI

for iproc in `seq 0 $[$NPROC-1]`;
do
      zeropad_iproc=$(printf "%04d" $iproc)
      `mkdir -p OUTPUT_FILES/DATABASES_MPIiproc$zeropad_iproc`
done




################################Four distcontinuities##################################
#For the IASP91 model, our SEM box (thickness of 120km) contains three discontinuities 
#(free surface, middle crust interface at 15 km and Moho at 35 km). You need to list the 
#topography files corresponding to each one in the file DATA/meshfem3D_files/interfaces.dat. 
#The topography files guide xmeshfem3D to generate meshes consistent with the 1D DSM model.
#In addition to that, we need another interface, ocean bottom. This interface does not 
#exist in the DSM model, but we need it to incorporate the ocean in the SEM box. 


############Moho and middle crust interface
#Here, we assume no topography on the middle crust interface and Moho, so you can easily 
#go to the directory DATA/meshfem3D_files/ and run get_inter.cmd to generate topo_Moho.dat 
#and topo_Midcrust.dat. 


############Free surface and ocean bottom

###ocean bottom
#make the file ./DATA/meshfem3D_files/topo_botOC.dat.
echo 
echo
echo
echo "Converte topogrpahy between two coordinate systems"
echo "Topography(longitude,latitude)-> topography(XI,ETA)"
echo "First, ./DATA/meshfem3D_files/topo_botOC.dat"

#../../bin/xcubedsphere_topo ith_interface topogrphy_as_function_latlon topography_basement min_allowed_topography max_allowed_topography
../../bin/xcubedsphere_topo 3 ./DATA/meshfem3D_files/SA_real_topo.dat -3878.0 -20000.0 -1500.0

#The first argument "3" indicates which interface we are working on. 
#There are four interfaces listed in the file ./DATA/meshfem3D_files/interfaces.dat.
#We are working on the third one, the ocean bottom surface.

#The second parameter "./DATA/meshfem3D_files/SA_real_topo.dat" shows the name of topogrphy 
#model in the longitude-latitude coordinate.o

#The third parameter "-3878.0" is the basement bathymetry of the ocean bottom. 
#The code ../../bin/xcubedsphere_topo tapers the bathymetry 
#to the basement value near the boundaries of SEM box.

#The last two parameters are the minnimum and maximum allowed topography. Any 
#topography/bathymetry out of them will be truncated.
#In this example, we set max_allowed_topography=-1500.0. By doing this, the layer between the 
#free surface and ocean bottom has a minimum thickness of 1500 meters. 
#Note that, this ocean bottom interface extends to land and is fixed at a depth of 1500 m.



echo
echo 
echo "next, the file ./DATA/meshfem3D_files/topo_surf.dat"
###free surface
#make the file ./DATA/meshfem3D_files/topo_surf.dat.
echo 
../../bin/xcubedsphere_topo 4 ./DATA/meshfem3D_files/SA_real_topo.dat 0.0 0.0 30000.0
#the above command is similar with the former command, but works on the fouth interface, free surface. 
#The basement topography is 0.0 meter and the minimum topography is specified as 0.0 meter.


#####################################################################################





#make the file DATA/meshfem3D_files/real_bathymetry_topography
#The below command produces the file DATA/meshfem3D_files/real_bathymetry_topography. 
#It is similar to what we deal with the free surface above.
#However, topography is not tapered at SEM box boundaries. xmeshfem3D will read this 
#file and determine each elements as fluid (ocean) or solid media type.
../../bin/xcubedsphere_topo_FindOcean ./DATA/meshfem3D_files/SA_real_topo.dat 740 740 -0.2        -0.2      0.01d0    0.01d0


#***************************************************************************************
#prepare PBS job  script
echo Preparing PBS scripts

#Meshing
                echo run_meshfem3D.pbs

                echo "#!/bin/bash" >run_meshfem3D.pbs
                echo "# Sample PBS for parallel jobs" >> run_meshfem3D.pbs
                echo "#PBS -l nodes=${NNODE}:ppn=16,walltime=00:20:00" >> run_meshfem3D.pbs
                echo "#PBS -N your_job_name" >>run_meshfem3D.pbs
                #specify your email address for report
                echo "#PBS -M your_email_address" >>run_meshfem3D.pbs
                echo "#PBS -e log.err" >>run_meshfem3D.pbs
                echo "#PBS -o log.log" >> run_meshfem3D.pbs

                echo "# Load all the needed modules (e.g. openmpi, intel et al.) ">>run_meshfem3D.pbs
                echo "module load openmpi intel" >> run_meshfem3D.pbs
                echo "ulimit -s unlimited " >> run_meshfem3D.pbs
                echo "cd $currentdir " >> run_meshfem3D.pbs
                echo "$mpiruncmd  ../../bin/xmeshfem3D " >> run_meshfem3D.pbs

#Generating database
                # this is a crazy line, but with pure integer division its hard to handle.
                #set nodes = `echo ${nodnum} | awk '{printf "%.0f\n", $1/16+0.49}'`
                echo run_generate_databases.pbs

                echo "#!/bin/bash" >run_generate_databases.pbs
                echo "# Sample PBS for parallel jobs" >> run_generate_databases.pbs
                echo "#PBS -l nodes=${NNODE}:ppn=16,walltime=00:30:00" >> run_generate_databases.pbs
                echo "#PBS -N your_job_name" >>run_generate_databases.pbs
                #generate_databasesify your email address for report
                echo "#PBS -M your_email_address" >>run_generate_databases.pbs
                echo "#PBS -e log.err" >>run_generate_databases.pbs
                echo "#PBS -o log.log" >> run_generate_databases.pbs

                echo "# Load all the needed modules (e.g. openmpi, intel et al.) ">>run_generate_databases.pbs
                echo "module load openmpi intel" >> run_generate_databases.pbs
                echo "ulimit -s unlimited " >> run_generate_databases.pbs
                echo "cd $currentdir " >> run_generate_databases.pbs
                echo "$mpiruncmd  ../../bin/xgenerate_databases " >> run_generate_databases.pbs

#Runing solver
                echo run_specfem3d.pbs
                echo "#!/bin/bash" >run_specfem3d.pbs
                echo "# Sample PBS for parallel jobs" >> run_specfem3d.pbs
                echo "#PBS -l nodes=${NNODE}:ppn=16,walltime=18:00:00" >> run_specfem3d.pbs
                echo "#PBS -N your_job_name" >>run_specfem3d.pbs
                #specfem3dify your email address for report
                echo "#PBS -M your_email_address" >>run_specfem3d.pbs
                echo "#PBS -e log.err" >>run_specfem3d.pbs
                echo "#PBS -o log.log" >> run_specfem3d.pbs

                echo "# Load all the needed modules (e.g. openmpi, intel et al.) ">> run_specfem3d.pbs
                echo "module load openmpi intel" >> run_specfem3d.pbs
                echo "ulimit -s unlimited " >> run_specfem3d.pbs
                echo "cd $currentdir " >> run_specfem3d.pbs
                echo "$mpiruncmd  ../../bin/xspecfem3D " >> run_specfem3d.pbs


#openmpi intel works for the tiger cluster at princeton. You may need to load your own module.
module load openmpi intel
echo 
echo
echo "submit a job of meshing"
#submit a meshing job and save the job id
#use qsub, sbatch or other job submition tool. Depends on the setting of your cluster.
jobid_mesh=`qsub run_meshfem3D.pbs`
echo "Meshing job submited, running this job takes about 5 minutes..."




echo
echo
echo "submit a job of generating databases"
#submit a generate_database job and launch it after meshing is finished.
#use qsub, sbatch or other job submition tool. Depends on the setting of your cluster.
jobid_gen=`qsub -W depend=afterok:${jobid_mesh} run_generate_databases.pbs`
echo "The job of generating databases submited, running this job takes about 10 minutes..."



echo
echo
echo "sub a job of solver"
#submit a job of solver and launch it after generating database is finished.
#use qsub, sbatch or other job submition tool. Depends on the setting of your cluster.
jobid_spec=`qsub -W depend=afterok:${jobid_gen} run_specfem3d.pbs`
echo "The job of solver submited, running this job takes about 14 hours..."


echo
echo
echo
echo
echo "Waiting for jobs Done! It takes about 9 hours."
echo "Next, move to the DSM step and run DSM for Green's Function."
echo `date`
