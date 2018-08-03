#!/bin/bash
#

echo "running 1D_iasp91 example: `date`"
currentdir=`pwd`
mpiruncmd="mpirun"

# sets up directory structure in current example directoy
echo
echo "setting up example..."
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

module load openmpi intel
echo "running mesher"
mpirun -np $NPROC ../../bin/xmeshfem3D
echo
echo

echo "running generate_databases"
mpirun -np $NPROC ../../bin/xgenerate_databases
echo
echo

echo "running solver"
echo "It takes about 5min-15min..."
mpirun -np $NPROC ../../bin/xspecfem3D
#sbatch run_specfem3d.pbs
#qsub run_specfem3d.pbs
echo 
echo


echo "Done! Next, move to the DSM step and run DSM for Green's Function."
echo `date`
