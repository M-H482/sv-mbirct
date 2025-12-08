#!/bin/bash

# This is a script for running the "sv-mbirct" program for 
# parallel beam computed tomography. Specifically, it will 
# reconstruct sample data available for download at 
# "http://github.com/HPImaging/mbir-demos.git".
#
# More information on the command line usage can be found
# in the Readme file, and also by running "./mbir_ct -help".

# # Set to number of physical cores
export OMP_NUM_THREADS=20
export OMP_DYNAMIC=true
# export OMP_DYNAMIC=false

cd "$(dirname $0)"

### Set executable and data locations
execdir="../bin"

dataDir="./exp_data"
dataName=${1:-s1a}
echo "dataName = $dataName"

# sample sub-folder organization
parName="$dataDir/$dataName/par/$dataName"

matDir="./sysmatrix/exp_data"
matName="$matDir/$dataName"

# create folders that hold pre-computed items and output if they don't exist
if [[ ! -d "$matDir" ]]; then
  echo "Creating directory $matDir"
  mkdir "$matDir"
fi


### MBIR calls

# Compute system matrix and write to file (only have to run this once for a given geometry)
#   -i specifies image parameter file
#   -j specifies sinogram parameter file
#   -m specifies output system matrix file basename
#   -v verbose level 0=quiet, 1=show progress (default), 2=show even more

time $execdir/mbir_ct -i $parName -j $parName -m $matName -v 2
size=$(du -sh "${matName}.2Dsvmatrix"| cut -f1)
echo "system matrix size: $size"

time $execdir/mbir_ct_stream -i $parName -j $parName -m ${matName}_stream -v 2
size=$(du -sh "${matName}_stream.2Dsvmatrix"| cut -f1)
echo "system matrix size: $size"

echo "Comparing results~"
g++ -O3 verify_sysmatrix_float.cpp -o verify_diff
Nx=$(cat $parName.imgparams|grep Nx|cut -d' ' -f2)
Ny=$(cat $parName.imgparams|grep Ny|cut -d' ' -f2)
./verify_diff "${matName}.2Dsvmatrix" "${matName}_stream.2Dsvmatrix" $Nx $Ny


echo "Executing rm ${matName}.2Dsvmatrix"
rm "${matName}.2Dsvmatrix"
echo "Executing rm ${matName}_stream.2Dsvmatrix"
rm "${matName}_stream.2Dsvmatrix"
