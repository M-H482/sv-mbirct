#!/bin/bash
#------- Program execution -------#
nvidia-smi
lscpu

bash ./runGeneratSysMat.sh case1
bash ./runGeneratSysMat.sh case2 
bash ./runGeneratSysMat.sh case3 

bash ./runGeneratSysMat.sh s4c use_stream
bash ./runGeneratSysMat.sh r4b use_stream
bash ./runGeneratSysMat.sh r4c use_stream
exit 1

for i in s{1..4}{a..c}; do
    	bash ./runGeneratSysMat.sh $i
done

for i in r{1..4}{a..c}; do
    	bash ./runGeneratSysMat.sh $i
done
