#!/bin/bash
#------- Program execution -------#
nvidia-smi
lscpu

bash ./runGeneratSysMat.sh s4c
bash ./runGeneratSysMat.sh r4b
bash ./runGeneratSysMat.sh r4c
exit 1

for i in s{1..4}{a..c}; do
    	bash ./runGeneratSysMat.sh $i
done

for i in r{1..4}{a..c}; do
    	bash ./runGeneratSysMat.sh $i
done
