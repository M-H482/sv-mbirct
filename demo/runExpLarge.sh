#!/bin/bash
#------ qsub option --------#
#PBS -q sg 
#PBS -l select=1:ngpus=1
#PBS -l walltime=240:00:00
#PBS -W group_list=c30636
#PBS -j oe
#------- Program execution -------#
nvidia-smi
lscpu

cd /work/c30636/jpeak/sv-mbirct/demo

bash ./runGeneratSysMat.sh s3c
for i in s4{a..c}; do
    	bash ./runGeneratSysMat.sh $i
done

bash ./runGeneratSysMat.sh r3c
for i in r4{a..c}; do
    	bash ./runGeneratSysMat.sh $i
done
