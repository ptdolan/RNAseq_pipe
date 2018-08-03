#!/bin/bash

#SBATCH --job-name=topHat
#
#SBATCH --time=10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH -o tophatStdOut.out
#SBATCH -e tophatError.out

sh /home/users/ptdolan/environs.sh
tophat /home/groups/jfrydman/Mus_musculus/Ensembl/NCBIM37/Sequence/Bowtie2Index /home/groups/jfrydman/TKW_RNAseq/*R1* /home/groups/jfrydman/TKW_RNAseq/*R2*


