#!/bin/bash

for i in $SCRATCH/*batchOutput/accepted_hits.bam; do
echo $i 

cat <<EOM > ${i/'batchOutput/accepted_hits.bam'/'cufflinks'}.slurm
#!/bin/bash
#SBATCH --job-name=cufflinksbatch
#
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G

set -o errexit

module purge
module load biology
module load bowtie2
module load tophat
module load cufflinks

cufflinks $i 

EOM
done;