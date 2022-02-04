#!/bin/bash

mkdir -p $2/Code/$1/Ch3_code
mkdir -p $2/CallResult/$1/
pushd $2/Code/$1/Ch3_code

echo -e "#!/bin/bash
#SBATCH -p normal_q
#SBATCH -A davidxie_lab
#SBATCH --nodes=1 --cpus-per-task=40
#SBATCH --exclusive
#SBATCH --time=48:00:00
#SBATCH --output=$1_meRanControlConv.txt
#SBATCH --mail-user zaustinj@vt.edu
#SBATCH --mail-type=END

module reset
module load SAMtools

## Genome
# Call sites
cd $2/Code/$1/Ch3_code
meRanCall -ccr -p 40 -o $2/CallResult/$1/$1_Genome10xCall_Controls.txt -bam $2/result/$1/$1_meRanGh_genomeMapMerged_sorted.bam -f /projects/epigenomicslab/Annotation/mm10.for.RNABS.fa -mBQ 30 -gref -rl 150 -cSeqID UnMethylated_Control

cd $2/Code/$1/Ch3_code
awk 'c&&!--c;/Analyzed: .*_meRanGh_genomeMap.*_sorted.bam/{c=18}' $1_meRanControlConv.txt > $1_GN_conv_rate.txt
sed -i 's/.*\t//' $1_GN_conv_rate.txt

" > $1_meRanCall_getConv.sbatch

sbatch $1_meRanCall_getConv.sbatch
