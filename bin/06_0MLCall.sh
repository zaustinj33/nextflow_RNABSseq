#!/bin/bash

pushd $2/Code/$1

GNconvrate=$(cat $2/Code/$1/$1_GN_conv_rate.txt)
GTF=/beegfs/projects/epigenomicslab/Annotation/Homo_sapiens.GRCh38.105.gtf
GNM=/beegfs/projects/epigenomicslab/Annotation/Homo_sapiens.GRCh38.dna.primary_assembly.fa

echo -e "#!/bin/bash
#SBATCH -p normal_q
#SBATCH -A davidxie_lab
#SBATCH --nodes=1 --cpus-per-task=40
#SBATCH --exclusive
#SBATCH --time=100:00:00
#SBATCH --mail-user zaustinj@vt.edu
#SBATCH --mail-type=END

mkdir -p \$TMPDIR/ZJ_tmp


meRanCall -p 40 -o \$TMPDIR/ZJ_tmp/$1_Genome10xCall_3_Cutoff_0ML.txt -bam $2/result/$1/$1_3_Ccutoff_PE.bam -f $GNM -mBQ 30 -gref -rl 150 -sc 10 -cr 1 -mr 0 -mcov 10 -regions $2/All_m5C.bed
mv \$TMPDIR/ZJ_tmp/$1_Genome10xCall_3_Cutoff_0ML.txt $2/CallResult/$1/$1_Genome10xCall_3_Cutoff_0ML.txt


" > $1_meRanCallCutoff_0ML.sbatch

sbatch $1_meRanCallCutoff_0ML.sbatch

