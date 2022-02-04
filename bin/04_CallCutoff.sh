#!/bin/bash

mkdir -p $2/Code/$1
mkdir -p $2/CallResult/$1/
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


## Genome
# Call sites
meRanCall -p 40 -fs5 6 -fs3 6 -rs5 6 -rs3 6 -o \$TMPDIR/ZJ_tmp/$1_Genome10xCall.txt -bam $2/result/$1/$1_meRanGh_genomeMap_sorted.bam -f $GNM -mBQ 30 -gref -rl 150 -sc 10 -cr 1 -mr 0.00001 -mcov 10
meRanCall -p 40 -fs5 6 -fs3 6 -rs5 6 -rs3 6 -o \$TMPDIR/ZJ_tmp/$1_Genome10xCall_3_Cutoff.txt -bam $2/result/$1/$1_3_Ccutoff_PE.bam -f $GNM -mBQ 30 -gref -rl 150 -sc 10 -cr 1 -mr 0.00001 -mcov 10
meRanCall -p 40 -fs5 6 -fs3 6 -rs5 6 -rs3 6 -o \$TMPDIR/ZJ_tmp/$1_Genome10xCall_3_Cutoff_FDR.txt -bam $2/result/$1/$1_3_Ccutoff_PE.bam -f $GNM -mBQ 30 -gref -rl 150 -sc 10 -cr $GNconvrate -fdr 0.05 -mr 0.00001 -mcov 10

mv \$TMPDIR/ZJ_tmp/$1_Genome10xCall.txt $2/CallResult/$1/$1_Genome10xCall.txt
mv \$TMPDIR/ZJ_tmp/$1_Genome10xCall_3_Cutoff.txt $2/CallResult/$1/$1_Genome10xCall_3_Cutoff.txt
mv \$TMPDIR/ZJ_tmp/$1_Genome10xCall_3_Cutoff_FDR_FDR_0.05.txt $2/CallResult/$1/$1_Genome10xCall_3_Cutoff_FDR_FDR_0.05.txt

# Annotate sites
sed -i 's/chrM/chrMT/g' $2/CallResult/$1/$1_Genome10xCall.txt
sed -i 's/chr//' $2/CallResult/$1/$1_Genome10xCall.txt

sed -i 's/chrM/chrMT/g' $2/CallResult/$1/$1_Genome10xCall_3_Cutoff.txt
sed -i 's/chr//' $2/CallResult/$1/$1_Genome10xCall_3_Cutoff.txt 

meRanAnnotate -t $2/CallResult/$1/$1_Genome10xCall.txt -ensGTF -g $GTF -o $2/CallResult/$1/$1_Genome10xCall_annotate.txt -f 'gene'
meRanAnnotate -t $2/CallResult/$1/$1_Genome10xCall_3_Cutoff.txt -ensGTF -g $GTF -o $2/CallResult/$1/$1_Genome10xCall_3_Cutoff_annotate.txt -f 'gene'

" > $1_meRanCallCutoff.sbatch

sbatch $1_meRanCallCutoff.sbatch
