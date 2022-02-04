#!/bin/bash

pushd $2/Code/$1/

echo -e "#!/bin/bash
#SBATCH -p normal_q
#SBATCH -A davidxie_lab
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --output=$1_Cutoff.out
#SBATCH --mail-user zaustinj@vt.edu
#SBATCH --mail-type=END

module reset
module load Pysam

cd $2
python 03b_readWrite_C_cutoff.py $1 $2

" > $1_cutoff.sbatch

sbatch $1_cutoff.sbatch