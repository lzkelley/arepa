#!/bin/sh
#SBATCH -p hernquist 
#SBATCH -J L75n910FP_tree 
#SBATCH -n 64 
#SBATCH -o OUTPUT.lsf
#SBATCH -e ERROR.lsf
#SBATCH --exclusive
#SBATCH --mail-user=mvogelsb@cfa.harvard.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mem-per-cpu=4000

minnum=0
maxnum=133

for ((i=$minnum; i<=$maxnum; i++))
do
./B-BaseTree param.txt $i 1>OUTPUT 2>ERROR 
done
