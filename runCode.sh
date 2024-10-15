#!/bin/bash

#SBATCH --partition=pm6-isw2,pm9-isw0,pm11-isw2
#SBATCH --job-name=coa25
#SBATCH --account=ehpc-reg-2023r03-178
#SBATCH --qos=ehpc-reg-2023r03-178
#SBATCH --time=96:00:00

#SBATCH --nodes=1
#SBATCH --ntasks=256
##SBATCH --ntasks-per-core=1
#SBATCH --mem=251G
#SBATCH -e job.%J.err
#SBATCH -o job.%J.out

##SBATCH --mail-type=ALL
##SBATCH --mail-user=v.sanjay@utwente.nl

source ~/.bash_shell

OhOut=("1.00e-02" "0.01778" "0.03162" "0.05623" "0.10000")
RhoIn="1e-3"

Rr="1e0"
LEVEL="12"

start="1025"
end="1029"

for i in `seq $start $end`;
do
cd $i
export OMP_NUM_THREADS=25
./coalescenceBubble ${OhOut[$i-$start]} $RhoIn $Rr $LEVEL 40.0 0.01
cd ..
done
