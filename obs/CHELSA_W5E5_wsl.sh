#!/bin/bash

#SBATCH --job-name=W5E5
#SBATCH --chdir=/home/karger/
#SBATCH -A node # Node account
#SBATCH --qos normal  # normal priority level
#SBATCH --mail-user=dirk.karger@wsl.ch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=3G
#SBATCH --time=696:24:15
#SBATCH --output=logs/w5e5_%A_%a.out

Y1=$1 #start year W5E5
Y2=$2 # end year W5E5
Y3=$3 #the year before the startyear

YEAR=$(date -d "$Y3-12-31 $SLURM_ARRAY_TASK_ID days" +%Y)
MONTH=$(date -d "$Y3-12-31 $SLURM_ARRAY_TASK_ID days" +%m)
DAY=$(date -d "$Y3-12-31 $SLURM_ARRAY_TASK_ID days" +%d)
INPUT=$'/storage/karger/chelsa_V2/INPUT/'
OUTPUT=$'/storage/karger/chelsa_V2/W5E5/'
TEMP=$'/home/karger/scratch/'
SRAD=$'/storage/karger/chelsa_V2/OUTPUT_DAILY/srad/'
W5E5=$'/storage/karger/W5E5/'
CC=$'/storage/karger/chelsa_V2/OUTPUT_DAILY/tcc/'
LAPSE=$'/storage/karger/chelsa_V2/OUTPUT_DAILY/tz/'
PR=pr_W5E5v2.0_${Y1}0101-${Y2}1231.nc
HU=hurs_W5E5v2.0_${Y1}0101-${Y2}1231.nc
TA=tas_W5E5v2.0_${Y1}0101-${Y2}1231.nc
TX=tasmax_W5E5v2.0_${Y1}0101-${Y2}1231.nc
TN=tasmin_W5E5v2.0_${Y1}0101-${Y2}1231.nc
SR=rsds_W5E5v2.0_${Y1}0101-${Y2}1231.nc

srun singularity exec -B /storage /storage/karger/singularity/chelsa_V2.1.cont7 python /home/karger/scripts/chelsa_exchelsa/CHELSA_5E5.py -b $YEAR -c $MONTH -d $DAY -t1 $SLURM_ARRAY_TASK_ID -i $INPUT -o $OUTPUT -t $TEMP -s $SRAD -w $W5E5 -cc $CC -l $LAPSE -pr $PR -hu $HU -ta $TA -tx $TX -tn $TN -sr $SR --calc_pr --calc_rsds --calc_tas --calc_tasmax --calc_tasmin
