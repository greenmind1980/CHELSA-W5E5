#!/bin/bash

#SBATCH --qos=short
#SBATCH --partition=standard
##SBATCH --qos=priority
##SBATCH --partition=priority
#SBATCH --job-name=chelsa_w5e5
#SBATCH --account=proclias
##SBATCH --mail-user=slange@pik-potsdam.de
##SBATCH --mail-type=ALL,TIME_LIMIT
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --output=/p/projects/proclias/1km/logs/chelsa_w5e5/downscaling/ch5_%A_%a
#SBATCH --error=/p/projects/proclias/1km/logs/chelsa_w5e5/downscaling/ch5_%A_%a

module purge
module load singularity


YEAR=$1
DAY_OF_THE_YEAR=$SLURM_ARRAY_TASK_ID

w5e5_version=2.0
[ "$w5e5_version" = "1.0" ] && YE=2016 || YE=2019
DAY=$(date -d "$(($YEAR-1))-12-31 $DAY_OF_THE_YEAR days" +%d)
MONTH=$(date -d "$(($YEAR-1))-12-31 $DAY_OF_THE_YEAR days" +%m)
YEAR=$(date -d "$(($YEAR-1))-12-31 $DAY_OF_THE_YEAR days" +%Y)  # increases the year if necessary (if DAY_OF_THE_YEAR is large)
Y1=$(expr $YEAR - $YEAR % 10 + 1)
[ $Y1 -gt $YEAR ] && Y1=$(expr $Y1 - 10)
[ $Y1 -eq 2011 ] && Y2=$YE || Y2=$(($Y1+9))
T1=$(( 1 + ($(date -d "$YEAR-$MONTH-$DAY UTC" +%s) - $(date -d "$Y1-01-01 UTC" +%s) ) / (60 * 60 * 24) ))
echo
echo Downscaling $YEAR-$MONTH-$DAY based on W5E5 v$w5e5_version $Y1-$Y2, timestep $T1 ...
echo

PR=pr_W5E5v${w5e5_version}_${Y1}0101-${Y2}1231.nc
HU=hurs_W5E5v${w5e5_version}_${Y1}0101-${Y2}1231.nc
TA=tas_W5E5v${w5e5_version}_${Y1}0101-${Y2}1231.nc
TX=tasmax_W5E5v${w5e5_version}_${Y1}0101-${Y2}1231.nc
TN=tasmin_W5E5v${w5e5_version}_${Y1}0101-${Y2}1231.nc
SR=rsds_W5E5v${w5e5_version}_${Y1}0101-${Y2}1231.nc
U=/p/projects/climate_data_central/reanalysis/ERA5/uas/uas_1hr_ECMWF-ERA5_observation_${YEAR}010100-${YEAR}123123.nc
V=/p/projects/climate_data_central/reanalysis/ERA5/vas/vas_1hr_ECMWF-ERA5_observation_${YEAR}010100-${YEAR}123123.nc

data_dir_root=/p/projects/proclias/1km/data/chelsa_w5e5
OUTPUT=$data_dir_root/output/
INPUT=$data_dir_root/input/input_org/
SRAD=$data_dir_root/input/rsds/
CC=$data_dir_root/input/tcc/
LAPSE=$data_dir_root/input/tz/
W5E5=$data_dir_root/input/w5e5/ #/p/projects/climate_data_central/observation/W5E5/v$w5e5_version/
TEMP=/p/tmp/slange/chelsa_w5e5/

# pause for 1-10 (random) seconds, to prevent too many jobs starting at the same time, to prevent "Couldn't determine user account information" error
sleep $[ ( $RANDOM % 10 ) + 1 ]s
mkdir -p /tmp/singularity/mnt/session
time singularity exec -B /p \
$data_dir_root/singularity/chelsa_V2.1.cont7 \
python -u CHELSA_W5E5.py \
-b $YEAR \
-c $MONTH \
-d $DAY \
-t1 $T1 \
-i $INPUT \
-o $OUTPUT \
-t $TEMP \
-s $SRAD \
-w $W5E5 \
-cc $CC \
-l $LAPSE \
-pr $PR \
-hu $HU \
-ta $TA \
-tx $TX \
-tn $TN \
-sr $SR \
-u $U \
-v $V \
--calc_pr \
--calc_tas \
#--calc_tasmax \
#--calc_tasmin \
#--calc_rsds
