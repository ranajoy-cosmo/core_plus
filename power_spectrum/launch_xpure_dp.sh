#!/bin/bash

################################################################################
# Input maps
# xpure will be run on each of the match of the following regular expression.
# If the match is a fits map, an auto spectrum will be computed.
# If the match is folder, the script will interpret the files therein as maps 
# and will compute all the cross spectra between them 
################################################################################
export map_expr='../full_sky/*map*'

################################################################################
# Parameters identified by the tag variable
################################################################################
export tag=ap30
export nside=2048
export apodized_length=30
export rundir="."
export bin_file=~/pb/files_pb/bin2_80_100_300_pb400.fits
export beam_file=~/pb/files_pb/beam_ziggy_2013.fits # the power spectrum will be corrected for this beam

## Intensity  
export binary_mask_I1="../full_sky/binary_mask.fits"
export weight_I1="../full_sky/weights.fits"

## Polarization
export binary_mask_P1=$binary_mask_I1
export weight_P1=$weight_I1

# Xpure mode: 0=xpol 1=pure 2=hybrid
export xpure_mode=2

################################################################################
# Job parameters
################################################################################
n_nodes=6
if [ "$NERSC_HOST" == "edison" ]
then
	proc_per_node=24
fi
if [ "$NERSC_HOST" == "cori" ]
then
	proc_per_node=32
fi
export nproc=$(( proc_per_node * n_nodes ))
export walltime="01:00:00"
export queue="regular"
n_jobs=2

################################################################################
# Xpure job
################################################################################
laucher="sbatch -p $queue --export=ALL -N $n_nodes -t $walltime -e ${rundir}/$tag.err -o ${rundir}/$tag.out xpure.sl"
echo $laucher
id=$(eval "$laucher")
id=${id##* }
for (( i=1; i<n_jobs; i++ )); do
    laucher="sbatch -p $queue --export=ALL -N $n_nodes -t $walltime -e ${rundir}/$tag.err -o ${rundir}/$tag.out --dependency=afternotok:$id xpure.sl"
    echo $laucher
    id=$(eval "$laucher")
    id=${id##* }
done

################################################################################
# Deleting the variables
################################################################################
unset map_files map_tags tag nproc nside apodized_length rundir bin_file beam_file nhit binary_mask_I1 weight_I1 binary_mask_P1 weight_P1 xpure_mode
