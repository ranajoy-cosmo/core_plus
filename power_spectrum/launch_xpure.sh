#!/bin/bash

################################################################################
# Input maps
# xpure will be run on each of the match of the following regular expression.
# If the match is a fits map, an auto spectrum will be computed.
# If the match is folder, the script will interpret the files therein as maps 
# and will compute all the cross spectra between them 
################################################################################

################################################################################
# Parameters identified by the tag variable
################################################################################
export tag=test1
export map_dir=$SCRATCH/core_long_term_output/reconstruction/first_set/2016_02_01__20_33_59
export map_expr=$map_dir/reconstructed_map.fits
export rundir=$map_dir/xpure_out
if [ ! -d "$rundir" ];then
    mkdir $rundir
fi


export nside=1024
export lmax=2000
export apodized_length=30
export fwhm=8.0
python make_prelims.py $tag $map_dir $map_expr $lmax $fwhm
export bin_file=$map_dir/$tag"_bins.fits"
export beam_file=$map_dir/$tag"_beam.fits" # the power spectrum will be corrected for this beam

## Intensity  
export binary_mask_I1=$map_dir/$tag"_binary_mask.fits"
export weight_I1=$map_dir/$tag"_weight_map.fits"

## Polarization
export binary_mask_P1=$binary_mask_I1
export weight_P1=$weight_I1

# Xpure mode: 0=xpol 1=pure 2=hybrid
export xpure_mode=1

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
export walltime="00:30:00"
export queue="debug"

################################################################################
# Xpure job
################################################################################
export time_stamp="$(date "+%Y_%m_%d__%H_%M_%S")"
sbatch -p $queue --export=ALL -N $n_nodes -t $walltime -e ${rundir}/%j.err -o ${rundir}/%j.out xpure.sl

################################################################################
# Deleting the variables
################################################################################
unset map_files map_tags tag nproc nside apodized_length rundir bin_file beam_file nhit binary_mask_I1 weight_I1 binary_mask_P1 weight_P1 xpure_mode time_stamp lmax nside
