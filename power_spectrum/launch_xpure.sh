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
export tag=input_scan
export map_dir=/global/homes/b/banerji/leap/output/2016-03-05--09-18-58_simulated_timestreams
export map_expr=$map_dir/scanned_map_250.fits
export rundir=$map_dir/xpure_out
mkdir $rundir


export nside=2048
export lmax=2000
export apodized_length=30
export fwhm=8.0
python make_prelims.py $map_dir $map_expr $lmax $fwhm
export bin_file=$map_dir/bins.fits
export beam_file=$map_dir/beam.fits # the power spectrum will be corrected for this beam

## Intensity  
export binary_mask_I1=$map_dir/binary_mask.fits
export weight_I1=$map_dir/weight_map.fits

## Polarization
export binary_mask_P1=$binary_mask_I1
export weight_P1=$weight_I1

# Xpure mode: 0=xpol 1=pure 2=hybrid
export xpure_mode=0

################################################################################
# Job parameters
################################################################################
n_nodes=1
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
n_jobs=2

################################################################################
# Xpure job
################################################################################

#sbatch -p $queue --export=ALL -N $n_nodes -t $walltime -e $tag.err -o $tag.out xpure.sl
sbatch -p $queue --export=ALL -N $n_nodes -t $walltime -e ${rundir}/$tag.err -o ${rundir}/$tag.out xpure.sl

################################################################################
# Deleting the variables
################################################################################
unset map_files map_tags tag nproc nside apodized_length rundir bin_file beam_file nhit binary_mask_I1 weight_I1 binary_mask_P1 weight_P1 xpure_mode
