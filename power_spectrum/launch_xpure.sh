#!/bin/bash

################################################################################
# Preprocessing of the explicit mapmaker map
################################################################################
input_maps_dir=
nhits_map=${input_maps_dir}/spectrum/nhits.fits
binary_mask_map=${input_maps_dir}/spectrum/binary_mask.fits
pol_weight_map=${input_maps_dir}/spectrum/pol_weight.fits
hdf5_maps=/scratch2/scratchdirs/dpoletti/explicit_mapmaking/full_season_ra23/files/map_wn_sprng*hdf5
nside=2048

#python ~/pb/cmb/src/expinv/em_power_spectrum.py --maps $hdf5_maps --nhits $nhits_map --binary $binary_mask_map --weight $pol_weight_map --nside $nside --opath $input_maps_dir/spectrum/ --vec_map /scratch2/scratchdirs/dpoletti/explicit_mapmaking/full_season_ra23/files/map_full_ra23.hdf5
#python ~/pb/cmb/src/expinv/em_power_spectrum.py --maps $hdf5_maps --nside $nside --opath $input_maps_dir/spectrum/maps/wn_sprng/ --vec_map /scratch2/scratchdirs/dpoletti/explicit_mapmaking/full_season_ra23/files/map_full_ra23.hdf5

################################################################################
# Maps
################################################################################
export map_expr='/global/scratch2/sd/dpoletti/explicit_mapmaking/full_season_ra23/spectrum/maps/wn_sprng/*.fits' #regular expression identifying the fits maps

################################################################################
# Parameters identified by the tag variable
################################################################################
export tag=ipp_ra23_ap30
export nside=$nside
export apodized_length=30
export rundir=$input_maps_dir/spectrum/ # directory where all the folders and files will be created
export bin_file=~/pb/files_pb/binpb400_alt.fits # bin of the power spectrum
export beam_file=~/pb/files_pb/beam_RA23.fits # the power spectrum will be corrected for this beam

# Input maps, masks
export nhit=$nhits_map

## Intensity  
export binary_mask_I1=$binary_mask_map
export weight_I1=$nhits_map

## Polarization
export binary_mask_P1=$binary_mask_map
export weight_P1=$pol_weight_map

# Xpure mode: 0=xpol 1=pure 2=hybrid
export xpure_mode=1

################################################################################
# Job parameters
################################################################################
export nproc=240
export walltime="00:30:00"
export queue="debug"

################################################################################
# Xpure job
################################################################################
qsub ~/workspace/xpure/xpure.pbs -q $queue -V -l mppwidth=$nproc -l walltime=$walltime -N $tag -e ${rundir}/$tag.err -o ${rundir}/$tag.out

################################################################################
# Deleting the variables
################################################################################
unset map_files map_tags tag nproc nside apodized_length rundir bin_file beam_file nhit binary_mask_I1 weight_I1 binary_mask_P1 weight_P1 xpure_mode

