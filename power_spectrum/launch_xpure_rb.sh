#!/bin/bash

################################################################################
# Preprocessing of the explicit mapmaker map
################################################################################
input_maps_dir_time_stamp=
input_maps_dir=/global/homes/b/banerji/simulation/output/reconstruction/${input_maps_dir_time_stamp}
nhits_map=${input_maps_dir}/hitmap_out.fits
binary_mask_map=${input_maps_dir}/binary_mask_map.fits
pol_weight_map=${input_maps_dir}/pol_weight_map.fits
#hdf5_maps=/scratch2/scratchdirs/dpoletti/explicit_mapmaking/full_season_ra23/files/map_wn_sprng*hdf5
nside=1024

ls ${nhits_map}
ls ${binary_mask_map}
ls ${pol_weight_map}

#python ~/pb/cmb/src/expinv/em_power_spectrum.py --maps $hdf5_maps --nhits $nhits_map --binary $binary_mask_map --weight $pol_weight_map --nside $nside --opath $input_maps_dir/spectrum/ --vec_map /scratch2/scratchdirs/dpoletti/explicit_mapmaking/full_season_ra23/files/map_full_ra23.hdf5
#python ~/pb/cmb/src/expinv/em_power_spectrum.py --maps $hdf5_maps --nside $nside --opath $input_maps_dir/spectrum/maps/wn_sprng/ --vec_map /scratch2/scratchdirs/dpoletti/explicit_mapmaking/full_season_ra23/files/map_full_ra23.hdf5

################################################################################
# Maps
################################################################################
export map_expr=${input_maps_dir}/sky_map.fits

################################################################################
# Parameters identified by the tag variable
################################################################################
export tag=b_001_pure_mik
export nside=${nside}
export apodized_length=30
export rundir=${input_maps_dir}/spectrum/ # directory where all the folders and files will be created
export bin_file=${input_maps_dir}/bin_file.fits # bin of the power spectrum
export beam_file=${input_maps_dir}/beam_file.fits # the power spectrum will be corrected for this beam

ls ${bin_file}
ls ${beam_file}

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
qsub ${input_maps_dir}/xpure.pbs -q $queue -V -l mppwidth=$nproc -l walltime=$walltime -N $tag -e ${rundir}/$tag.err -o ${rundir}/$tag.out

################################################################################
# Deleting the variables
################################################################################
unset map_files map_tags tag nproc nside apodized_length rundir bin_file beam_file nhit binary_mask_I1 weight_I1 binary_mask_P1 weight_P1 xpure_mode

