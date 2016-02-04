#! /bin/bash

output_folder="/global/homes/b/banerji/leap/apps/simulated_timestreams/bolo/maps_and_spectra"
#ln -s ${output_folder}/masks/mask_ebex.fits nhits_map.fits
#ln -s ${output_folder}/masks/mask_ebex.fits pol_weight_map.fits
#ln -s ${output_folder}/masks/mask_ebex.fits binary_mask_map.fits
#ln -s ${output_folder}/maps_b_001/full_sky_1024_8.fits sky_map.fits

ln -s /global/homes/s/stolpovs/temp/cov_map_coverage_angspeed2.6_delta_az15.0_hp.fits nhits_map.fits
ln -s /global/homes/s/stolpovs/temp/weight_map_coverage_angspeed2.6_delta_az15.0_hp.fits pol_weight_map.fits
ln -s /global/homes/s/stolpovs/temp/binary_map_coverage_angspeed2.6_delta_az15.0_hp.fits binary_mask_map.fits
ln -s /global/homes/s/stolpovs/temp/map_0.fits sky_map.fits


#output_folder=/global/scratch2/sd/banerji/leap_output/2015-08-10--13-38-01_iqu_maps

#ln -s ${output_folder}/signal_map.fits sky_map.fits
#ln -s ${output_folder}/hit_map.fits nhits_map.fits 

python make_bins_and_beam.py
