#!/bin/bash

module switch PrgEnv-gnu PrgEnv-intel

####################################################
# Setting variables and creating directories
####################################################

# Directories
bindir=/global/homes/p/peloton/mapmaker/xpure/xpure/trunk/build/$NERSC_HOST
maskdir=${rundir}/mask/${tag} && mkdir -p $maskdir
outputdir=${rundir}/cls/${tag} && mkdir -p $outputdir
pardir=${rundir}/parameters/${tag} && mkdir -p $pardir
mlldir=${rundir}/mll/${tag} && mkdir -p $mlldir


# Variables
noise_for_T=0 # muK/pix                      
noise_for_P=0 # muK/pi                      
lmax=$(( nside*2 ))
nmask=1

# Mask
apodized_mask_I1=${maskdir}/apodized_mask_I_${tag}.fits
apodized_mask_P1=${maskdir}/apodized_mask_P_${tag}.fits

# Spin maps
## Intensity
output_spin0_I1=${maskdir}/spin0_I_${tag}.fits
output_spin1_I1=${maskdir}/spin1_I_${tag}.fits
output_spin2_I1=${maskdir}/spin2_I_${tag}.fits
## Polarization
output_spin0_P1=${maskdir}/spin0_P_${tag}.fits
output_spin1_P1=${maskdir}/spin1_P_${tag}.fits
output_spin2_P1=${maskdir}/spin2_P_${tag}.fits


(( n_cpus = SLURM_CPUS_ON_NODE / 2 * SLURM_NNODES ))
echo "n_cpus are$n_cpus, SLURM_CPUS_ON_NODE $SLURM_CPUS_ON_NODE, SLURM_NNODES $SLURM_NNODES"

#########################################################################################
# From binary mask to binary apodized mask
#########################################################################################
if [ -f $apodized_mask_I1 ] ; then 
    echo $apodized_mask_I1 found. Skipping the mask apodization
else
    echo "*** Apodization ***"
    echo "Apodizing temperature mask"
    time srun -n 1 ${bindir}/myapodizemask ${binary_mask_I1} ${apodized_mask_I1} -minpix 1 -inside 1 -radius ${apodized_length} &
    echo "Apodizing polarization mask"
    time srun -n 1 ${bindir}/myapodizemask ${binary_mask_P1} ${apodized_mask_P1} -minpix 1 -inside 1 -radius ${apodized_length} &
    wait
fi


#########################################################################################
# Computation of spin window function
#########################################################################################
if [ -f $output_spin0_P1 ] ; then 
    echo ${output_spin0_P1} found. Skipping the spin functions computation
else
    # Intensity
    echo "*** Spin functions ***"
    cat > ${pardir}/param_all_I${tag}.par << EOF
##healpix parameters
###################
nside = $nside
lmax = $lmax                            #should not exceed 2*Nside

##input file parameters
#######################
maskBinary = ${binary_mask_I1}           #same input as for 'myapodizemask' 
window_spin0 = ${apodized_mask_I1}       #output of 'myapodizemask'
inverseNoise = ${weight_I1}      	#inverse noise weighting

##output file parameters
########################
output_spin0 = ${output_spin0_I1}
output_spin1 = ${output_spin1_I1}
output_spin2 = ${output_spin2_I1}
EOF

# Polarization
cat > ${pardir}/param_all_P${tag}.par << EOF
##healpix parameters
###################
nside = $nside
lmax = $lmax 				#should not exceed 2*Nside

##input file parameters
#######################
maskBinary = ${binary_mask_P1}  		#same input as for 'myapodizemask' 
window_spin0 = ${apodized_mask_P1}    	#output of 'myapodizemask'
inverseNoise = ${weight_P1}     		#inverse noise weighting

##output file parameters
########################
output_spin0 = ${output_spin0_P1}
output_spin1 = ${output_spin1_P1}
output_spin2 = ${output_spin2_P1}
EOF

if [ $(( SLURM_NNODES % 2 )) -eq 0 ]; then
    echo "Intensity spin functions, in parallel"
    time srun -n $((n_cpus / 2)) ${bindir}/scalar2spin ${pardir}/param_all_I${tag}.par >& ${pardir}/output_scalar2spinI${tag} &

    echo "Polarization spin functions, in parallel"
    time srun -n $((n_cpus / 2)) ${bindir}/scalar2spin ${pardir}/param_all_P${tag}.par >& ${pardir}/output_scalar2spinP${tag} &
else
    echo "Intensity spin functions"
    time srun -n $((n_cpus)) ${bindir}/scalar2spin ${pardir}/param_all_I${tag}.par >& ${pardir}/output_scalar2spinI${tag}
    wait

    echo "Polarization spin functions"
    time srun -n $((n_cpus)) ${bindir}/scalar2spin ${pardir}/param_all_P${tag}.par >& ${pardir}/output_scalar2spinP${tag} &
fi

    wait
fi

#########################################################################################
# X2PURE Mode-mode mixing matrix
#########################################################################################
if [ -f ${mlldir}/mll_TT_TT_BinMask1.fits ] ; then 
    echo ${mlldir}/mll_TT_TT_BinMask1.fits found. Skipping the mixing matrix computation
else
    echo "*** Mode mixing matrix ***"
    cat > ${pardir}/createMll.par << EOF

######### MODE #############
# 0 : Standard formalism
# 1 : Pure formalism
# 2 : Hybrid formalism
############################
mode = $xpure_mode

############ SETUP #########
nside = $nside
lmax = $lmax
nmask = $nmask

EOF

    spinI[0]=${output_spin0_I1}

    spinP0[0]=${output_spin0_P1}
    spinP1[0]=${output_spin1_P1}
    spinP2[0]=${output_spin2_P1}

    for(( i=1; i<=$nmask; i++)); do
            ind=$(($i - 1))
            cat >> ${pardir}/createMll.par << EOF

maskfile${i}_T  = ${spinI[${ind}]}

maskfile${i}_E_spin0 = ${spinP0[${ind}]}
maskfile${i}_E_spin1 = ${spinP1[${ind}]}
maskfile${i}_E_spin2 = ${spinP2[${ind}]} 

maskfile${i}_B_spin0 = ${spinP0[${ind}]}
maskfile${i}_B_spin1 = ${spinP1[${ind}]}
maskfile${i}_B_spin2 = ${spinP2[${ind}]}

mllfile_TT_TT_${i} = ${mlldir}/mll_TT_TT_BinMask${i}.fits

mllfile_EE_EE_${i} = ${mlldir}/mll_spinEE_EE_${i}.fits
mllfile_EE_BB_${i} = ${mlldir}/mll_spinEE_BB_${i}.fits
mllfile_EE_EB_${i} = ${mlldir}/mll_spinEE_EB_${i}.fits
mllfile_BB_BB_${i} = ${mlldir}/mll_spinBB_BB_${i}.fits
mllfile_BB_EE_${i} = ${mlldir}/mll_spinBB_EE_${i}.fits
mllfile_BB_EB_${i} = ${mlldir}/mll_spinBB_EB_${i}.fits

mllfile_TE_TE_${i} = ${mlldir}/mll_spinTE_TE_${i}.fits
mllfile_TE_TB_${i} = ${mlldir}/mll_spinTE_TB_${i}.fits
mllfile_TB_TE_${i} = ${mlldir}/mll_spinTB_TE_${i}.fits
mllfile_TB_TB_${i} = ${mlldir}/mll_spinTB_TB_${i}.fits

mllfile_EB_EB_${i} = ${mlldir}/mll_spinEB_EB_${i}.fits
mllfile_EB_EE_${i} = ${mlldir}/mll_spinEB_EE_${i}.fits
mllfile_EB_BB_${i} = ${mlldir}/mll_spinEB_BB_${i}.fits

EOF

    done

    time srun -n $n_cpus ${bindir}/x2pure_create_mll ${pardir}/createMll.par
fi

#########################################################################################
# XPURE spectrum computation
#########################################################################################
echo "*** Spectra computation ***"

# Number of xpure which will run in parallel:

#Create parameter file
spectrum_dir_or_files=($map_expr)
n_spectra=${#spectrum_dir_or_files[@]}
echo Number of spectra to be computed is $n_spectra
n_spectra_actually_run=0
for ((num=0; num<${n_spectra}; num++)); do
	spectrum_dir_or_file=${spectrum_dir_or_files[$num]}
	map_tag=${spectrum_dir_or_file##*/}
	map_tag=${map_tag%.*}  # Assumes no dot in in the name other than the extension
	map_files=($(ls $spectrum_dir_or_file))
	nmaps=${#map_files[@]}
	if [ "$nmaps" == 1 ] && [ -f ${outputdir}/*${map_tag}* ]; then 
		echo Skipping ${map_tag}
		continue
	fi
	if [ -d ${outputdir}/${map_tag} ]; then
		n_cross=(${outputdir}/${map_tag}/*.fits)
		n_cross=${#n_cross[@]}
		if (( n_cross == nmaps * (nmaps + 1) / 2 )); then
			echo Skipping ${map_tag}
			continue
		fi
	fi
	echo "************************ Spectrum ${map_tag} ************************"
	cat > ${pardir}/xpure_${map_tag}.par << _EOF_

mode = ${xpure_mode}

nside = ${nside}
nmaps = $nmaps
nmasks = $nmask


#nhitfileT=$nhit
#nhitfileP=$nhit
sigmaT1 = $noise_for_T ##in muK/pix, = noise polar/sqrt(2)
sigmaP1 = $noise_for_P

#inpBellfile = ${beam_file}

#CHOICE OF INPUT SPECTRUM
#inpCellfile = ${cl_input}


lmaxSim = ${lmax}

_EOF_


	for ((i_map=1; i_map<=nmaps; i_map++)); do
		map_file=${map_files[$((i_map-1))]}
		if [ -d $spectrum_dir_or_file ]; then
			mkdir -p ${outputdir}/$map_tag
			map_file=$spectrum_dir_or_file/$map_file
		fi
		cat >> ${pardir}/xpure_${map_tag}.par << _EOF_ 

bellfile${i_map} = ${beam_file}
mapfile${i_map} = ${map_file}

_EOF_

	done



	# Only one mask
	i=1
	ind=$(($i - 1))
	cat >> ${pardir}/xpure_${map_tag}.par << _EOF_

maskfile${i}_T  = $output_spin0_I1

maskfile${i}_E_spin0 = $output_spin0_P1
maskfile${i}_E_spin1 = $output_spin1_P1
maskfile${i}_E_spin2 = $output_spin2_P1

maskfile${i}_B_spin0 = $output_spin0_P1
maskfile${i}_B_spin1 = $output_spin1_P1
maskfile${i}_B_spin2 = $output_spin2_P1

mllfile_TT_TT_${i} = ${mlldir}/mll_TT_TT_BinMask${i}.fits

mllfile_EE_EE_${i} = ${mlldir}/mll_spinEE_EE_${i}.fits
mllfile_EE_BB_${i} = ${mlldir}/mll_spinEE_BB_${i}.fits
mllfile_EE_EB_${i} = ${mlldir}/mll_spinEE_EB_${i}.fits
mllfile_BB_BB_${i} = ${mlldir}/mll_spinBB_BB_${i}.fits
mllfile_BB_EE_${i} = ${mlldir}/mll_spinBB_EE_${i}.fits
mllfile_BB_EB_${i} = ${mlldir}/mll_spinBB_EB_${i}.fits

mllfile_TE_TE_${i} = ${mlldir}/mll_spinTE_TE_${i}.fits
mllfile_TE_TB_${i} = ${mlldir}/mll_spinTE_TB_${i}.fits
mllfile_TB_TE_${i} = ${mlldir}/mll_spinTB_TE_${i}.fits
mllfile_TB_TB_${i} = ${mlldir}/mll_spinTB_TB_${i}.fits

mllfile_EB_EB_${i} = ${mlldir}/mll_spinEB_EB_${i}.fits
mllfile_EB_EE_${i} = ${mlldir}/mll_spinEB_EE_${i}.fits
mllfile_EB_BB_${i} = ${mlldir}/mll_spinEB_BB_${i}.fits

noise_biasT_1 = 0.0
noise_biasP_1 = 0.0 
 
bintab = ${bin_file}

pseudofile = ${outputdir}/${map_tag}/pseudopure_${map_tag}
cellfile = ${outputdir}/${map_tag}/cellpure_${map_tag}

lmax = ${lmax}

_EOF_

	#RUN
	#if (( n_nodes > (n_maps - num % n_nodes) )); then 
	#    n_proc_xpure=$(( n_nodes / (n_maps - (num / n_nodes) * n_nodes) * 24 ))
	#else
	#fi
	n_nodes_per_xpure=2
	while (( SLURM_NNODES % n_nodes_per_xpure )); do
		n_nodes_per_xpure--
	done
	(( n_cpus_xpure = SLURM_CPUS_ON_NODE / 2 * n_nodes_per_xpure ))
	(( n_simulataneous_xpure = SLURM_NNODES / n_nodes_per_xpure ))

	while [ "$(jobs | grep Running | wc -l)" -eq "$n_simulataneous_xpure" ]; do
	   sleep 5
	done 

	>&2 echo launching $n_spectra_actually_run th xpure on $n_nodes_per_xpure  processors. nodes are $SLURM_NNODES
	>&2 jobs
	time srun -n $n_cpus_xpure --nodes $n_nodes_per_xpure ${bindir}/x2pure ${pardir}/xpure_${map_tag}.par &
	>&2 echo $!
	(( n_map_spectra_run++ ))
done
wait

mv $bin_file $beam_file $binary_mask_I1 $weight_I1 $rundir 
