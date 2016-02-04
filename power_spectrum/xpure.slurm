#!/bin/bash

#SBATCH -p debug
#SBATCH -N 240 
#SBATCH -t 00:30:00 

####################################################
# Setting variables and creating directories
####################################################

# Directories
bindir=/global/homes/p/peloton/mapmaker/xpure/xpure/trunk/build/$NERSC_HOST
maskdir=${rundir}/mask/${tag} && mkdir -p $maskdir
outputdir=${rundir}/cls/${tag} && mkdir -p $outputdir
pardir=${rundir}/parameters/${tag} && mkdir -p $pardir
mlldir=${rundir}/mll/${tag} && mkdir -p $mlldir
cd $mlldir #cludge for storing the binned mll in the correct directory


# Variables
noise_for_T=0 # muK/pix                      
noise_for_P=0 # muK/pi                      
#lmax=$(( nside*2 ))
lmax=1500
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


#########################################################################################
# From binary mask to binary apodized mask
#########################################################################################
if [ -f $apodized_mask_I1 ] ; then 
    echo $apodized_mask_I1 found. Skipping the mask apodization
else
    echo "*** Apodization ***"
    echo "Apodizing temperature mask"
    time srun -n 1 ${bindir}/myapodizemask ${binary_mask_I1} ${apodized_mask_I1} -minpix 1 -inside 1 -radius ${apodized_length} &
    echo "Apodizing temperature mask"
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

    echo "Intensity spin functions"
    time aprun -n $((nproc / 2)) ${bindir}/scalar2spin ${pardir}/param_all_I${tag}.par >& ${pardir}/output_scalar2spinI${tag} &

    echo "Polarization spin functions"
    time aprun -n $((nproc / 2)) ${bindir}/scalar2spin ${pardir}/param_all_P${tag}.par >& ${pardir}/output_scalar2spinP${tag} &

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

    #time aprun -n $nproc ${bindir}/x2pure_create_mll ${pardir}/createMll.par
    time srun -n $nproc ${bindir}/x2pure_create_mll ${pardir}/createMll.par
fi

#########################################################################################
# XPURE spectrum computation
#########################################################################################
echo "*** Spectra computation ***"

# Number of xpure which will run in parallel:
(( n_nodes = nproc / 24 ))

#Create parameter file
n_maps=$(ls -1 $map_expr | wc -l)
map_files=($(ls $map_expr))
echo $n_maps
n_map_actually_run=0
for ((num=0; num<${n_maps}; num++)); do
    map_file=${map_files[num]}
    map_tag=${map_file##*/}
    map_tag=${map_tag%.*}
    if [ -f ${outputdir}/cellpure_${map_tag}*.fits ]; then 
        echo Skipping ${map_tag}
    else
        echo "************************ Spectrum ${map_tag} ************************"
        cat > ${pardir}/xpure_${map_tag}.par << _EOF_

mode = ${xpure_mode}

nside = ${nside}
nmaps = 1
nmasks = $nmask

bellfile1 = ${beam_file}

mapfile1 = ${map_file}

nhitfileT=$nhit
nhitfileP=$nhit
sigmaT1 = $noise_for_T ##in muK/pix, = noise polar/sqrt(2)
sigmaP1 = $noise_for_P

#inpBellfile = ${beam_file}

#CHOICE OF INPUT SPECTRUM
#inpCellfile = ${cl_input}


lmaxSim = ${lmax}

_EOF_

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

pseudofile = ${outputdir}/pseudopure_${map_tag}
cellfile = ${outputdir}/cellpure_${map_tag}

lmax = ${lmax}

_EOF_

        #RUN
        #if (( n_nodes > (n_maps - num % n_nodes) )); then 
        #    n_proc_xpure=$(( n_nodes / (n_maps - (num / n_nodes) * n_nodes) * 24 ))
        #else
            n_proc_xpure=24
        #fi
        while [ "$(jobs | grep Running | wc -l)" -eq "$n_nodes" ]; do
           sleep 5
        done 

        >&2 echo launching $n_map_actually_run th xpure on $n_proc_xpure processors. nodes are $n_nodes
        >&2 jobs
        #time aprun -n $n_proc_xpure ${bindir}/x2pure ${pardir}/xpure_${map_tag}.par &
        time srun -n $n_proc_xpure ${bindir}/x2pure ${pardir}/xpure_${map_tag}.par &
        >&2 echo $!
        (( n_map_actually_run++ ))
    fi ## To avoid to recompute the cell already done
done
wait
