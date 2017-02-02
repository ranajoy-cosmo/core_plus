#!/usr/bin/env python

import numpy as np
import healpy as hp
import os
import sys
import shutil
import time
import importlib
from memory_profiler import profile
from mpi4py import MPI
from simulation.lib.data_management.data_utilities import get_local_bolo_segment_list
from simulation.timestream_simulation.bolo import Bolo
import simulation.lib.utilities.prompter as prompter
from simulation.lib.utilities.time_util import get_time_stamp 
import simulation.map_maker.covariance_matrix_utils as cov_ut

def run_mpi():
    print "Rank :", rank, "started"
    #Defining the dimensions of the arrays and matrices
    npix = hp.nside2npix(config.nside_out)
    dim, ind_elements = cov_ut.get_dim(config.pol_type)

    recon_dir = get_recon_dir()
    if rank == 0:
        make_data_dirs() 

    #The matrices local to the process
    inv_cov_matrix_local = np.zeros((npix, ind_elements), dtype=np.float)
    b_matrix_local = np.zeros((npix, dim), dtype=np.float)
    hitmap_local = np.zeros(npix, dtype=np.float)

    #Getting list of bolos and segments for the particular process
    bolo_segment_dict = get_local_bolo_segment_list(rank, size, config.bolo_list, config.segment_list)

    tot_seg = 0
    for keys in bolo_segment_dict.keys():
        tot_seg += len(bolo_segment_dict[keys])

    time.sleep(0.1*rank)
    print "Rank :", rank, ", Bolos and Segments :", bolo_segment_dict
    comm.Barrier()

    #Iterating over the bolos
    for bolo_name in bolo_segment_dict.keys():
        if config.subtract_template:
            TEMPLATE_name_0 = config.template_dict[bolo_name][0]
            TEMPLATE_name_1 = config.template_dict[bolo_name][1]
            bolo_TEMPLATE_0 = Bolo(TEMPLATE_name_0, config)
            bolo_TEMPLATE_1 = Bolo(TEMPLATE_name_1, config)
            estimated_y = np.load(os.path.join(config.general_data_dir, config.sim_tag, bolo_name+"_estimated_y.npy"))
        #If I am taking a differenced signal
        if config.take_diff_signal:
            bolo_a = Bolo(bolo_name + 'a', config)
            bolo_b = Bolo(bolo_name + 'b', config)
        else:
            bolo = Bolo(bolo_name, config)
        #Iterating over the segments in the bolo
        for segment in bolo_segment_dict[bolo_name]:
            start_seg = time.time()
            prompter.prompt("Rank : %d doing Bolo : %s and segment : %d" % (rank, bolo_name, segment))
            #Acquiring the signal
            if config.take_diff_signal:
                signal, v, pol_ang = acquire_difference_signal(bolo_a, bolo_b, segment, config.noise_only_map)
            else:
                signal, v, pol_ang = acquire_signal(bolo, segment, config.noise_only_map)
            if config.subtract_template:
                signal_TEMPLATE_0 = bolo_TEMPLATE_0.read_timestream(segment, read_list=["signal"])["signal"]
                signal_TEMPLATE_1 = bolo_TEMPLATE_1.read_timestream(segment, read_list=["signal"])["signal"]
                signal -= estimated_y[0]*signal_TEMPLATE_0 + estimated_y[1]*signal_TEMPLATE_1
            hitpix = hp.vec2pix(config.nside_out, v[...,0], v[...,1], v[...,2])
            del v
            #Generating the inverse covariance matrix
            cov_ut.get_inv_cov_matrix(hitpix, pol_ang, signal, inv_cov_matrix_local, b_matrix_local, hitmap_local, npix, config.pol_type)
            stop_seg = time.time()
            prompter.prompt("Rank : " + str(rank) + " Time taken : " + str(stop_seg - start_seg) + ". Projected time : " + str((stop_seg - start_seg)*tot_seg))

    #Saving space
    if config.subtract_template:
        del signal_TEMPLATE_0
        del signal_TEMPLATE_1
    del signal
    del pol_ang
    del hitpix
        
    start_dist_inv = time.time()

    #Distributing and gathering the segments of the matrices in the proper process
    #The sky pixels are chunked and are handled by individual processes
    inv_cov_matrix_local_segment = distribute_matrix(inv_cov_matrix_local, "cov_matrix")
    del inv_cov_matrix_local
    b_matrix_local_segment = distribute_matrix(b_matrix_local, "b_matrix")
    del b_matrix_local
    hitmap_local_segment = distribute_matrix(hitmap_local, "hitmap")
    del hitmap_local

    #Inverting the local segment of the inverse covariance matrix
    cov_matrix_local_segment = cov_ut.get_covariance_matrix(inv_cov_matrix_local_segment, hitmap_local_segment, config.pol_type)

    #Estimating the local sky segment
    sky_map_local_segment = cov_ut.get_sky_map(cov_matrix_local_segment, b_matrix_local_segment, hitmap_local_segment, config.pol_type)

    stop_dist_inv = time.time()
    prompter.prompt("Rank : " + str(rank) + " Time taken to distribute and invert: " + str(stop_dist_inv - start_dist_inv))

    write_segments(hitmap_local_segment, "hitmap", recon_dir) 
    write_segments(inv_cov_matrix_local_segment, "inverse_covariance_matrix", recon_dir) 
    write_segments(cov_matrix_local_segment, "covariance_matrix", recon_dir) 
    write_segments(sky_map_local_segment, "sky_map", recon_dir) 


def acquire_signal(bolo, segment, noise_only=False):
    if config.simulate_ts:
        t_stream = bolo.simulate_timestream(segment)
    else:
        t_stream = bolo.read_timestream(segment, noise_only=noise_only)

    return t_stream["signal"], t_stream["v"], t_stream["pol_ang"]
 

def acquire_difference_signal(bolo_a, bolo_b, segment, noise_only=False):
    if config.simulate_ts:
        t_stream_a = bolo_a.simulate_timestream(segment)
        t_stream_b = bolo_b.simulate_timestream(segment, noise_only=noise_only)
    else:
        t_stream_a = bolo_a.read_timestream(segment)
        t_stream_b = bolo_b.read_timestream(segment, read_list=["signal"], noise_only=noise_only)

    signal = 0.5*(t_stream_a["signal"] - t_stream_b["signal"])

    return signal, t_stream_a["v"], t_stream_a["pol_ang"]
        
#Where all the distribution and gathering takes place
def distribute_matrix(local_full_matrix, matrix_type):
    npix = hp.nside2npix(config.nside_out)
    dim, ind_elements = cov_ut.get_dim(config.pol_type)
    segment_length = npix/size 

    if matrix_type == "hitmap":
        local_segmented_matrix = np.zeros(segment_length) 
    else:
        inner_dim = {"cov_matrix" : ind_elements, "b_matrix" : dim}
        local_segmented_matrix = np.zeros((segment_length, inner_dim[matrix_type])) 

    for i in range(size):
        start = i*segment_length
        stop = (i+1)*segment_length
        comm.Reduce(local_full_matrix[start:stop], local_segmented_matrix, MPI.SUM, root=i)
        
    return local_segmented_matrix


def write_segments(maps, map_name, recon_dir):
    np.save(os.path.join(recon_dir, map_name + "_segments", str(rank).zfill(4)), maps) 


def get_recon_dir():
    sim_dir = os.path.join(config.general_data_dir, config.sim_tag)
    recon_dir = os.path.join(sim_dir, config.map_making_tag)

    return recon_dir

def get_segment_dirs():
    recon_dir = get_recon_dir()

    map_segment_dir = os.path.join(recon_dir, "sky_map_segments")
    hitmap_segment_dir = os.path.join(recon_dir, "hitmap_segments")
    cov_matrix_segment_dir = os.path.join(recon_dir, "covariance_matrix_segments")
    inv_cov_matrix_segment_dir = os.path.join(recon_dir, "inverse_covariance_matrix_segments")

    return map_segment_dir, hitmap_segment_dir, cov_matrix_segment_dir, inv_cov_matrix_segment_dir


def make_data_dirs():
    sim_dir = os.path.join(config.general_data_dir, config.sim_tag)
    if not os.path.exists(sim_dir):
        os.makedirs(sim_dir)

    time_stamp = get_time_stamp()

    recon_dir = os.path.join(sim_dir, config.map_making_tag)
    config_dir = os.path.join(recon_dir, "config_files")
    #default_config_file = os.path.join(config.base_dir, "map_maker/config_files/default_config.py")
    #current_config_file = os.path.join(config.base_dir, "map_maker/config_files", config_file + ".py")
    map_segment_dir, hitmap_segment_dir, cov_matrix_segment_dir, inv_cov_matrix_segment_dir = get_segment_dirs()

    if os.path.exists(recon_dir):
        shutil.rmtree(recon_dir)
    else:
        pass

    os.makedirs(recon_dir)
    os.makedirs(config_dir)
    os.makedirs(map_segment_dir)
    os.makedirs(hitmap_segment_dir)
    os.makedirs(inv_cov_matrix_segment_dir)
    os.makedirs(cov_matrix_segment_dir)

    """
    #shutil.copy(default_config_file, config_dir)
    #shutil.copyfile(current_config_file, os.path.join(config_dir, "config_file_" + time_stamp + ".py"))
    
    if config.timestream_data_products and config.simulate_ts:
        scan_dir = os.path.join(sim_dir, config.scan_tag)
        if not os.path.exists(scan_dir):
            os.makedirs(scan_dir)
        config_dir = os.path.join(scan_dir, "config_files")
        if not os.path.exists(config_dir):
            os.makedirs(config_dir)
        default_config_file = "/global/homes/b/banerji/simulation/map_maker/config_files/default_config.py"
        shutil.copy(default_config_file, config_dir)
        current_config_file = os.path.join("/global/homes/b/banerji/simulation/map_maker/config_files", config_file + '.py') 
        shutil.copy(current_config_file, config_dir)
    """

#The different sky chunks are brought together to form the final maps
def accumulate_segments(size):
    recon_dir = get_recon_dir()
    map_segment_dir, hitmap_segment_dir, cov_matrix_segment_dir, inv_cov_matrix_segment_dir = get_segment_dirs()

    npix = hp.nside2npix(config.nside_out)
    npix_segment = npix/size
    dim, ind_elements = cov_ut.get_dim(config.pol_type)

    sky_map = np.empty((dim, npix))

    for i in range(size):
        start = i*npix_segment
        stop = (i+1)*npix_segment
        #print os.path.join(recon_dir, map_segment_dir, str(i).zfill(4) + '.npy'), os.path.exists(os.path.join(recon_dir, map_segment_dir, str(i).zfill(4) + '.npy'))
        sky_map_segment = np.load(os.path.join(recon_dir, map_segment_dir, str(i).zfill(4) + '.npy')) 
        sky_map[..., start:stop] = sky_map_segment

    #if config.pol_type == 'T':
    #    sky_map[np.isnan(sky_map)] = 0.0
    #else:
    #    sky_map[..., np.isnan(sky_map[0])] = 0.0

    hp.write_map(os.path.join(recon_dir, "sky_map.fits"), sky_map)
    del sky_map
    del sky_map_segment
    shutil.rmtree(os.path.join(recon_dir, map_segment_dir))

    hitmap = np.empty(npix)

    for i in range(size):
        start = i*npix_segment
        stop = (i+1)*npix_segment
        hitmap_segment = np.load(os.path.join(recon_dir, hitmap_segment_dir, str(i).zfill(4) + '.npy')) 
        hitmap[..., start:stop] = hitmap_segment

    hp.write_map(os.path.join(recon_dir, "hitmap.fits"), hitmap)
    del hitmap
    del hitmap_segment
    shutil.rmtree(os.path.join(recon_dir, hitmap_segment_dir))

    inverse_cov_matrix = np.empty((npix, ind_elements))

    for i in range(size):
        start = i*npix_segment
        stop = (i+1)*npix_segment
        inverse_cov_matrix_segment = np.load(os.path.join(recon_dir, inv_cov_matrix_segment_dir, str(i).zfill(4) + '.npy')) 
        inverse_cov_matrix[start:stop] = inverse_cov_matrix_segment

    hp.write_map(os.path.join(recon_dir, "inverse_covariance_maps.fits"), inverse_cov_matrix.T)
    del inverse_cov_matrix
    del inverse_cov_matrix_segment
    shutil.rmtree(os.path.join(recon_dir, inv_cov_matrix_segment_dir))

    cov_matrix = np.empty((npix, ind_elements))

    for i in range(size):
        start = i*npix_segment
        stop = (i+1)*npix_segment
        cov_matrix_segment = np.load(os.path.join(recon_dir, cov_matrix_segment_dir, str(i).zfill(4) + '.npy')) 
        if config.pol_type == "TQU":
            cov_matrix_segment = cov_matrix_segment.reshape((npix_segment, dim**2))[..., np.array([0,1,2,4,5,8])]
        elif config.pol_type =="QU":
            cov_matrix_segment = cov_matrix_segment.reshape((npix_segment, dim**2))[..., np.array([0,1,3])]
        else: 
            cov_matrix_segment = cov_matrix_segment.reshape((npix_segment, dim**2))
        cov_matrix[start:stop] = cov_matrix_segment

    hp.write_map(os.path.join(recon_dir, "covariance_maps.fits"), cov_matrix.T)
    del cov_matrix
    del cov_matrix_segment
    shutil.rmtree(os.path.join(recon_dir, cov_matrix_segment_dir))


if __name__=="__main__":
    config_file = sys.argv[1]
    run_type = sys.argv[2]

    config = importlib.import_module(config_file).config

    if run_type=='accumulate_segments':
        size = int(sys.argv[3])
        accumulate_segments(size)

    if run_type=='run_mpi':
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        rank = comm.Get_rank()
        run_mpi()
        comm.Barrier()
        if rank==0:
            accumulate_segments(size)

    if run_type=='run_serial':
        run_serial()
