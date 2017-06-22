#!/usr/bin/env python

import numpy as np
import healpy as hp
import os
import sys
import shutil
import time
import importlib
import pickle as pkl
from simulation.lib.utilities.generic_class import Generic
from memory_profiler import profile
from simulation.lib.data_management.data_utilities import get_local_bolo_segment_list
from simulation.lib.utilities.time_util import get_time_stamp
from simulation.timestream_simulation.bolo import Bolo
from simulation.lib.utilities.prompter import prompt
from simulation.lib.utilities.time_util import get_time_stamp 
import simulation.map_maker.covariance_matrix_utils as cov_ut


#@profile
def run_mpi():
    #Getting list of bolos and segments for the particular process
    bolo_segment_dict = get_local_bolo_segment_list(rank, size, config.bolo_list, config.segment_list)

    tot_seg = 0
    for keys in bolo_segment_dict.keys():
        tot_seg += len(bolo_segment_dict[keys])

    time.sleep(0.1*rank)
    prompt("Rank : {} \nBolos and Segments : {} \n# of segments : {}\n".format(rank, bolo_segment_dict, tot_seg))
    comm.Barrier()

    #Defining the dimensions of the arrays and matrices
    npix = hp.nside2npix(config.nside_out)
    dim, ind_elements = cov_ut.get_dim(config.pol_type)

    dir_names = get_dir_names()

    #The matrices local to the process
    inv_cov_matrix_local = np.zeros((npix, ind_elements), dtype=np.float)
    b_matrix_local = np.zeros((npix, dim), dtype=np.float)
    hitmap_local = np.zeros(npix, dtype=np.float)

    count = 0
    time_taken = 0
    #Iterating over the bolos
    for bolo_name in bolo_segment_dict.keys():
        if config.subtract_template:
            estimated_y = np.load(os.path.join(config.general_output_dir, config.sim_tag, bolo_name+"_estimated_y.npy"))
        #If I am taking a differenced signal
        if config.take_diff_signal:
            bolo_a = Bolo(bolo_name + 'a', config)
            bolo_b = Bolo(bolo_name + 'b', config)
        else:
            bolo = Bolo(bolo_name, config)
        #Iterating over the segments in the bolo
        for segment in bolo_segment_dict[bolo_name]:
            count += 1
            start_seg = time.time()
            prompt("Rank : {} doing Bolo : {} and segment : {}, {} of {}\n".format(rank, bolo_name, segment, count, tot_seg))
            #Acquiring the signal
            if config.sim_type == "signal":
                if config.take_diff_signal:
                    signal, hitpix, pol_ang = acquire_difference_signal(bolo_a, bolo_b, segment, config.noise_only_map)
                else:
                    signal, hitpix, pol_ang = acquire_signal(bolo, segment, config.noise_only_map)
                if config.subtract_template:
                    for i in range(len(config.TEMPLATE_list)):  
                        TEMPLATE_name = config.TEMPLATE_list[i]
                        TEMPLATE_signal = bolo_a.read_timestream(segment, return_field=[TEMPLATE_name])[TEMPLATE_name]
                        signal -= estimated_y[i] * TEMPLATE_signal
            if config.sim_type == "template":
                signal, hitpix, pol_ang = acquire_signal_template(bolo, segment)
            #Generating the inverse covariance matrix
            cov_ut.get_inv_cov_matrix(hitpix, pol_ang, signal, inv_cov_matrix_local, b_matrix_local, hitmap_local, npix, config.pol_type)
            stop_seg = time.time()
            time_taken += stop_seg - start_seg
            prompt("Rank : {}, Time taken : {}. Total time take : {}, Projected time : {}, Finished {} of {}\n".format(rank, stop_seg - start_seg, time_taken, time_taken*tot_seg/count, count, tot_seg))

    #Saving space
    if config.subtract_template:
        del TEMPLATE_signal
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
    write_segments(inv_cov_matrix_local_segment, "inverse_covariance_matrix", dir_names.recon_dir) 
    del inv_cov_matrix_local_segment

    #Estimating the local sky segment
    sky_map_local_segment = cov_ut.get_sky_map(cov_matrix_local_segment, b_matrix_local_segment, hitmap_local_segment, config.pol_type)

    stop_dist_inv = time.time()
    prompt("Rank : {}, Time taken to distribute and invert : {}\n".format(rank, stop_dist_inv - start_dist_inv))

    write_segments(hitmap_local_segment, "hitmap", dir_names.recon_dir) 
    write_segments(cov_matrix_local_segment, "covariance_matrix", dir_names.recon_dir) 
    write_segments(sky_map_local_segment, "sky_map", dir_names.recon_dir) 


def acquire_signal(bolo, segment, noise_only=False):
    if config.simulate_ts:
        t_stream = bolo.simulate_timestream_signal(segment, return_field=["signal", "pointing_vec", "pol_ang"])
    else:
        t_stream = bolo.read_timestream(segment, return_field=["signal", "pointing_vec", "pol_ang"], noise_only=noise_only)

    hitpix = hp.vec2pix(config.nside_out, t_stream["pointing_vec"][...,0], t_stream["pointing_vec"][...,1], t_stream["pointing_vec"][...,2])
    return t_stream["signal"], hitpix, t_stream["pol_ang"]

def acquire_signal_template(bolo, segment, tm_type=None):
    if config.template_type == "tm_bandpass":
        signal_type = config.tm_bandpass_type[0]
        t_stream = bolo.read_timestream(segment, return_field=[signal_type, "pointing_vec", "pol_ang"])
    if config.template_type == "tm_gradient":
        signal_type = config.tm_gradient_type[0]
        t_stream = bolo.read_timestream(segment, return_field=[signal_type, "pointing_vec", "pol_ang"])

    hitpix = hp.vec2pix(config.nside_out, t_stream["pointing_vec"][...,0], t_stream["pointing_vec"][...,1], t_stream["pointing_vec"][...,2])
    return t_stream[signal_type], hitpix, t_stream["pol_ang"]


def acquire_difference_signal(bolo_a, bolo_b, segment, noise_only=False):
    if config.simulate_ts:
        t_stream_a = bolo_a.simulate_timestream_signal(segment, return_field=["signal", "pointing_vec", "pol_ang"])
        t_stream_b = bolo_b.simulate_timestream_signal(segment, return_field=["signal"])
    else:
        t_stream_a = bolo_a.read_timestream(segment, return_field=["signal", "pointing_vec", "pol_ang"], noise_only=noise_only)
        t_stream_b = bolo_b.read_timestream(segment, return_field=["signal"], noise_only=noise_only)

    signal = 0.5*(t_stream_a["signal"] - t_stream_b["signal"])

    hitpix = hp.vec2pix(config.nside_out, t_stream_a["pointing_vec"][...,0], t_stream_a["pointing_vec"][...,1], t_stream_a["pointing_vec"][...,2])
    return signal, hitpix, t_stream_a["pol_ang"]
        
#Where all the distribution and gathering takes place
#@profile
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


def get_dir_names():
    dir_names = Generic()
    dir_names.sim_dir = os.path.join(config.general_output_dir, config.sim_tag)
    dir_names.scan_dir = os.path.join(dir_names.sim_dir, config.scan_tag)
    dir_names.recon_dir = os.path.join(dir_names.sim_dir, config.map_making_tag)
    dir_names.map_segment_dir = os.path.join(dir_names.recon_dir, "sky_map_segments")
    dir_names.hitmap_segment_dir = os.path.join(dir_names.recon_dir, "hitmap_segments")
    dir_names.cov_matrix_segment_dir = os.path.join(dir_names.recon_dir, "covariance_matrix_segments")
    dir_names.inv_cov_matrix_segment_dir = os.path.join(dir_names.recon_dir, "inverse_covariance_matrix_segments")

    return dir_names 


#The different sky chunks are brought together to form the final maps
#@profile
def accumulate_segments(size):
    acc_start_time = time.time()
    dir_names = get_dir_names()

    npix = hp.nside2npix(config.nside_out)
    npix_segment = npix/size
    dim, ind_elements = cov_ut.get_dim(config.pol_type)

    sky_map = np.empty((dim, npix))

    for i in range(size):
        start = i*npix_segment
        stop = (i+1)*npix_segment
        sky_map_segment = np.load(os.path.join(dir_names.recon_dir, dir_names.map_segment_dir, str(i).zfill(4) + '.npy')) 
        sky_map[..., start:stop] = sky_map_segment

    #if config.pol_type == 'T':
    #    sky_map[np.isnan(sky_map)] = 0.0
    #else:
    #    sky_map[..., np.isnan(sky_map[0])] = 0.0

    hp.write_map(os.path.join(dir_names.recon_dir, "sky_map.fits"), sky_map)
    del sky_map
    del sky_map_segment
    shutil.rmtree(os.path.join(dir_names.recon_dir, dir_names.map_segment_dir))

    hitmap = np.empty(npix)

    for i in range(size):
        start = i*npix_segment
        stop = (i+1)*npix_segment
        hitmap_segment = np.load(os.path.join(dir_names.recon_dir, dir_names.hitmap_segment_dir, str(i).zfill(4) + '.npy')) 
        hitmap[..., start:stop] = hitmap_segment

    hp.write_map(os.path.join(dir_names.recon_dir, "hitmap.fits"), hitmap)
    if config.pol_type == "TQU":
        mask_map = hitmap > 3
    elif config.pol_type == "QU":
        mask_map = hitmap > 2
    else:
        mask_map = hitmap > 1
    hp.write_map(os.path.join(dir_names.recon_dir, "mask.fits"), mask_map)

    del hitmap
    del hitmap_segment
    shutil.rmtree(os.path.join(dir_names.recon_dir, dir_names.hitmap_segment_dir))

    inverse_cov_matrix = np.empty((npix, ind_elements))

    for i in range(size):
        start = i*npix_segment
        stop = (i+1)*npix_segment
        inverse_cov_matrix_segment = np.load(os.path.join(dir_names.recon_dir, dir_names.inv_cov_matrix_segment_dir, str(i).zfill(4) + '.npy')) 
        inverse_cov_matrix[start:stop] = inverse_cov_matrix_segment

    hp.write_map(os.path.join(dir_names.recon_dir, "inverse_covariance_maps.fits"), inverse_cov_matrix.T)
    del inverse_cov_matrix
    del inverse_cov_matrix_segment
    shutil.rmtree(os.path.join(dir_names.recon_dir, dir_names.inv_cov_matrix_segment_dir))

    cov_matrix = np.empty((npix, ind_elements))

    for i in range(size):
        start = i*npix_segment
        stop = (i+1)*npix_segment
        cov_matrix_segment = np.load(os.path.join(dir_names.recon_dir, dir_names.cov_matrix_segment_dir, str(i).zfill(4) + '.npy')) 
        if config.pol_type == "TQU":
            cov_matrix_segment = cov_matrix_segment.reshape((npix_segment, dim**2))[..., np.array([0,1,2,4,5,8])]
        elif config.pol_type =="QU":
            cov_matrix_segment = cov_matrix_segment.reshape((npix_segment, dim**2))[..., np.array([0,1,3])]
        else: 
            cov_matrix_segment = cov_matrix_segment.reshape((npix_segment, dim**2))
        cov_matrix[start:stop] = cov_matrix_segment

    hp.write_map(os.path.join(dir_names.recon_dir, "covariance_maps.fits"), cov_matrix.T)
    del cov_matrix
    del cov_matrix_segment
    shutil.rmtree(os.path.join(dir_names.recon_dir, dir_names.cov_matrix_segment_dir))
    acc_stop_time = time.time()
    prompt("Time taken to accumulate segments : {}s\n".format(acc_stop_time - acc_start_time))


def make_data_dirs():
    dir_names = get_dir_names()
    if not os.path.exists(dir_names.sim_dir):
        try:
            os.makedirs(dir_names.sim_dir)
        except OSError:
            pass

    if not os.path.exists(dir_names.scan_dir):
        try:
            os.makedirs(dir_names.scan_dir)
        except OSError:
            pass

    if not os.path.exists(dir_names.recon_dir):
        os.makedirs(dir_names.recon_dir)
    else:
        shutil.rmtree(dir_names.recon_dir)
        os.makedirs(dir_names.recon_dir)

    os.makedirs(dir_names.map_segment_dir)
    os.makedirs(dir_names.hitmap_segment_dir)
    os.makedirs(dir_names.inv_cov_matrix_segment_dir)
    os.makedirs(dir_names.cov_matrix_segment_dir)

    config_dir = os.path.join(dir_names.sim_dir, "config_files")
    if not os.path.exists(config_dir):
        try:
            os.makedirs(config_dir)
        except OSError:
            pass
    this_config_dir = os.path.join(config_dir, time_stamp)
    try:
        os.makedirs(this_config_dir)
    except OSError:
        pass
    with open(os.path.join(this_config_dir, "config_file.pkl"), "w") as outfile:
        pkl.dump(config, outfile)
    bolo_config = importlib.import_module(config.bolo_config_file).bolo_config
    with open(os.path.join(this_config_dir, "bolo_config_file.pkl"), "w") as outfile:
        pkl.dump(bolo_config, outfile)

    code_dir = os.path.join(dir_names.sim_dir, "code_files")
    if not os.path.exists(code_dir):
        try:
            os.makedirs(code_dir)
        except OSError:
            pass
    this_code_dir = os.path.join(code_dir, time_stamp)
    try:
        os.makedirs(this_code_dir)
    except OSError:
        pass
    shutil.copyfile(os.path.join(config.base_dir, "map_maker", "map_maker.py"), os.path.join(this_code_dir, "map_maker.py"))
    shutil.copyfile(os.path.join(config.base_dir, "timestream_simulation", "bolo.py"), os.path.join(this_code_dir, "bolo.py"))
    shutil.copyfile(os.path.join(config.base_dir, "timestream_simulation", "beam_kernel.py"), os.path.join(this_code_dir, "beam_kernel.py"))


def start_message():
    display_string = "\n#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n"
    display_string += "#* BEGINNING SIMULATION\n"
    display_string += "#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n"
    display_string += "TIME STAMP : {}\n".format(time_stamp)
    display_string += "RUN TYPE : {}\n".format(run_type)
    sim_duration = config.t_segment*len(config.segment_list)
    display_string += "SIMULATION DURATION : {}s, {}years\n".format(sim_duration, sim_duration/365.25/24.0/60.0/60.0)
    display_string += "SCAN STRATEGY : {}\n".format(config.scan_strategy_name)
    display_string += "COORDINATE SYSTEM : {}\n".format(config.coordinate_system)
    display_string += "SIMULATE / READ TIMESTREAM DATA : {}\n".format("simulate" if config.simulate_ts else "read")
    display_string += "No. OF PROCESSES : {}\n".format(size)
    display_string += "DETECTOR LIST : {}\n".format(config.bolo_list)
    display_string += "No. OF DETECTORS : {}\n".format(len(config.bolo_list))
    display_string += "SEGMENT LIST : {}\n".format(config.segment_list)
    display_string += "No. OF SEGMENTS : {}\n".format(len(config.segment_list))
    display_string += "Segment length : {}s, {}h\n".format(config.t_segment, config.t_segment/60.0/60.0)
    display_string += "SIMULATION TYPE : {}\n".format(config.sim_type)
    if config.sim_type == "signal":
        if config.simulate_ts:
            display_string += "BEAM TYPE : {}\n".format(config.beam_type)
            display_string += "WRITE FIELD : {}\n".format(config.timestream_data_products)
    else:
        if config.template_type == "tm_gradient":
            display_string += "GRADIENT TYPES : {}\n".format(config.tm_gradient_type)
    display_string += "NOTES : {}\n".format(config.notes)
    display_string += "#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n\n"
    prompt(display_string, sys.stdout)

    if config.simulate_ts:
        bolo = Bolo(config.bolo_list[0], config)
        bolo.display_params()
        for bolo_name in config.bolo_list:
            bolo = Bolo(bolo_name, config)
            bolo.beam.display_beam_settings()

def run_check(verbose=True):
    if 12*config.nside_out**2 % size:
        if rank == 0 and verbose:
            prompt("# of processors is not compatible with the distribution of pixels on the segmented maps. Exiting") 
        sys.exit()

    if config.simulate_ts:
        if config.sim_type == "gradient":
            if config.sim_pol_type != 'T':
                if rank == 0 and verbose:
                    prompt("When computing gradients, the simulation polarisation is set to T only. Changing")
                config.sim_pol_type = 'T'


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* Main function definition. This is where the code begins when executed
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

if __name__=="__main__":
    config_file = sys.argv[1]
    run_type = sys.argv[2]

    config = importlib.import_module(config_file).config

    time_stamp = get_time_stamp()

    if run_type == 'accumulate_segments':
        size = int(sys.argv[3])
        accumulate_segments(size)

    if run_type == 'run_mpi':
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        rank = comm.Get_rank()
        comm.bcast(time_stamp, root=0)
        if rank == 0:
            make_data_dirs()
            start_message()
        comm.Barrier()
        run_check()
        run_mpi()
        comm.Barrier()
        if rank==0:
            accumulate_segments(size)

    if run_type == 'run_serial':
        make_data_dirs()
        start_message()
        run_serial()
