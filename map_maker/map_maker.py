#!/usr/bin/env python

import numpy as np
import healpy as hp
import os
import sys
import shutil
import time
import importlib
<<<<<<< HEAD
from memory_profiler import profile
from mpi4py import MPI
=======
#from memory_profiler import profile
from pysimulators.sparse import FSRBlockMatrix
from pysimulators import ProjectionOperator
>>>>>>> ebefaac9d61e24bda10f077c00b7a153e9574f50
from simulation.lib.data_management.data_utilities import get_local_bolo_segment_list
from simulation.timestream_simulation.bolo import Bolo
import simulation.lib.utilities.prompter as prompter
from simulation.lib.utilities.time_util import get_time_stamp 
<<<<<<< HEAD
import simulation.map_maker.covariance_matrix_utils as cov_ut
=======

def get_inv_cov_matrix(v, pol, ts, inv_cov_matrix, b_matrix):
    nsamples = hitpix.size
    npix = hp.nside2npix(config.nside_out)

    hitpix = hp.vec2pix(config.nside_out, v[...,0], v[...,1], v[...,2])

    n = np.bincount(hitpix, minlength=npix)
    cos_4 = np.bincount(hitpix, weights=np.cos(4*pol), minlength=npix)
    cos_2_pol = np.cos(2*pol)
    sin_2_pol = np.sin(2*pol)

    if config.pol_type == "TQU"
        inv_cov_matrix[..., 0] += 0.25*n
        inv_cov_matrix[..., 1] += 0.25*np.bincount(hitpix, weights=cos_2_pol, minlength=npix)
        inv_cov_matrix[..., 2] += 0.25*np.bincount(hitpix, weights=sin_2_pol, minlength=npix)
        inv_cov_matrix[..., 3] += 0.25*0.5*(n + cos_4) 
        inv_cov_matrix[..., 4] += 0.25*0.5*np.bincount(hitpix, weights=np.sin(4*pol), minlength=npix)
        inv_cov_matrix[..., 5] += 0.25*0.5*(n - cos_4) 

        b_matrix[..., 0] += np.bincount(hitpix, 0.5*ts, minlength=npix)
        b_matrix[..., 1] += np.bincount(hitpix, 0.5*cos_2_pol*ts, minlength=npix)
        b_matrix[..., 2] += np.bincount(hitpix, 0.5*sin_2_pol*ts, minlength=npix)
    else:
        inv_cov_matrix[..., 0] += 0.25*n
        inv_cov_matrix[..., 1] += 0.25*0.5*(n + cos_4) 
        inv_cov_matrix[..., 2] += 0.25*0.5*np.bincount(hitpix, weights=np.sin(4*pol), minlength=npix)
        inv_cov_matrix[..., 3] += 0.25*0.5*(n - cos_4) 

        b_matrix[..., 0] += np.bincount(hitpix, 0.5*cos_2_pol*ts, minlength=npix)
        b_matrix[..., 1] += np.bincount(hitpix, 0.5*sin_2_pol*ts, minlength=npix)


def write_maps_and_config(sky_map, hitmap, bad_pix, recon_dir):
    hp.write_map(os.path.join(recon_dir, "reconstructed_map.fits"), sky_map)
    hitmap[bad_pix] = 0
    hp.write_map(os.path.join(recon_dir, "hitmap.fits"), hitmap)


def write_covariance_maps(maps, map_type, recon_dir):
    out_dir = os.path.join(recon_dir, map_type)
    os.makedirs(out_dir)

    if map_type == "partial_covariance_maps" || config.pol_type == "pol_only":
        map_legends = {"QQ" : 0, "QU" : 1, "UU" : 2}
    else:
        map_legends = {"TT" : 0, "TQ" : 1, "TU" : 2, "QQ" : 3, "QU" : 4, "UU" : 5}

    for leg in map_legends.keys():
        hp.write_map(os.path.join(out_dir, "map_" + leg + ".fits"), maps[..., map_legends[leg][0], map_legends[leg][1]])
        #np.save(os.path.join(out_dir, "map_" + leg), maps[..., map_legends[leg][0], map_legends[leg][1]])
>>>>>>> ebefaac9d61e24bda10f077c00b7a153e9574f50


def get_signal(segment, bolo_a, bolo_b):
    if config.mode == "single_bolo"
        if config.simulate_ts:
            signal, v, pol_ang = bolo_a.simulate_timestream(segment)
        else:
            signal, v, pol_ang = bolo_a.read_timestream(segment)
    else:
        if config.simulate_ts:
            signal, v, pol_ang = bolo_a.simulate_timestream(segment)
            signal -= bolo_b.simulate_timestream(segment, return_field=["signal"])
        else:
            signal, v, pol_ang = bolo_a.read_timestream(segment)
            signal -= bolo_b.read_timestream(segment, return_field=["signal"])

    return signal, v, pol_ang


def run_mpi():
    print "Rank :", rank, "started"
    npix = hp.nside2npix(config.nside_out)
    dim, ind_elements = cov_ut.get_dim(config.pol_type)

<<<<<<< HEAD
    inv_cov_matrix_local = np.zeros((npix, ind_elements), dtype=np.float)
    b_matrix_local = np.zeros((npix, dim), dtype=np.float)
    hitmap_local = np.zeros(npix, dtype=np.float)
=======
    if config.pol_type == "TQU"
        inv_cov_matrix_local = np.zeros((npix, 6))
        b_matrix_local = np.zeros((npix, 3))
    else:
        inv_cov_matrix_local = np.zeros((npix, 3))
        b_matrix_local = np.zeros((npix, 2))
>>>>>>> ebefaac9d61e24bda10f077c00b7a153e9574f50

    bolo_segment_dict = get_local_bolo_segment_list(rank, size, config.bolo_list, config.segment_list)

    time.sleep(0.1*rank)
    print "Rank :", rank, ", Bolos and Segments :", bolo_segment_dict
    comm.Barrier()

    recon_dir = get_recon_dir()
    if rank == 0:
        make_data_dirs() 

    if config.subtract_template:
        bolo_TEMPLATE = Bolo("bolo_TEMPLATE", config)
        estimated_y = np.load("estimated_y.npy")

    for bolo_name in bolo_segment_dict.keys():
<<<<<<< HEAD
        print "Rank :", rank, "Bolos class being generated"
        if config.take_diff_signal:
            bolo_a = Bolo(bolo_name + 'a', config)
            bolo_b = Bolo(bolo_name + 'b', config)
        else:
            bolo = Bolo(bolo_name, config)
=======
        if config.mode = "single_bolo":
            bolo_1 = Bolo(bolo_name, config)
            bolo_2 = None
        else:
            bolo_1 = Bolo(bolo_name + 'a', config)
            bolo_2 = Bolo(bolo_name + 'b', config)
>>>>>>> ebefaac9d61e24bda10f077c00b7a153e9574f50
        for segment in bolo_segment_dict[bolo_name]:
            prompter.prompt("Rank : %d doing Bolo : %s and segment : %d" % (rank, bolo_name, segment))
<<<<<<< HEAD
            if config.take_diff_signal:
                signal, v, pol_ang = acquire_difference_signal(bolo_a, bolo_b, segment)
            else:
                signal, v, pol_ang = acquire_signal(bolo, segment)
            if config.subtract_template:
                signal_TEMPLATE = bolo_TEMPLATE.read_timestream(segment, read_list=["signal"])["signal"]
                signal -= estimated_y*signal_TEMPLATE
            print "Rank :", rank, "Bolos signal read"
            hitpix = hp.vec2pix(config.nside_out, v[...,0], v[...,1], v[...,2])
            del v
            cov_ut.get_inv_cov_matrix(hitpix, pol_ang, signal, inv_cov_matrix_local, b_matrix_local, hitmap_local, npix, config.pol_type)
            print "Rank :", rank, "Inverse covariance matrix generated"

    if config.subtract_template:
        del signal_TEMPLATE
    del signal
    del pol_ang
    del hitpix

    inv_cov_matrix_local_segment = distribute_matrix(inv_cov_matrix_local, "cov_matrix")
    del inv_cov_matrix_local
    b_matrix_local_segment = distribute_matrix(b_matrix_local, "b_matrix")
    del b_matrix_local
    hitmap_local_segment = distribute_matrix(hitmap_local, "hitmap")
    del hitmap_local
=======
            signal, v, pol_ang = get_signal(segment, bolo_1, bolo_2)
            hitpix = hp.vec2pix(config.nside_out, v[...,0], v[...,1], v[...,2])
            get_inv_cov_matrix(hitpix, pol_ang, signal, inv_cov_matrix_local, b_matrix_local)
            segment_stop = time.time()
            prompter.prompt("Rank : %d doing Bolo : %s and segment : %d and time taken : %d" % (rank, bolo_name, segment, segment_stop - segment_start))

    if config.pol_type == "TQU"
        inv_cov_matrix= np.zeros((npix, 6))
        b_matrix= np.zeros((npix, 3))
    else:
        inv_cov_matrix= np.zeros((npix, 3))
        b_matrix= np.zeros((npix, 2))
>>>>>>> ebefaac9d61e24bda10f077c00b7a153e9574f50

    cov_matrix_local_segment = cov_ut.get_covariance_matrix(inv_cov_matrix_local_segment, hitmap_local_segment, config.pol_type)

<<<<<<< HEAD
    sky_map_local_segment = cov_ut.get_sky_map(cov_matrix_local_segment, b_matrix_local_segment, hitmap_local_segment, config.pol_type)

    write_segments(hitmap_local_segment, "hitmap", recon_dir) 
    write_segments(inv_cov_matrix_local_segment, "inverse_covariance_matrix", recon_dir) 
    write_segments(cov_matrix_local_segment, "covariance_matrix", recon_dir) 
    write_segments(sky_map_local_segment, "sky_map", recon_dir) 


def acquire_signal(bolo, segment):
    if config.simulate_ts:
        t_stream = bolo.simulate_timestream(segment)
    else:
        t_stream = bolo.read_timestream(segment)
=======
    del b_matrix_local
    del inv_cov_matrix_local

    hitpix = 4*inv_cov_matrix[..., 0]

    if config.pol_type == "TQU"
        inv_cov_matrix[hitpix<3] = np.array([1.0, 0.0, 0.0, 1.0, 0.0, 1.0])
    else:
        inv_cov_matrix[hitpix<3] = np.array([1.0, 0.0, 1.0])
    write_covariance_maps(inv_cov_matrix, "inverse_covariance_maps", recon_dir)


        cov_matrix = np.linalg.inv(inv_cov_matrix)
        write_covariance_maps(cov_matrix, "covariance_maps", recon_dir)

        if "partial_covariance_maps" in config.map_making_data_products:
            cov_matrix_partial = np.linalg.inv(inv_cov_matrix[..., 1:, 1:])
            write_covariance_maps(cov_matrix_partial, "partial_covariance_maps", recon_dir)

        sky_rec = np.sum(cov_matrix*b_matrix[..., None], axis=1).T

        sky_rec[..., bad_pix] = np.nan

        write_maps_and_config(sky_rec, 4*inv_cov_matrix[..., 0, 0], bad_pix, recon_dir)
        
        stop = time.time()

        prompter.prompt("Total time taken : %d" % (stop - start))


def run_serial():
    start = time.time()
    npix = hp.nside2npix(config.nside_out)

    if config.pol_type == "TQU"
        inv_cov_matrix_unravel = np.zeros((npix, 6))
        b_matrix = np.zeros((npix, 3))
    else:
        inv_cov_matrix_unravel = np.zeros((npix, 3))
        b_matrix = np.zeros((npix, 2))

    recon_dir = make_data_dirs() 

    for bolo_name in config.bolo_list: 
        if config.mode = "single_bolo":
            bolo_1 = Bolo(bolo_name, config)
            bolo_2 = None
        else:
            bolo_1 = Bolo(bolo_name + 'a', config)
            bolo_2 = Bolo(bolo_name + 'b', config)
        for segment in config.segment_list: 
            segment_start = time.time()
            prompter.prompt("Bolo : %s and segment : %d" % (bolo_name, segment))
            signal, v, pol_ang = get_signal(segment, bolo_1, bolo_2)
            hitpix = hp.vec2pix(config.nside_out, v[...,0], v[...,1], v[...,2])
            get_inv_cov_matrix(hitpix, pol_ang, signal, inv_cov_matrix_unravel, b_matrix)
            segment_stop = time.time()
            prompter.prompt("Bolo : %s and segment : %d and time taken : %d" % (bolo_name, segment, segment_stop - segment_start))

    if config.pol_type == "TQU"
        inv_cov_matrix= np.zeros((npix, 3, 3))
    else:
        inv_cov_matrix= np.zeros((npix, 3))

    hitpix = 4*inv_cov_matrix[..., 0]

    if config.pol_type == "TQU"
        inv_cov_matrix[hitpix<3] = np.array([1.0, 0.0, 0.0, 1.0, 0.0, 1.0])
    else:
        inv_cov_matrix[hitpix<3] = np.array([1.0, 0.0, 1.0])
    write_covariance_maps(inv_cov_matrix, "inverse_covariance_maps", recon_dir)
>>>>>>> ebefaac9d61e24bda10f077c00b7a153e9574f50

    return t_stream["signal"], t_stream["v"], t_stream["pol_ang"]
 

def acquire_difference_signal(bolo_a, bolo_b, segment):
    if config.simulate_ts:
        t_stream_a = bolo_a.simulate_timestream(segment)
        t_stream_b = bolo_b.simulate_timestream(segment)
    else:
        t_stream_a = bolo_a.read_timestream(segment)
        t_stream_b = bolo_b.read_timestream(segment, read_list=["signal"])

    signal = 0.5*(t_stream_a["signal"] - t_stream_b["signal"])

    return signal, t_stream_a["v"], t_stream_a["pol_ang"]
        

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

def accumulate_segments(size):
    recon_dir = get_recon_dir()
    map_segment_dir, hitmap_segment_dir, cov_matrix_segment_dir, inv_cov_matrix_segment_dir = get_segment_dirs()

    npix = hp.nside2npix(config.nside_out)
    npix_segment = npix/size
    dim, ind_elements = cov_ut.get_dim(config.pol_type)

    sky_map = np.empty((dim, npix))

<<<<<<< HEAD
    for i in range(size):
        start = i*npix_segment
        stop = (i+1)*npix_segment
        print os.path.join(recon_dir, map_segment_dir, str(i).zfill(4) + '.npy'), os.path.exists(os.path.join(recon_dir, map_segment_dir, str(i).zfill(4) + '.npy'))
        sky_map_segment = np.load(os.path.join(recon_dir, map_segment_dir, str(i).zfill(4) + '.npy')) 
        sky_map[..., start:stop] = sky_map_segment

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


=======
>>>>>>> ebefaac9d61e24bda10f077c00b7a153e9574f50
if __name__=="__main__":
    config_file = sys.argv[1]
    run_type = sys.argv[2]

    config = importlib.import_module(config_file).config

<<<<<<< HEAD
    if run_type=='accumulate_segments':
        size = int(sys.argv[3])
        accumulate_segments(size)

=======
>>>>>>> ebefaac9d61e24bda10f077c00b7a153e9574f50
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
