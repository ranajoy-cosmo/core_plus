import numpy as np
import healpy as hp
import os
import sys


def get_dim(pol_type):
    if pol_type == "TQU":
        dim = 3
    elif pol_type =="QU":
        dim = 2
    else:
        dim = 1

    ind_elements = dim*(dim+1)/2

    return dim, ind_elements


def make_block_matrix(flat_matrix, pol_type):
    func_dict = {'TQU' : make_block_matrix_TQU, 'QU' : make_block_matrix_QU, 'T' : make_block_matrix_T}
    block_matrix = func_dict[pol_type](flat_matrix)
    return block_matrix


def make_block_matrix_TQU(flat_matrix):
    num_blocks = flat_matrix.shape[0]

    block_matrix = np.empty((num_blocks, 3, 3))
    
    block_matrix[..., 0, 0] = flat_matrix[..., 0]
    block_matrix[..., 0, 1] = flat_matrix[..., 1]
    block_matrix[..., 0, 2] = flat_matrix[..., 2]
    block_matrix[..., 1, 0] = block_matrix[..., 0, 1] 
    block_matrix[..., 1, 1] = flat_matrix[..., 3]
    block_matrix[..., 1, 2] = flat_matrix[..., 4]
    block_matrix[..., 2, 0] = block_matrix[..., 0, 2] 
    block_matrix[..., 2, 1] = block_matrix[..., 1, 2] 
    block_matrix[..., 2, 2] = flat_matrix[..., 5]

    return block_matrix


def make_block_matrix_QU(flat_matrix):
    num_blocks = flat_matrix.shape[0]

    block_matrix = np.empty((num_blocks, 2, 2))
    
    block_matrix[..., 0, 0] = flat_matrix[..., 0]
    block_matrix[..., 0, 1] = flat_matrix[..., 1]
    block_matrix[..., 1, 0] = block_matrix[..., 0, 1] 
    block_matrix[..., 1, 1] = flat_matrix[..., 2]

    return block_matrix


def make_block_matrix_T(flat_matrix):
    num_blocks = flat_matrix.shape[0]

    block_matrix = np.empty((num_blocks, 1, 1))

    block_matrix[..., 0, 0] = flat_matrix[..., 0]

    return block_matrix


def get_inv_cov_matrix(hitpix, pol, ts, inv_cov_matrix, b_matrix, hitmap, npix, pol_type):
    func_dict = {'TQU' : get_inv_cov_matrix_TQU, 'QU' : get_inv_cov_matrix_QU, 'T' : get_inv_cov_matrix_T}
    func_dict[pol_type](hitpix, pol, ts, inv_cov_matrix, b_matrix, hitmap, npix)
    

def get_inv_cov_matrix_TQU(hitpix, pol, ts, inv_cov_matrix, b_matrix, hitmap, npix):
    n = np.bincount(hitpix, minlength=npix)
    cos_4 = np.bincount(hitpix, weights=np.cos(4*pol), minlength=npix)

    hitmap += n

    inv_cov_matrix[..., 0] += 0.25*n
    inv_cov_matrix[..., 1] += 0.25*np.bincount(hitpix, weights=np.cos(2*pol), minlength=npix)
    inv_cov_matrix[..., 2] += 0.25*np.bincount(hitpix, weights=np.sin(2*pol), minlength=npix)
    inv_cov_matrix[..., 3] += 0.25*0.5*(n + cos_4)
    inv_cov_matrix[..., 4] += 0.25*0.5*np.bincount(hitpix, weights=np.sin(4*pol), minlength=npix)
    inv_cov_matrix[..., 5] += 0.25*0.5*(n - cos_4)

    b_matrix[..., 0] += np.bincount(hitpix, weights=0.5*ts, minlength=npix)
    b_matrix[..., 1] += np.bincount(hitpix, weights=0.5*ts*np.cos(2*pol), minlength=npix)
    b_matrix[..., 2] += np.bincount(hitpix, weights=0.5*ts*np.sin(2*pol), minlength=npix)


def get_inv_cov_matrix_QU(hitpix, pol, ts, inv_cov_matrix, b_matrix, hitmap, npix):
    n = np.bincount(hitpix, minlength=npix)
    cos_4 = np.bincount(hitpix, weights=np.cos(4*pol), minlength=npix)

    hitmap += n

    inv_cov_matrix[..., 0] += 0.25*0.5*(n + cos_4)
    inv_cov_matrix[..., 1] += 0.25*0.5*np.bincount(hitpix, weights=np.sin(4*pol), minlength=npix)
    inv_cov_matrix[..., 2] += 0.25*0.5*(n - cos_4)

    b_matrix[..., 0] += np.bincount(hitpix, weights=0.5*ts*np.cos(2*pol), minlength=npix)
    b_matrix[..., 1] += np.bincount(hitpix, weights=0.5*ts*np.sin(2*pol), minlength=npix)


def get_inv_cov_matrix_T(hitpix, pol, ts, inv_cov_matrix, b_matrix, hitmap, npix):
    n = np.bincount(hitpix, minlength=npix)

    hitmap += n

    inv_cov_matrix[..., 0] += 0.25*n

    b_matrix[..., 0] += np.bincount(hitpix, weights=0.5*ts, minlength=npix)


def get_covariance_matrix(inv_cov_matrix, hitmap, pol_type):
    dim, ind_elements = get_dim(pol_type)
    inv_cov_matrix_block = make_block_matrix(inv_cov_matrix, pol_type)
    mask_bad_pix(inv_cov_matrix_block, hitmap, pol_type)
    
    cov_matrix = np.linalg.inv(inv_cov_matrix_block)
    
    return cov_matrix


def get_sky_map(cov_matrix, b_matrix, hitmap, pol_type):
    sky_map = np.sum(cov_matrix*b_matrix[..., None], axis=1).T
    if pol_type == "TQU":
        sky_map[..., hitmap<3] = np.nan
    elif pol_type == "QU":
        sky_map[..., hitmap<2] = np.nan
    else:
        sky_map[..., hitmap<1] = np.nan
    return sky_map


def mask_bad_pix(inv_cov_matrix, hitmap, pol_type):
    if pol_type == "TQU":
        inv_cov_matrix[hitmap<3] = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    elif pol_type == "QU":
        inv_cov_matrix[hitmap<2] = np.array([[1.0, 0.0], [0.0, 1.0]])
    else:
        inv_cov_matrix[hitmap<1] = np.array([[1.0]])


def write_maps(maps, map_type, field_names, recon_dir):

    hp.write_map(os.path.join(recon_dir, map_type+'.fits'), maps, column_names=field_names)

    if pol_type == "QU":
        hp.write_map(os.path.join(recon_dir, map_type+'.fits'), maps, column_names=['QQ', 'QU', 'UU'])
        #map_legends = {"QQ" : (0,0), "QU" : (0,1), "UU" : (1,1)}
    else:
        hp.write_map(os.path.join(recon_dir, map_type+'.fits'), maps, column_names=['TT', 'TQ', 'TU', 'QQ', 'QU', 'UU'])
        #map_legends = {"TT" : (0,0), "TQ" : (0,1), "TU" : (0,2), "QQ" : (1,1), "QU" : (1,2), "UU" : (2,2)}

    #for leg in map_legends.keys():
    #    hp.write_map(os.path.join(out_dir, "map_" + leg + ".fits"), maps[..., map_legends[leg][0], map_legends[leg][1]])

