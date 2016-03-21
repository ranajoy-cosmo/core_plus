#!/usr/bin/env python

import numpy as np

alpha = [0, 45, 90, 135]

n = alpha.size

cov_mat = np.zeros((3,3))

cov_mat[0,0] = n
cov_mat[0,1] = np.sum(np.cos(2*np.radians(alpha)))
cov_mat[0,2] = np.sum(np.sin(2*np.radians(alpha)))
cov_mat[1,0] = cov_mat[0,1]
cov_mat[1,1] = 0.5*(n + np.sum(np.cos(4*np.radians(alpha))))
cov_mat[1,2] = 0.5*np.sum(np.sin(4*np.radians(alpha)))
cov_mat[2,0] = cov_mat[0,2]
cov_mat[2,1] = cov_mat[1,2]
cov_mat[2,2] = 0.5*(n - np.sum(np.cos(4*np.radians(alpha))))

cov_mat /= 4

print cov_mat
