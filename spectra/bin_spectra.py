#!/usr/bin/env python

import numpy as np

Cl = np.load("test_4000.npy")[0][:1999]
ell = np.arange(2, 2001)

Dl = ell*(ell + 1)*Cl/2/np.pi

bins = np.arange(2, 2001, 10)

Dl_sum = np.zeros(bins.size-1)
ell_sum = np.zeros(bins.size-1)

for i in range(bins.size - 1):
    Dl_sum[i] = np.sum(Dl[bins[i]-2:bins[i+1]-2])
    ell_sum[i] = np.sum(np.arange(bins[i], bins[i+1]))

Dl_binned = Dl_sum/ell_sum
