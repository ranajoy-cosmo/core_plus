#!/usr/bin/env python

import numpy as np
import sys

def generate_frequency_bands(freq_centre, band_width, edge_variation, num_bands): #GHz, fraction, fraction, number
    freq_low = freq_centre*(1.0 - band_width/2.0)
    freq_high = freq_centre*(1.0 + band_width/2.0)

    freq_low_list = np.random.normal(loc=freq_low, scale=edge_variation*freq_centre, size=num_bands)
    freq_high_list = np.random.normal(loc=freq_high, scale=edge_variation*freq_centre, size=num_bands)

    freq_centre_list = 0.5*(freq_low_list + freq_high_list)
    freq_width_list = (freq_high_list - freq_low_list)/freq_centre_list

    freq_low_mean = np.mean(freq_low_list)
    freq_high_mean = np.mean(freq_high_list)
    freq_centre_mean = np.mean(freq_centre_list)
    freq_width_mean = np.mean(freq_width_list)

    freq_low_variation = np.std(freq_low_list)/freq_centre_mean
    freq_high_variation = np.std(freq_high_list)/freq_centre_mean
    freq_centre_variation = np.std(freq_centre_list)/freq_centre_mean
    freq_width_variation = np.std(freq_width_list)

    print "#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*"
    print "#* INPUT BAND PROPERTIES"
    print "Central frequency : {0}".format(freq_centre)
    print "Band width : {0}".format(band_width)
    print "Low frequency : {0}".format(freq_low)
    print "High frequency : {0}".format(freq_high)
    print "Edge variation : {0}%".format(100*edge_variation)
    print "#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n"

    print "#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*"
    print "#* BAND PROPERTIES"
    print "Mean central frequency : {0:.3f}".format(freq_centre_mean)
    print "Mean low frequency : {0:.3f}".format(freq_low_mean)
    print "Mean high frequency : {0:.3f}".format(freq_high_mean)
    print "Mean band width : {0:.3f}".format(freq_width_mean)
    print "Variation low frequency : {0:.3f}%".format(100*freq_low_variation)
    print "Variation high frequency : {0:.3f}%".format(100*freq_high_variation)
    print "Variation centre frequency : {0:.3f}%".format(100*freq_centre_variation)
    print "Variation band width : {0:.3f}%".format(100*freq_width_variation)
    print "#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*"

    return freq_centre_list, freq_width_list

def print_frequency_bands(central_freq_list, band_width_list):
    central_freq_rounded = central_freq_list.round(2)
    band_width_rounded = band_width_list.round(2)
    print "Common elements :", len(central_freq_rounded) - len(set(central_freq_rounded))
    print central_freq_rounded
    print band_width_rounded


if __name__=="__main__":
    freq_centre = float(sys.argv[1])
    band_width = float(sys.argv[2])
    edge_variation = float(sys.argv[3])
    num_bands = float(sys.argv[4])

    freq_centre_list, freq_width_list = generate_frequency_bands(freq_centre, band_width, edge_variation, num_bands)
    print_frequency_bands(freq_centre_list, freq_width_list)
