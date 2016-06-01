#!/usr/bin/env python

import numpy as np
import healpy as hp
import simulation.lib.quaternion.quaternion as qt

class Pointing:


    def __init__(self, bolo_params, scan_params):
        self._set_pointing_params(bolo_params, scan_params)
        self.get_initial_axes()

    def _set_pointing_params(self, bolo_params, scan_params):
        self._alpha = np.radians(scan_params.alpha)                         #radians
        self._beta = np.radians(scan_params.beta)                           #radians
        self._t_year = scan_params.t_year                                   #seconds
        self._t_prec = scan_params.t_prec                                  #seconds   
        self._t_spin = scan_params.t_spin                                   #seconds
        self._t_segment = scan_params.t_segment                             #seconds
        self._oversampling_rate = scan_params.oversampling_rate
        self._pol_phase_ini = np.radians(bolo_params.pol_phase_ini)         #radians
        self._pointing_offset_x = bolo_params.pointing_offset_x
        self._pointing_offset_y = bolo_params.pointing_offset_y

    def get_initial_axes(self):
        self._x_axis = np.array([1.0, 0.0, 0.0])
        self._y_axis = np.array([0.0, 1.0, 0.0])
        self._z_axis = np.array([0.0, 0.0, 1.0])
        self._axis_rev = self._z_axis
        self._axis_prec = self._x_axis
        self._axis_spin = np.array([np.cos(self._alpha), 0.0, np.sin(self._alpha)])

    def get_initial_vector(self, beta_offset):
        beta_offset_rad = np.radians(beta_offset)
        theta = self._alpha + self._beta + beta_offset_rad
        v_ini = np.array([np.cos(theta), 0.0, np.sin(theta)])

        return v_ini

    def get_local_axes_S2(self, v_pointing):
        theta, phi = hp.vec2dir(v_pointing)
        x_local = np.array(zip(np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi), -np.sin(theta)))
        y_local = np.array(zip(-np.sin(phi), np.cos(phi), np.zeros(phi.size)))

        return x_local, y_local

    def get_boresight_initial_vectors(self):
        theta = self._alpha + self._beta

        self._v_boresight_ini = get_initial_vectors(0)
        self._pol_vector_ini = np.array([np.sin(theta), 0.0, -np.cos(theta)])

    def generate_quaternion_boresight(pad_length=0):
        n_steps = self._t_segment*self._sampling_rate*self._oversampling_rate + 2*pad_length
        t_steps = t_start + np.arange(-pad_length, n_steps - pad_length, dtype=np.float)/(self._sampling_rate*self._oversampling_rate)

        w_spin = 2*np.pi/self._t_spin
        w_prec = 2*np.pi/self._t_prec
        w_rev = 2*np.pi/self._t_year

        self.quat_tot = qt.multiply(qt.make_quaternion(), qt.multiply(qt.make_quaternion(), qt.make_quaternion()))
