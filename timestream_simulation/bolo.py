import numpy as np
import healpy as hp
import os
import shutil
import importlib
import sys
import time
from simulation.timestream_simulation.beam_kernel import Beam
from simulation.timestream_simulation.noise import Noise 
from simulation.lib.utilities.generic_class import Generic
from simulation.lib.utilities.prompter import prompt
import simulation.lib.quaternion.quaternion as quaternion
from memory_profiler import profile

"""
This module contains the Bolo class. The Bolo class is a generator of instances of individual bolometers with their independent properties given in the config file. The member functions of the class generate timestream signals for the particular bolometer configuration and 
"""

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#  
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

class Bolo:

    def __init__(self, bolo_name, config):
        self.config = Generic()
        bolo_config = importlib.import_module(config.bolo_config_file).bolo_config.bolos[bolo_name]
        self.config.name = bolo_name
        self.config.__dict__.update(config.__dict__)
        self.config.__dict__.update(bolo_config.__dict__)
        self.set_bolo_dirs()
        if config.simulate_ts:
            self.calculate_params()
            self.beam = Beam(self.config, bolo_config)
            self.noise_class = Noise(self.config)
            #sys.exit()
            if not config.sim_pol_type == "noise_only":
                self.get_sky_map()
            self.get_initial_axes()
            self.get_nsamples()
            if self.config.write_beam:
                self.beam.write_beam(self.bolo_dir)



#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Reading the timestream data for the bolo that was already simulated. 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def read_timestream(self, segment, return_field, noise_only=False):
        segment_dir = self.get_segment_dir(segment)
        t_stream = {}

        if "signal" in return_field:
            if noise_only:
                t_stream["signal"] = np.load(os.path.join(segment_dir, "noise.npy"))
            else:
                t_stream["signal"] = np.load(os.path.join(segment_dir, "signal.npy"))
            read_list.pop("signal")
        for return_item in return_field:
            t_stream[read_item] = np.load(os.path.join(segment_dir, read_item + ".npy"))

        return t_stream 

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Simulating the time-ordered data for a given bolo with any beam
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#    @profile
    def simulate_timestream(self, segment, return_field):
        t_return_stream = dict.fromkeys(return_field)
        write_field = self.config.timestream_data_products
        if write_field:
            self.make_write_dir(segment)

        #Generating the quaternion
        rot_quaternion = self.generate_quaternion(segment)

        #Generating the pointing and orientation vectors for the central line
        prompt("0.0\n", sys.stdout)
        
        vec_init = self.get_initial_vec(0.0)
        pointing_vec = self.get_vec_obv(vec_init, rot_quaternion)
        hitpix = hp.vec2pix(self.config.nside_in, pointing_vec[...,0], pointing_vec[...,1], pointing_vec[...,2])

        #Generating the polarisation angle of the detector on the sky
        pol_ang = self.get_pol_ang(rot_quaternion, pointing_vec)
        if self.config.sim_pol_type == "T" or self.config.sim_pol_type == "noise_only":
            cos2=None
            sin2=None
        else:
            cos2 = np.cos(2*pol_ang)
            sin2 = np.sin(2*pol_ang)

        if "pointing_vec" in write_field:
            self.write_timestream_data(pointing_vec, "pointing_vec", segment)
        if "pointing_vec" in return_field:
            t_return_stream["pointing_vec"] = pointing_vec
        elif (self.beam.del_beta) > 1:
            pass
        else:
            del pointing_vec

        if "pol_ang" in write_field:
            self.write_timestream_data(pol_ang, "pol_ang", segment)
        if "pol_ang" in return_field:
            t_return_stream["pol_ang"] = pol_ang
        else:
            del pol_ang


        if self.config.sim_type == "signal":
            if self.config.beam_type == "pencil":
                signal = self.generate_signal_pencil_beam(hitpix, cos2, sin2)
            elif self.config.beam_type in ["full_simulated", "from_file"]:
                beam_kernel_row = self.beam.get_beam_row(0.0)                       #The input argument is the beam offset from the centre
                signal = self.generate_signal_full_beam(hitpix, beam_kernel_row, cos2, sin2)
                #Iterating over the beam map and integrating
                for del_beta in self.beam.del_beta:
                    if del_beta == 0.0:
                        continue
                    prompt(str(del_beta), sys.stdout)
                    beam_kernel_row = self.beam.get_beam_row(del_beta)
                    vec_init = self.get_initial_vec(del_beta)
                    pointing_vec = quaternion.transform(rot_qt, v_init)
                    hitpix = hp.vec2pix(self.config.nside_in, pointing_vec[...,0], pointing_vec[...,1], pointing_vec[...,2])
                    signal += self.generate_signal_full_beam(hitpix, beam_kernel_row, cos2, sin2)
                beam_sum = np.sum(self.beam.beam_kernel[0])
                signal /= beam_sum
            else:
                prompt("Beam type not recognised", sys.stdout)
                sys.exit()
                    

            if self.config.noise_type != "none":
                noise = self.noise_class.simulate_timestream_noise_from_parameters()
                if "noise" in write_field:
                    self.write_timestream_data(noise, "noise", segment)
                if "noise" in return_field:
                    t_return_stream["noise"] = signal
                signal[::self.config.oversampling_rate] += noise 

            if "signal" in write_field:
                self.write_timestream_data(signal, "signal", segment)
            if "signal" in return_field:
                t_return_stream["signal"] = signal

        if self.config.sim_type == "gradient":
            signal = dict.fromkeys(["signal_top", "signal_central", "signal_bottom"])

            if list(set(["gradient_co", "gradient_coXgradient_co", "gradient_crossXgradient_cross"]) & set(self.config.gradient_type)):
                vec_init_central = self.get_initial_vec(0.0)
                pointing_vec = self.get_vec_obv(vec_init_central, rot_quaternion)
                hitpix = hp.vec2pix(self.config.nside_in, pointing_vec[...,0], pointing_vec[...,1], pointing_vec[...,2])
                signal["signal_central"] = self.generate_signal_pencil_beam(hitpix, None, None)
            if list(set(["gradient_cross", "gradient_crossXgradient_cross", "gradient_coXgradient_cross"]) & set(self.config.gradient_type)):
                vec_init_top = self.get_initial_vec(-self.config.scan_resolution)
                pointing_vec = self.get_vec_obv(vec_init_top, rot_quaternion)
                hitpix = hp.vec2pix(self.config.nside_in, pointing_vec[...,0], pointing_vec[...,1], pointing_vec[...,2])
                signal["signal_top"] = self.generate_signal_pencil_beam(hitpix, None, None)
                vec_init_bottom = self.get_initial_vec(self.config.scan_resolution)
                pointing_vec = self.get_vec_obv(vec_init_bottom, rot_quaternion)
                hitpix = hp.vec2pix(self.config.nside_in, pointing_vec[...,0], pointing_vec[...,1], pointing_vec[...,2])
                signal["signal_bottom"] = self.generate_signal_pencil_beam(hitpix, None, None)

                for gradient_type in list(set(write_field) | set(return_field)):
                    gradient = self.generate_signal_gradient(gradient_type, signal["signal_top"], signal["signal_central"], signal["signal_bottom"]) 
                    if gradient_type in write_field:
                        self.write_timestream_data(gradient, gradient_type, segment)
                    if gradient_type in return_field:
                        t_stream_gradient[gradient_type] = gradient

        if return_field:
            return t_return_stream


    #Generating timestream signal, no noise. Option of Intensity only, Polarisation only or all components
    #Although this promotes code duplicity, the signal generation subroutine is done separately for pencil beams, full beams and gradient calculation to not make the single subroutine too congested and difficult to understand and modify
    def generate_signal_pencil_beam(self, hit_pix, cos2, sin2):
        if self.config.sim_pol_type == "noise_only":
            signal = np.zeros(self.nsamples - 2*self.pad) 

        elif self.config.sim_pol_type == "T":
            signal = 0.5*self.sky_map[hit_pix]

        elif self.config.sim_pol_type in ["QU", "_QU"]:
            signal = 0.5*(self.sky_map[0][hit_pix]*cos2 + self.sky_map[1][hit_pix]*sin2) 

        else:
            signal = 0.5*(self.sky_map[0][hit_pix] + self.sky_map[1][hit_pix]*cos2 + self.sky_map[2][hit_pix]*sin2) 

        return signal

    def generate_signal_full_beam(self, hit_pix, beam_kernel_row, cos2, sin2):
        if self.config.sim_pol_type == "noise_only":
            signal = np.zeros(self.nsamples - 2*self.pad) 

        elif self.config.sim_pol_type == "T":
            signal = np.convolve(0.5*self.sky_map[hit_pix], beam_kernel_row[0], mode='valid')

        elif self.config.sim_pol_type in ["QU", "_QU"]:
            signal = np.convolve(0.5*(self.sky_map[0][hit_pix]*cos2 + self.sky_map[1][hit_pix]*sin2), beam_kernel_row[1], mode='valid')
            signal += np.convolve(0.5*(-1.0*self.sky_map[0][hit_pix]*sin2 + self.sky_map[1][hit_pix]*cos2), beam_kernel_row[2], mode='valid')

        else:
            signal = np.convolve(0.5*self.sky_map[0][hit_pix], beam_kernel_row[0], mode='valid')
            signal += np.convolve(0.5*(self.sky_map[1][hit_pix]*cos2 + self.sky_map[2][hit_pix]*sin2), beam_kernel_row[1], mode='valid')
            signal += np.convolve(0.5*(-1.0*self.sky_map[1][hit_pix]*sin2 + self.sky_map[2][hit_pix]*cos2), beam_kernel_row[2], mode='valid')

        return signal

    def generate_signal_gradient(self, gradient_type, signal_top=None, signal_central=None, signal_bottom=None): 
        if gradient_type == "gradient_co":
            gradient = (np.roll(signal_central, -1) - np.roll(signal_central, 1)) / (2*self.config.scan_resolution)
        elif gradient_type == "gradient_cross":
            gradient = (signal_top - signal_bottom) / (2*self.config.scan_resolution) 
        elif gradient_type == "gradient_coXgradient_co":
            gradient = (np.roll(signal_central, -1) - 2*signal_central + np.roll(signal_central, 1)) / self.config.scan_resolution**2
        elif gradient_type == "gradient_crossXgradient_cross":
            gradient = (signal_top - 2*signal_central + signal_bottom) / self.config.scan_resolution**2
        elif gradient_type == "gradient_coXgradient_cross":
            gradient = (np.roll(signal_top, -1) - np.roll(signal_bottom, -1) - np.roll(signal_top, 1) + np.roll(signal_bottom, 1)) / (4*self.config.scan_resolution**2)

            return gradient[1:-1]


    def get_nsamples(self):
        if self.config.sim_type == "signal":
            if self.config.beam_type == "pencil_beam":
                self.pad = 0
            else:
                self.pad = self.beam.del_beta.size/2
        else:
            self.pad = 1

        self.nsamples = int(self.config.t_segment*self.config.sampling_rate)*self.config.oversampling_rate + 2*self.pad


    #Calculating additional scan parameters
    def calculate_params(self):
        #Cross scan resolution
        self.config.theta_cross = 360.0*60.0*np.sin(np.radians(self.config.alpha))*self.config.t_spin/self.config.t_prec
        #Co scan resolution
        self.config.theta_co = 360*60*np.sin(np.radians(self.config.beta))/self.config.sampling_rate/self.config.t_spin

        self.config.scan_resolution = self.config.theta_co/self.config.oversampling_rate


    #Setting the initial axes
    def get_initial_axes(self):
        alpha = np.deg2rad(self.config.alpha)                                   #radians
        beta = np.deg2rad(self.config.beta)                                     #radians

        self.axis_spin = np.array([np.cos(alpha), 0.0, np.sin(alpha)])
        self.axis_prec = np.array([1.0, 0.0, 0.0])
        self.axis_rev = np.array([0.0, 0.0, 1.0])


    #Setting the initial positioning of pointing vectors
    def get_initial_vec(self, del_beta):
        alpha = np.deg2rad(self.config.alpha)                                   #radians
        beta = np.deg2rad(self.config.beta)                                     #radians
        if self.config.beam_type == "pencil":
            del_x = np.deg2rad(self.config.offset_x/60.0/60.0)    #radians
            del_y = np.deg2rad(self.config.offset_y/60.0/60.0)    #radians
        else:
            del_x = 0.0
            del_y = 0.0
        del_beta_rad = np.deg2rad(del_beta/60.0)                                #radians
        total_opening = alpha + beta + del_beta_rad + np.radians(self.config.focal_plane_del_beta/60.0)

        u_view = np.array([np.cos(total_opening), 0.0, np.sin(total_opening)])
        
        x_roll_axis = np.array([0.0, 1.0, 0.0])
        y_roll_axis = np.array([-np.sin(total_opening), 0.0, np.cos(total_opening)])

        q_x_roll = quaternion.make_quaternion(del_x, x_roll_axis)
        q_y_roll = quaternion.make_quaternion(del_y, y_roll_axis)
        q_offset = quaternion.multiply(q_x_roll, q_y_roll)
        u_view = quaternion.transform(q_offset, u_view)

        return u_view

    
#    @profile
    def generate_quaternion(self, segment):
        t_start = self.config.t_segment*segment

        t_steps = t_start + (1.0/self.config.sampling_rate/self.config.oversampling_rate)*np.arange(-self.pad, self.nsamples - self.pad)

        w_spin = 2*np.pi/self.config.t_spin
        w_prec = 2*np.pi/self.config.t_prec
        w_rev = 2*np.pi/self.config.t_year

        r_total = quaternion.multiply(quaternion.make_quaternion(w_rev*t_steps, self.axis_rev), quaternion.multiply(quaternion.make_quaternion(w_prec*t_steps, self.axis_prec), quaternion.make_quaternion(w_spin*t_steps, self.axis_spin)))

        return r_total


#    @profile
    #Simulating the pointing
    def get_vec_obv(self, v_init, rot_qt):
        v = quaternion.transform(rot_qt, v_init)

        if self.config.coordinate_system == "galactic":
            v = self.transform_to_gal_coords(v)

        return v


    def transform_to_gal_coords(self, v):
        rot = hp.Rotator(coord=['E', 'G'])
        theta, phi = hp.vec2ang(v)
        theta_gal, phi_gal = rot(theta, phi)
        return hp.ang2vec(theta_gal, phi_gal)


#    @profile
    def get_pol_ang(self, rot_qt, v_dir):
        alpha = np.deg2rad(self.config.alpha)                                   #radians
        beta = np.deg2rad(self.config.beta)                                     #radians
        total_opening = alpha + beta

        pol_ini = np.deg2rad(self.config.pol_phase_ini)
        pol_vec_ini = np.array([0.0, 1.0, 0.0])

        pol_vec = quaternion.transform(rot_qt, np.tile(pol_vec_ini, self.nsamples).reshape(-1,3))
        if self.config.coordinate_system == "galactic":
            pol_vec = self.transform_to_gal_coords(pol_vec)

        theta, phi = hp.vec2ang(v_dir)

        x_local = np.array(zip(np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi), -np.sin(theta)))
        y_local = np.array(zip(-np.sin(phi), np.cos(phi), np.zeros(phi.size)))

        proj_x = np.sum(pol_vec*x_local, axis=-1)
        proj_y = np.sum(pol_vec*y_local, axis=-1)

        #pol_ang = np.pi - (np.arctan2(proj_y, proj_x) + pol_ini) % np.pi 
        pol_ang = (np.arctan2(proj_y, proj_x) + pol_ini) % np.pi 
        del x_local, y_local, proj_x, proj_y

        return pol_ang 

    
    def add_noise(self):
        if self.config.noise_type == "white":
            return np.random.normal(scale=self.config.noise_sigma, size=self.nsamples - 2*self.pad)


    def get_sky_map(self):
        if self.config.sim_pol_type == "noise_only":
            self.sky_map = None
        elif self.config.sim_pol_type == "T":
            self.sky_map = hp.read_map(self.config.input_map, verbose=False)
        elif self.config.sim_pol_type == "QU":
            self.sky_map = hp.read_map(self.config.input_map, field=(0,1), verbose=False)
        elif self.config.sim_pol_type == "_QU":
            self.sky_map = hp.read_map(self.config.input_map, field=(1,2), verbose=False)
        else:
            self.sky_map = hp.read_map(self.config.input_map, field=(0,1,2), verbose=False)

        if not self.config.sim_pol_type == "noise_only":
            map_nside = hp.get_nside(self.sky_map)
            if map_nside != self.config.nside_in:
                prompt("NSIDE of config does not match NSIDE of map", sys.stdout)
                sys.exit()
        self.npix = hp.nside2npix(self.config.nside_in)


    def make_write_dir(self, segment):
        if not os.path.exists(self.bolo_dir):
            try:
                os.makedirs(self.bolo_dir)
            except OSError:
                pass

        segment_dir = self.get_segment_dir(segment) 
        #This will overwrite old data written for that particular segment
        if os.path.exists(segment_dir):
            if self.config.overwrite:
                shutil.rmtree(segment_dir)
                os.makedirs(segment_dir)
            else:
                pass
        else:
            os.makedirs(segment_dir)


    def set_bolo_dirs(self):
        self.sim_dir = os.path.join(self.config.general_output_dir, self.config.sim_tag)
        self.scan_dir = os.path.join(self.sim_dir, self.config.scan_tag)
        self.bolo_dir = os.path.join(self.scan_dir, self.config.name)


    def get_segment_dir(self, segment):
        segment_name = str(segment+1).zfill(4)
        return os.path.join(self.bolo_dir, segment_name)


    def write_timestream_data(self, ts_data, data_name, segment):
        write_dir = self.get_segment_dir(segment)
        if self.config.sim_type == "signal":
            if data_name == "signal":
                np.save(os.path.join(write_dir, data_name), ts_data[::self.config.oversampling_rate])
            elif data_name == "noise":
                np.save(os.path.join(write_dir, data_name), ts_data)
            else:
                if self.config.beam_type == "pencil":
                    np.save(os.path.join(write_dir, data_name), ts_data)
                else:
                    np.save(os.path.join(write_dir, data_name), ts_data[self.pad:-self.pad][::self.config.oversampling_rate])
        else:
            np.save(os.path.join(write_dir, data_name), ts_data)


    def display_params(self):
        display_string = "\n#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n"
        display_string += "#* SCAN PARAMETERS\n"
        display_string += "#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n"
        display_string += "Alpha : %f degrees\n" % (self.config.alpha)
        display_string += "Beta : %f degrees\n" % (self.config.beta)
        display_string += "T precession : %f hours\n" % (self.config.t_prec/60.0/60.0)
        display_string += "T spin : %f seconds\n" % (self.config.t_spin)
        display_string += "Sampling rate : %f Hz\n" % (self.config.sampling_rate)
        display_string += "Theta co : %f arcmin\n" % (self.config.theta_co)
        display_string += "Theta cross : %f arcmin\n" % (self.config.theta_cross)
        display_string += "Oversampling rate : %d\n" % (self.config.oversampling_rate)
        display_string += "Scan resolution for beam integration : %f arcmin\n" % (self.config.scan_resolution)
        if self.config.sim_pol_type != "noise_only":
            display_string += "Pixel size for NSIDE = %d : %f arcmin\n" % (self.config.nside_in, hp.nside2resol(self.config.nside_in, arcmin=True))
        n_steps = int(self.config.t_segment*self.config.sampling_rate)*self.config.oversampling_rate 
        display_string += "#Samples per segment : %d\n" %(n_steps)
        display_string += "#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n\n"
        prompt(display_string, sys.stdout)
