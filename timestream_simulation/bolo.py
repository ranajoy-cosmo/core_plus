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
import simulation.lib.utilities.prompter as prompter
import simulation.lib.quaternion.quaternion as quaternion

class Bolo:

    def __init__(self, bolo_name, config):
        self.name = bolo_name
        self.config = Generic()
        bolo_config = importlib.import_module(config.bolo_config_file).bolo_config.bolos[bolo_name]
        self.config.__dict__.update(config.__dict__)
        self.config.__dict__.update(bolo_config.__dict__)
        if config.simulate_ts:
            self.beam = Beam(self.config, bolo_config)
            self.noise_class = Noise(self.config)
            if not config.sim_pol_type == "noise_only":
                self.get_sky_map()
            self.calculate_params()
            self.get_initial_axes()
            self.get_nsamples()
        self.set_bolo_dirs()

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Simulating the time-ordered data for a given bolo with any beam
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    #@profile
    def simulate_timestream(self, segment, return_field=["signal", "v", "pol_ang"]):
        if segment == 0:
            #self.display_params()
            #self.beam.display_beam_settings()
            if self.config.write_beam:
                self.beam.write_beam(self.bolo_dir)
        rot_qt = self.generate_quaternion(segment)

        t_stream = {"signal" : None, "v" : None, "pol_ang" : None}

        prompter.prompt("0.0")
        #Simulating the scan along the centre of the FOV
        v_init = self.get_initial_vec(0.0)
        v_central = self.get_v_obv(v_init, rot_qt)

        pol_ang = self.get_pol_ang(rot_qt, v_central) 
        if self.config.sim_pol_type == "T":
            cos2=None
            sin2=None
        else:
            cos2 = np.cos(2*pol_ang)
            sin2 = np.sin(2*pol_ang)
 
        if "timestream_data" in self.config.timestream_data_products:
            self.make_write_dir(segment)
            self.write_timestream_data(v_central, "pointing_vec", segment)
            self.write_timestream_data(pol_ang, "pol_ang", segment)
        
        if not self.config.pipe_with_map_maker:
            del pol_ang

        beam_kernel_row = self.beam.get_beam_row(0.0)                       #The input argument is the beam offset from the centre
        hit_pix = hp.vec2pix(self.config.nside_in, v_central[...,0], v_central[...,1], v_central[...,2])
        signal = self.get_signal(hit_pix, beam_kernel_row, cos2, sin2)

        for del_beta in self.beam.del_beta:
            if del_beta == 0.0:
                continue
            prompter.prompt(str(del_beta))
            beam_kernel_row = self.beam.get_beam_row(del_beta)
            v_init = self.get_initial_vec(del_beta)
            v = quaternion.transform(rot_qt, v_init)
            hit_pix = hp.vec2pix(self.config.nside_in, v[...,0], v[...,1], v[...,2])
            signal += self.get_signal(hit_pix, beam_kernel_row, cos2, sin2)

        beam_sum = np.sum(self.beam.beam_kernel[0])
        signal /= beam_sum

        if self.config.add_noise:
            noise = self.noise_class.simulate_timestream_noise_from_parameters()
            if "noise" in self.config.timestream_data_products:
                self.write_timestream_data(noise, "noise", segment)
            signal[::self.config.oversampling_rate] += noise 

        if "timestream_data" in self.config.timestream_data_products:
            self.write_timestream_data(signal, "signal", segment)

        if self.config.pipe_with_map_maker:
            if self.config.do_pencil_beam:
                t_stream["signal"] = signal
                t_stream["v"] = v_central
                t_stream["pol_ang"] = pol_ang
                return t_stream
            else:
                t_stream["signal"] = signal[::self.config.oversampling_rate]
                t_stream["v"] = v_central[self.pad:-self.pad][::self.config.oversampling_rate]
                t_stream["pol_ang"] = pol_ang[self.pad:-self.pad][::self.config.oversampling_rate]
                return t_stream

        if "hitmap" in self.config.timestream_data_products:
            del signal
            hit_pix = hp.vec2pix(self.config.nside_in, v_central[...,0], v_central[...,1], v_central[...,2])
            hitmap = self.get_hitmap(hit_pix)
            return hitmap 
    

    def get_signal(self, hit_pix, beam_kernel_row, cos2, sin2):
        if self.config.sim_pol_type == "noise_only":
            signal = np.zeros(self.nsamples - 2*self.pad) 

        elif self.config.sim_pol_type == "T":
            if self.config.do_pencil_beam:
                signal = 0.5*self.sky_map[hit_pix]
            else: 
                signal = np.convolve(0.5*self.sky_map[hit_pix], beam_kernel_row[0], mode='valid')

        elif self.config.sim_pol_type in ["QU", "_QU"]:
            if self.config.do_pencil_beam:
                signal = 0.5*(self.sky_map[0][hit_pix]*cos2 + self.sky_map[1][hit_pix]*sin2) 
            else:
                signal = np.convolve(0.5*(self.sky_map[0][hit_pix]*cos2 + self.sky_map[1][hit_pix]*sin2), -1.0*beam_kernel_row[1], mode='valid')
                signal += np.convolve(0.5*(self.sky_map[0][hit_pix]*sin2 + self.sky_map[1][hit_pix]*cos2), -1.0*beam_kernel_row[2], mode='valid')

        elif self.config.sim_pol_type == "pol_only":
            if self.config.do_pencil_beam:
                signal = 0.5*(self.sky_map[1][hit_pix]*cos2 + self.sky_map[2][hit_pix]*sin2)
            else:
                signal = np.convolve(0.5*(self.sky_map[1][hit_pix]*cos2 + self.sky_map[2][hit_pix]*sin2), -1.0*beam_kernel_row[1], mode='valid')
                signal += np.convolve(0.5*(self.sky_map[1][hit_pix]*sin2 + self.sky_map[2][hit_pix]*cos2), -1.0*beam_kernel_row[2], mode='valid')

        else:
            if self.config.do_pencil_beam:
                signal = 0.5*(self.sky_map[0][hit_pix] + self.sky_map[1][hit_pix]*cos2 + self.sky_map[2][hit_pix]*sin2) 
            else:
                signal = np.convolve(0.5*self.sky_map[0][hit_pix], beam_kernel_row[0], mode='valid')
                signal += np.convolve(0.5*(self.sky_map[1][hit_pix]*cos2 + self.sky_map[2][hit_pix]*sin2), -1.0*beam_kernel_row[1], mode='valid')
                signal += np.convolve(0.5*(self.sky_map[1][hit_pix]*sin2 + self.sky_map[2][hit_pix]*cos2), -1.0*beam_kernel_row[2], mode='valid')

        return signal


    def get_nsamples(self):
        self.pad = self.beam.del_beta.size/2
        self.nsamples = int(self.config.t_segment*self.config.sampling_rate)*self.config.oversampling_rate + 2*self.pad


    def calculate_params(self):

        self.config.theta_cross = 360.0*60.0*np.sin(np.radians(self.config.alpha))*self.config.t_spin/self.config.t_prec
        self.config.theta_co = 360*60*np.sin(np.radians(self.config.beta))/self.config.sampling_rate/self.config.t_spin

        self.config.scan_resolution = self.config.theta_co/self.config.oversampling_rate


    def get_initial_axes(self):
        alpha = np.deg2rad(self.config.alpha)                                   #radians
        beta = np.deg2rad(self.config.beta)                                     #radians

        self.axis_spin = np.array([np.cos(alpha), 0.0, np.sin(alpha)])
        self.axis_prec = np.array([1.0, 0.0, 0.0])
        self.axis_rev = np.array([0.0, 0.0, 1.0])


    def get_initial_vec(self, del_beta):
        alpha = np.deg2rad(self.config.alpha)                                   #radians
        beta = np.deg2rad(self.config.beta)                                     #radians
        if self.config.do_pencil_beam:
            del_x = np.deg2rad(self.config.offset_x/60.0/60.0)    #radians
            del_y = np.deg2rad(self.config.offset_y/60.0/60.0)    #radians
        else:
            del_x = 0.0
            del_y = 0.0
        del_beta_rad = np.deg2rad(del_beta/60.0)                                #radians
        total_opening = alpha + beta + del_beta_rad

        u_view = np.array([np.cos(total_opening), 0.0, np.sin(total_opening)])
        
        x_roll_axis = np.array([0.0, 1.0, 0.0])
        y_roll_axis = np.array([-np.sin(total_opening), 0.0, np.cos(total_opening)])

        q_x_roll = quaternion.make_quaternion(del_x, x_roll_axis)
        q_y_roll = quaternion.make_quaternion(del_y, y_roll_axis)
        q_offset = quaternion.multiply(q_x_roll, q_y_roll)
        u_view = quaternion.transform(q_offset, u_view)

        return u_view


    def generate_quaternion(self, segment):

        t_start = self.config.t_segment*segment

        t_steps = t_start + (1.0/self.config.sampling_rate/self.config.oversampling_rate)*np.arange(-self.pad, self.nsamples - self.pad)

        w_spin = -2*np.pi/self.config.t_spin
        w_prec = -2*np.pi/self.config.t_prec
        w_rev = -2*np.pi/self.config.t_year

        r_total = quaternion.multiply(quaternion.make_quaternion(w_rev*t_steps, self.axis_rev), quaternion.multiply(quaternion.make_quaternion(w_prec*t_steps, self.axis_prec), quaternion.make_quaternion(w_spin*t_steps, self.axis_spin)))
        #r_total = quaternion.multiply(quaternion.make_quaternion(w_prec*t_steps, self.axis_prec), quaternion.make_quaternion(w_spin*t_steps, self.axis_spin))

        return r_total


    def get_v_obv(self, v_init, rot_qt):
        v = quaternion.transform(rot_qt, v_init)

        if self.config.gal_coords:
            v = self.transform_to_gal_coords(v)

        return v


    def transform_to_gal_coords(self, v):
        rot = hp.Rotator(coord=['E', 'G'])
        theta, phi = hp.vec2ang(v)
        theta_gal, phi_gal = rot(theta, phi)
        return hp.ang2vec(theta_gal, phi_gal)


    def get_pol_ang(self, rot_qt, v_dir):
        alpha = np.deg2rad(self.config.alpha)                                   #radians
        beta = np.deg2rad(self.config.beta)                                     #radians
        total_opening = alpha + beta

        pol_ini = np.deg2rad(self.config.pol_phase_ini)
        pol_vec_ini = np.array([0.0, 1.0, 0.0])

        pol_vec = quaternion.transform(rot_qt, np.tile(pol_vec_ini, self.nsamples).reshape(-1,3))
        if self.config.gal_coords:
            pol_vec = self.transform_to_gal_coords(pol_vec)

        theta, phi = hp.vec2ang(v_dir)

        x_local = np.array(zip(np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi), -np.sin(theta)))
        y_local = np.array(zip(-np.sin(phi), np.cos(phi), np.zeros(phi.size)))

        proj_x = np.sum(pol_vec*x_local, axis=-1)
        proj_y = np.sum(pol_vec*y_local, axis=-1)

        pol_ang = np.pi - (np.arctan2(proj_y, proj_x) + pol_ini) % np.pi 

        return pol_ang 

    
    def add_noise(self):
        if self.config.noise_type == "white":
            return np.random.normal(scale=self.config.noise_sigma, size=self.nsamples - 2*self.pad)


    def get_hitmap(self, hitpix):
        hitmap = np.bincount(hitpix, minlength=self.npix)

        return hitmap
    

    def get_sky_map(self):
        if self.config.sim_pol_type == "noise_only":
            self.sky_map = None
        elif self.config.sim_pol_type == "T":
            self.sky_map = hp.read_map(self.config.input_map)
        elif self.config.sim_pol_type == "QU":
            self.sky_map = hp.read_map(self.config.input_map, field=(0,1))
        elif self.config.sim_pol_type == "_QU":
            self.sky_map = hp.read_map(self.config.input_map, field=(1,2))
        else:
            self.sky_map = hp.read_map(self.config.input_map, field=(0,1,2))

        if not self.config.sim_pol_type == "noise_only":
            map_nside = hp.get_nside(self.sky_map)
            if map_nside != self.config.nside_in:
                prompter.prompt_warning("NSIDE of config does not match NSIDE of map")
                sys.exit()
        self.npix = hp.nside2npix(self.config.nside_in)


    def make_write_dir(self, segment):
        if not os.path.exists(self.bolo_dir):
            os.makedirs(self.bolo_dir)

        segment_dir = self.get_segment_dir(segment) 
        if os.path.exists(segment_dir):
            shutil.rmtree(segment_dir)

        os.makedirs(segment_dir)


    def set_bolo_dirs(self):
        self.sim_dir = os.path.join(self.config.general_data_dir, self.config.sim_tag)
        self.scan_dir = os.path.join(self.sim_dir, self.config.scan_tag)
        self.bolo_dir = os.path.join(self.scan_dir, self.name)


    def get_segment_dir(self, segment):
        segment_name = str(segment+1).zfill(4)
        return os.path.join(self.bolo_dir, segment_name)


    def write_timestream_data(self, ts_data, data_name, segment):
        write_dir = self.get_segment_dir(segment)
        if data_name == "signal":
            np.save(os.path.join(write_dir, data_name), ts_data[::self.config.oversampling_rate])
        elif data_name == "noise":
            np.save(os.path.join(write_dir, data_name), ts_data)
        else:
            if self.config.do_pencil_beam:
                np.save(os.path.join(write_dir, data_name), ts_data)
            else:
                np.save(os.path.join(write_dir, data_name), ts_data[self.pad:-self.pad][::self.config.oversampling_rate])


    def read_timestream(self, segment, noise_only=False, read_list=["signal", "v", "pol_ang"]):
        segment_dir = self.get_segment_dir(segment)
        t_stream = {"signal" : None, "v" : None, "pol_ang" : None}  
        if "signal" in read_list:
            if noise_only:
                t_stream["signal"] = np.load(os.path.join(segment_dir, "noise.npy"))
            else:
                t_stream["signal"] = np.load(os.path.join(segment_dir, "signal.npy"))
        if "v" in read_list:
            t_stream["v"] = np.load(os.path.join(segment_dir, "pointing_vec.npy"))
        if "pol_ang" in read_list:
            t_stream["pol_ang"] = np.load(os.path.join(segment_dir, "pol_ang.npy"))

        return t_stream 


    def display_params(self):
        display_string = ""
        display_string += "Alpha : %f degrees\n" % (self.config.alpha)
        display_string += "Beta : %f degrees\n" % (self.config.beta)
        t_flight = self.config.t_segment*len(self.config.segment_list)
        display_string += "T flight : %f hours / %f days\n" % (t_flight/60.0/60.0, t_flight/60.0/60.0/24.0)
        display_string += "T segment : %f hours / %f days\n" % (self.config.t_segment/60.0/60.0, self.config.t_segment/60.0/60.0/24)
        display_string += "T precession : %f hours\n" % (self.config.t_prec/60.0/60.0)
        display_string += "T spin : %f seconds\n" % (self.config.t_spin)
        display_string += "Scan sampling rate : %f Hz\n" % (self.config.sampling_rate)
        display_string += "Theta co : %f arcmin\n" % (self.config.theta_co)
        display_string += "Theta cross : %f arcmin\n" % (self.config.theta_cross)
        display_string += "Oversampling rate : %d\n" % (self.config.oversampling_rate)
        display_string += "Scan resolution for beam integration : %f arcmin\n" % (self.config.scan_resolution)
        display_string += "Pixel size for NSIDE = %d : %f arcmin\n" % (self.config.nside_in, hp.nside2resol(self.config.nside_in, arcmin=True))
        n_steps = int(self.config.t_segment*self.config.sampling_rate)*self.config.oversampling_rate 
        display_string += "#Samples per segment : %d\n" %(n_steps)

        prompter.prompt(display_string, True)
        #prompter.prompt(display_string, False if not self.config.action=="display_params" else True)

