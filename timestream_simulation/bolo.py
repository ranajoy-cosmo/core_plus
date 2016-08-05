import numpy as np
import healpy as hp
import os
import shutil
from simulation.timestream_simulation.beam_kernel import Beam

class Bolo:

    def __init__(self, bolo_name, config):
        self.name = bolo_name
        bolo_config = importlib.import_module("simulation.timestream_simulation.bolo_config_files." + 
                config.bolo_config_file).bolos[bolo_name]
        self.config = Generic()
        self.config.__dict__.update(config.__dict__)
        self.config.__dict__.update(bolo_config.__dict__)
        self.calculate_parameters()
        self.beam = Beam(config, bolo_config)
        self.get_sky_map()
        self.get_initial_axes()
        self.get_nsamples()

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Simulating the time-ordered data for a given bolo with any beam
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    #@profile
    def simulate_timestream(self, segment):
        rot_qt = self.generate_quaternion(segment)

        #Simulating the scan along the centre of the FOV
        v_init = self.get_initial_vec(0.0)
        v = quaternion.transform(rot_qt, v_init)
        hit_pix_central, v = self.get_hitpix(v, ret_v=True)

        pol_ang = self.get_pol_ang(rot_qt, v) 
        if self.config.sim_type == "TQU":
            cos2 = np.cos(2*pol_ang)
            sin2 = np.sin(2*pol_ang)

        if "timestream_data" in self.config.timestream_data_products:
            write_dir = self.make_write_dir(segment)
            if self.config.do_pencil_beam:
                self.write_timestream_data(v, "pointing_vec", write_dir)
                self.write_timestream_data(pol_ang, "pol_ang", write_dir)
            else:
                self.write_timestream_data(v[self.pad:-self.pad][::self.config.oversampling_rate], "pointing_vec", write_dir)
                self.write_timestream_data(pol_ang[self.pad:-self.pad][::self.config.oversampling_rate], "pol_ang", write_dir)

        beam_kernel_row = self.beam.get_beam_row(0.0)                       #The input argument is the beam offset from the centre

        signal = self.get_signal(hit_pix_central, beam_kernel_row, cos2, sin2)

        for del_beta in self.beam.del_beta:
            if del_beta == 0.0:
                continue
            beam_kernel_row = self.beam.get_beam_row(del_beta)
            v_init = self.get_initial_vec(del_beta)
            v = quaternion.transform(rot_qt, v_init)
            hit_pix = self.get_hitpix(v)
            signal += self.get_signal(hit_pix, beam_kernel_row, cos2, sin2)

        beam_sum = np.sum(self.beam.beam_kernel)
        signal /= beam_sum


        if self.config.add_noise:
            signal += self.add_noise(nsamples) 

        if "timestream_data" in self.config.timestream_data_products:
            self.write_timestream_data(signal[::self.config.oversampling_rate], "signal", write_dir)

        if self.config.pipe_with_map_maker:
            return signal, v, pol_angle

        if self.config.make_scanned_map:
            hitmap = self.get_hitmap(hit_pix_central)
            return hitmap 

    
    def get_signal(self, hit_pix, beam_kernel_row, cos2, sin2):
        if self.config.pol_type == "noise_only":
            signal = np.zeros(self.nsamples - 2*self.pad) 

        elif self.config.pol_type == "T_only":
            signal = 0.5*self.sky_map[hit_pix]
            if not self.config.do_pencil_beam:
                signal = np.convolve(signal, beam_kernel_row, mode = 'valid')

        else:
            signal = 0.5*(sky_map[0][hitpix] + sky_map[1][hit_pix]*cos2 + sky_map[2][hitpix]*sin2)
            if not self.config.do_pencil_beam:
                signal = np.convolve(signal, beam_kernel_row, mode = 'valid')

        return signal

        

    def calculate_params(self):

        self.config.theta_cross = 360.0*60.0*np.sin(np.radians(self.config.alpha))*self.config.t_spin/self.config.t_prec
        self.config.theta_co = 360*60*np.sin(np.radians(self.config.beta))/self.config.sampling_rate/self.config.t_spin

        self.config.scan_resolution = self.config.theta_co/self.config.oversampling_rate


    def get_hitpix(self, v, ret_v=False):
        if self.config.gal_coords:
            rot = hp.Rotator(coord=['E', 'G'])
            theta, phi = hp.vec2ang(v)
            theta_gal, phi_gal = rot(theta, phi)
            hit_pix = hp.ang2pix(self.config.nside_in, theta_gal, phi_gal)
            if ret_v:
                v = hp.ang2vec(theta_gal, phi_gal)
        else:
            hit_pix = hp.vec2pix(self.config.nside, v[...,0], v[...,1], v[...,2])

        if ret_v:
            return hit_pix, v
        else:
            return hit_pix


    def get_nsamples(self):
        self.pad = self.beam.del_beta.size/2
        self.nsamples = int(self.config.t_segment*self.config.sampling_rate)*self.config.oversampling_rate + 2*self.pad

    
    def add_noise(self):
        if self.config.noise_type == "white":
            return np.random.normal(scale=self.config.noise_level, size=self.nsamples - 2*self.pad)


    #@profile
    def get_hitmap(self, hitpix):
        hitmap = np.bincount(hitpix, minlength=self.npix)

        return hitmap
    

    #@profile
    def generate_quaternion(self, segment):

        t_start = config.t_segment*segment

        t_steps = t_start + (1.0/config.sampling_rate/config.oversampling_rate)*np.arange(-self.pad, self.nsamples - self.pad)

        w_spin = -2*np.pi/config.t_spin
        w_prec = -2*np.pi/config.t_prec
        w_rev = -2*np.pi/config.t_year

        r_total = quaternion.multiply(quaternion.make_quaternion(w_rev*t_steps, self.axis_rev), quaternion.multiply(quaternion.make_quaternion(w_prec*t_steps, self.axis_prec), quaternion.make_quaternion(w_spin*t_steps, self.axis_spin)))

        return r_total

    #@profile
    def get_pol_ang(self, rot_qt, v_dir=None):

        pol_init = np.deg2rad(self.bolo_params.pol_phase_ini)
        x_axis = np.array([0.0, 1.0, 0.0])

        pol_vec = quaternion.transform(rot_qt, np.tile(x_axis, n_steps).reshape(-1,3))

        if v_dir is None:
            v_init = self.get_initial_vec(0.0)
            v_dir = quaternion.transform(rot_qt, v_init)

        theta, phi = hp.vec2ang(v_dir)

        x_local = np.array(zip(np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi), -np.sin(theta)))
        y_local = np.array(zip(-np.sin(phi), np.cos(phi), np.zeros(phi.size)))

        proj_x = np.sum(pol_vec*x_local, axis=-1)
        proj_y = np.sum(pol_vec*y_local, axis=-1)

        pol_ang = (np.arctan2(proj_y, proj_x) + pol_init) % np.pi 
        #pol_ang *= 2.0

        return pol_ang 

    def get_initial_vec(self, del_beta):
        alpha = np.deg2rad(self.config.alpha)                                   #radians
        beta = np.deg2rad(self.config.beta)                                     #radians
        if self.config.do_pencil_beam:
            del_x = np.deg2rad(self.bolo_params.pointing_offset_x/60.0/60.0)    #radians
            del_y = np.deg2rad(self.bolo_params.pointing_offset_y/60.0/60.0)    #radians
        else:
            del_x = 0.0
            del_y = 0.0
        del_beta_rad = np.deg2rad(del_beta/60.0)                                #radians

        self.u_view = np.array([np.cos(alpha + beta + del_beta_rad), 0.0, np.sin(alpha + beta + del_beta_rad)])


    def get_initial_axes(self):
        alpha = np.deg2rad(self.config.alpha)                                   #radians
        beta = np.deg2rad(self.config.beta)                                     #radians

        self.axis_spin = np.array([np.cos(alpha), 0.0, np.sin(alpha)])
        self.axis_prec = np.array([1.0, 0.0, 0.0])
        self.axis_rev = np.array([0.0, 0.0, 1.0])


    def get_sky_map(self):
        if self.config.sim_pol_type == "noise_only":
            self.sky_map = None
        elif self.config.sim_pol_type == "T_only":
            self.sky_map = hp.read_map(self.config.map_file)
        else:
            self.sky_map = hp.read_map(self.config.map_file, field=(0,1,2))

        self.npix = hp.nside2npix(self.config.nside_in)


    def make_write_dir(self, segment):
        sim_dir = os.path.join(self.config.general_data_dir, self.config.sim_tag)
        if not os.path.exists(sim_dir):
            os.makedirs(sim_dir)
        scan_dir = os.path.join(sim_dir, self.config.scan_tag)
        if not os.path.exists(scan_dir):
            os.makedirs(scan_dir)

        bolo_dir = os.path.join(sim_dir, self.config.scan_tag, self.bolo_name)
        if not os.path.exists(bolo_dir):
            os.makedirs(bolo_dir)

        segment_name = str(segment+1).zfill(4)
        segment_dir = os.path.join(bolo_dir, segment_name)
        if os.path.exists(segment_dir):
            shutil.rmtree(segment_dir)

        os.path.makedirs(segment_dir)

        return segment_dir


    def write_timestream_data(self, ts_data, data_name, write_dir):
        np.save(os.path.join(write_dir, data_name), ts_data)
