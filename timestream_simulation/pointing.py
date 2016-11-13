import numpy as np
import healpy as hp
import sys
from simulation.lib.utilities.generic_class import Generic
import simulation.lib.utilities.prompter as prompter
import simulation.lib.quaternion.quaternion as quaternion

class Pointing():
    
    def __init__(self, config, segment, beam_width):
        self.config = Generic()
        self.config.__dict__.update(config.__dict__)
        self.set_sim_time_and_samples(segment, beam_width)

    def set_sim_time_and_samples(self, segment, beam_width):
        self.t_start = segment*self.config.t_segment
        self.t_stop = (segment + 1)*self.config.t_segment
        self.delta_t = 1.0/self.config.sampling_rate/self.config.oversampling_rate
        self.pad = beam_width/2
        self.n_samples = self.config.t_segment*self.config.sampling_rate*self.config.oversampling_rate + 2*self.pad

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

        t_steps = t_start + (1.0/self.config.sampling_rate/self.config.oversampling_rate)*np.arange(-self.pad, self.nsamples - self.pad)

        w_spin = -2*np.pi/self.config.t_spin
        w_prec = -2*np.pi/self.config.t_prec
        w_rev = -2*np.pi/self.config.t_year

        r_total = quaternion.multiply(quaternion.make_quaternion(w_rev*t_steps, self.axis_rev), quaternion.multiply(quaternion.make_quaternion(w_prec*t_steps, self.axis_prec), quaternion.make_quaternion(w_spin*t_steps, self.axis_spin)))
        #r_total = quaternion.multiply(quaternion.make_quaternion(w_prec*t_steps, self.axis_prec), quaternion.make_quaternion(w_spin*t_steps, self.axis_spin))

        return r_total


