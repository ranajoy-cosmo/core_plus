import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import os
import simulation.power_spectrum.noise_power_spectra as nps
import simulation.lib.plotting.plot_th as pt
import simulation.power_spectrum.spectra_tools as st

class Plot_Maker():

    def __init__(self, map_dir, plot_dir, config):
        self.map_dir = map_dir
        self.plot_dir = plot_dir
        self.config = config
        if not os.path.exists(self.plot_dir):
            os.makedirs(self.plot_dir)
        else:
            pass
        
        self.load_maps_and_spectra()
        self.get_spectra()
        self.rescale_maps_and_spectra()


    def load_maps_and_spectra(self):
        self.noise_map = np.array(hp.read_map(os.path.join(self.map_dir, "sky_map.fits"), field=(0,1,2)))
        self.cov_map = np.array(hp.read_map(os.path.join(self.map_dir, "covariance_maps.fits"), field=(0,1,2,3,4,5)))
        self.hitmap = hp.read_map(os.path.join(self.map_dir, "hitmap.fits"))
        self.nside = hp.get_nside(self.noise_map)
        self.f_sky = float(np.sum(self.hitmap > 3)) / self.hitmap.size


    def get_spectra(self):
        if self.config['new_spectra']:
            self.spectra_obv = st.estimate_cl(self.noise_map, lmax = 3*self.nside - 1, fwhm=self.config['beam_fwhm'])
            np.save(os.path.join(self.map_dir, "spectra_noise"), self.spectra_obv)
        else:
            self.spectra_obv = np.load(os.path.join(self.map_dir, "spectra_noise.npy"))
        self.spectra_th = nps.sensitivity_to_noise_mean_Cl(self.config['det_sens'], self.config['beam_fwhm'], self.config['lmax_plot'], self.config['t_mission'], self.config['n_det_rescale'], self.f_sky)


    def rescale_maps_and_spectra(self):
        pix_area = hp.nside2resol(self.nside, arcmin=True)**2
        n_det_ratio = float(self.config['n_det'])/self.config['n_det_rescale']
        noise_rms = self.config['det_sens'] * np.sqrt(self.config['scan_freq'])
        self.noise_map *= np.sqrt(pix_area) * np.sqrt(n_det_ratio)
        self.cov_map *= 0.25 * noise_rms**2 * pix_area * n_det_ratio
        self.spectra_obv *= n_det_ratio


    def plot_noise_spectra(self):
        lmax = self.config['lmax_plot']
        pt.plot_theoretical(lmax)
        pt.plot_spectra(self.spectra_th[0], lmax=lmax, color='b', label="TT : Theoretical")
        pt.plot_spectra(self.spectra_th[1], lmax=lmax, color='g', label="EE/BB : Theoretical")
        pt.plot_spectra(self.spectra_obv[0], lmax=lmax, color='r', label="TT : Observed")
        pt.plot_spectra(self.spectra_obv[1], lmax=lmax, color='y', label="EE : Observed")
        pt.plot_spectra(self.spectra_obv[2], lmax=lmax, color='c', label="BB : Observed")
        pt.make_decorations(ylim=[1e-3,100])
        plt.title("Noise spectra : " + self.config['title'])
        plt.savefig(os.path.join(self.plot_dir, "noise_spectra.png"))


    def plot_noise_histogram(self):
        weights = np.ones(12*self.nside**2) / (12*self.nside**2)
        range = self.config['noise_hist_range']
        plt.hist(self.noise_map[0], bins=100, range=range, label="T", alpha=0.5, weights=weights)
        plt.hist(self.noise_map[1], bins=100, range=range, label="Q", alpha=0.5, weights=weights)
        plt.hist(self.noise_map[2], bins=100, range=range, label="U", alpha=0.5, weights=weights)
        plt.xlabel("$[\mu K.arcmin]$", fontsize=15)
        plt.ylabel("$f_{sky}$", fontsize=15)
        plt.legend()
        plt.title("Noise map histogram : " + self.config['title'])
        plt.savefig(os.path.join(self.plot_dir, "noise_map_histogram.png"))


    def plot_noise_histogram_cumulative(self):
        weights = np.ones(12*self.nside**2) / (12*self.nside**2)
        range = self.config['noise_hist_range']
        range = (0, range[1])
        plt.hist(np.abs(self.noise_map[0]), bins=100, range=range, label="T", alpha=0.5, cumulative=True, histtype='step', weights=weights)
        plt.hist(np.abs(self.noise_map[1]), bins=100, range=range, label="Q", alpha=0.5, cumulative=True, histtype='step', weights=weights)
        plt.hist(np.abs(self.noise_map[2]), bins=100, range=range, label="U", alpha=0.5, cumulative=True, histtype='step', weights=weights)
        plt.xlabel("$|\sigma| \, [\mu K.arcmin]$", fontsize=15)
        plt.ylabel("$f_{sky} \,<\, |\sigma|$", fontsize=15)
        plt.legend(loc="upper left")
        plt.title("Noise map histogram : " + self.config['title'])
        plt.savefig(os.path.join(self.plot_dir, "cumulative_noise_map_histogram.png"))


    def plot_cov_diagonal(self):
        weights = np.ones(12*self.nside**2) / (12*self.nside**2)
        range = self.config['cov_diag_range']
        bins = self.config['cov_diag_bins']
        bin_width = float(range[1] - range[0]) / bins
        plt.hist(np.sqrt(self.cov_map[0]), bins=bins, range=range, label="TT", alpha=0.5, weights=weights)
        plt.hist(np.sqrt(self.cov_map[3]), bins=bins, range=range, label="QQ", alpha=0.5, weights=weights)
        plt.hist(np.sqrt(self.cov_map[5]), bins=bins, range=range, label="UU", alpha=0.5, weights=weights)
        plt.xlabel("$\sigma \, [\mu K.arcmin]$", fontsize=15)
        plt.ylabel("$f_{sky}$", fontsize=15)
        plt.text(12, 0.095, r'$0.1 \mu K.arcmin \, bins$', fontsize=15)
        plt.legend(loc="upper left")
        plt.title(self.config['title'].format("Polarisation Sensitivity"))
        plt.savefig(os.path.join(self.plot_dir, "cov_diagonal_histogram.png"))


    def plot_cov_diagonal_cumulative(self):
        weights = np.ones(12*self.nside**2) / (12*self.nside**2)
        range = self.config['cov_diag_range']
        plt.hist(np.sqrt(self.cov_map[0]), bins=100, range=range, label="TT", alpha=0.5, cumulative=True, histtype='step', weights=weights)
        plt.hist(np.sqrt(self.cov_map[3]), bins=100, range=range, label="QQ", alpha=0.5, cumulative=True, histtype='step', weights=weights)
        plt.hist(np.sqrt(self.cov_map[5]), bins=100, range=range, label="UU", alpha=0.5, cumulative=True, histtype='step', weights=weights)
        plt.xlabel("$\sigma \, [\mu K.arcmin]$", fontsize=15)
        plt.ylabel("$f_{sky} \,<\, \sigma$", fontsize=15)
        plt.legend(loc="upper left")
        plt.title(self.config['title'].format("Polarisation Sensitivity"))
        plt.savefig(os.path.join(self.plot_dir, "cov_diagonal_cumulative_histogram.png"))


    def plot_cov_off_diagonal(self):
        weights = np.ones(12*self.nside**2) / (12*self.nside**2)
        range = self.config['cov_off_diag_range']
        plt.hist(np.abs(np.sqrt(self.cov_map[1] + 0j)), bins=100, range=range, label="TQ", alpha=0.5, weights=weights)
        plt.hist(np.abs(np.sqrt(self.cov_map[2] + 0j)), bins=100, range=range, label="TU", alpha=0.5, weights=weights)
        plt.hist(np.abs(np.sqrt(self.cov_map[4] + 0j)), bins=100, range=range, label="QU", alpha=0.5, weights=weights)
#        plt.hist(np.sqrt(self.cov_map[1]), bins=100, range=range, label="TQ", alpha=0.5, weights=weights)
#        plt.hist(np.sqrt(self.cov_map[2]), bins=100, range=range, label="TU", alpha=0.5, weights=weights)
#        plt.hist(np.sqrt(self.cov_map[4]), bins=100, range=range, label="QU", alpha=0.5, weights=weights)
        plt.xlabel("$|\sigma| \, [\mu K.arcmin]$", fontsize=15)
        plt.ylabel("$f_{sky}$", fontsize=15)
        plt.legend()
        plt.title("Polarisation Sensitivity : " + self.config['title'])
        plt.savefig(os.path.join(self.plot_dir, "cov_off_diagonal_histogram.png"))

    def plot_cov_off_diagonal_cumulative(self):
        weights = np.ones(12*self.nside**2) / (12*self.nside**2)
        range = self.config['cov_off_diag_range']
        plt.hist(np.abs(np.sqrt(self.cov_map[1] + 0j)), bins=100, range=range, label="TQ", alpha=0.5, cumulative=True, histtype='step', weights=weights)
        plt.hist(np.abs(np.sqrt(self.cov_map[2] + 0j)), bins=100, range=range, label="TU", alpha=0.5, cumulative=True, histtype='step', weights=weights)
        plt.hist(np.abs(np.sqrt(self.cov_map[4] + 0j)), bins=100, range=range, label="QU", alpha=0.5, cumulative=True, histtype='step', weights=weights)
#        plt.hist(np.sqrt(self.cov_map[1]), bins=100, range=range, label="TQ", alpha=0.5, cumulative=True, histtype='step', weights=weights)
#        plt.hist(np.sqrt(self.cov_map[2]), bins=100, range=range, label="TU", alpha=0.5, cumulative=True, histtype='step', weights=weights)
#        plt.hist(np.sqrt(self.cov_map[4]), bins=100, range=range, label="QU", alpha=0.5, cumulative=True, histtype='step', weights=weights)
        plt.xlabel("$|\sigma| \, [\mu K.arcmin]$", fontsize=15)
        plt.ylabel("$f_{sky} \,<\, |\sigma|$", fontsize=15)
        plt.legend(loc="upper left")
        plt.title("Polarisation Sensitivity : " + self.config['title'])
        plt.savefig(os.path.join(self.plot_dir, "cov_off_diagonal_cumulative_histogram.png"))


    def plot_cov_maps(self, plot_list, range=None):
        label_dict = {'TT' : 0, 'QQ' : 3, 'UU' : 5, 'TQ' : 1, 'TU' : 2, 'QU' : 4}
        for plot_name in list(set(label_dict.keys()) & set(plot_list)):
            if range:
                hp.mollview(np.abs(np.sqrt(self.cov_map[label_dict[plot_name]] + 0j)), min=range[0], max=range[1], unit="$\mu K.arcmin$")
            else:
                hp.mollview(np.abs(np.sqrt(self.cov_map[label_dict[plot_name]] + 0j)), unit="$\mu K.arcmin$")
#            hp.mollview(np.abs(np.sqrt(self.cov_map[label_dict[plot_name]] + 0j)), unit="$\mu K.arcmin$", title=title)
            plt.savefig(os.path.join(self.plot_dir, "cov_" + plot_name +".png"))



t_year = 365.25*24*60*60

map_dir = os.path.join("/global/homes/b/banerji/simulation/output/focal_plane_noise", "36_bolo_4_year")
plot_dir = os.path.join("/global/homes/b/banerji/simulation/output/focal_plane_noise/mission_paper_plots_2", "36_det_4_years")

config = {}
config['title'] = "$145$ GHz {} : $144$ Detectors : $4$ years"
config['t_mission'] = 4*t_year
config['n_det'] = 36.0
config['n_det_rescale'] = 144.0
config['det_sens'] = 39.9
config['scan_freq'] = 85.0
config['beam_fwhm'] = np.radians(7.7/60.0)
config['lmax_plot'] = 2000
config['noise_hist_range'] = (-20,20)
config['noise_hist_bins'] = 100
config['cov_diag_range'] = (1,7)
config['cov_diag_bins'] = 120
config['new_spectra'] = False
