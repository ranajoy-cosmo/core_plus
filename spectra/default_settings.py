from simulation.lib.utilities.generic_class import Generic
from simulation.settings.global_settings import global_scanning, global_paths
import os

settings = Generic()

settings.tag = global_paths.tag

settings.lmax = 4000

settings.cosmo_params = {}
settings.cosmo_params['source'] = "pycamb_default"

settings.normalise_spectra = False

settings.make_alm = True

out_folder = global_paths.spectra_folder
file_name = settings.tag + '_' + str(settings.lmax)
settings.spectra_file = os.path.join(out_folder, file_name)
settings.write_spectra = True
settings.plot_spectra = True
settings.return_spectra = False
