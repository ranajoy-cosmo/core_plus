from simulation.lib.utilities.generic_class import Generic

settings = Generic()

settings.tag = 'test'

settings.lmax = 3000

settings.cosmo_params = {}
settings.cosmo_params['source'] = "pycamb_default"

settings.normalise_spectra = False

out_folder = "/Users/banerji/CORE+/simulation/maps_and_spectra/spectra/"
file_name = settings.tag + '_' + str(settings.lmax)
settings.spectra_file = out_folder + file_name
settings.write_spectra = True
settings.plot_spectra = False
settings.return_spectra = False
