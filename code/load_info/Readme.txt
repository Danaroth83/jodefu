How to add new datasets:
In order:
load_cut: Include your image tag and the cut you desire for the image
load_default_sensor: Assign a sensor to your image tag
load_Bands_to_sharpen: Save data on the bands to select and display
load_resolution: Info on the bits and if available on GSD
load_wavelength: Assign central wavelengths and bandwidth to your sensor
load_MTF and load_MTF_PAN: saving the gain at Nyquist frequency (if available)

Add Spectral Responses in:
../../data/relative_spectral_responses/[sensorname]_Spectral_Responses.mat
If not available, change manually:
load_minmaxspectralresponses
load_spectralweights


For all functions that rely on spectra, it is required to have their 
Spectral Responses in the folder:

../../data/relative_spectral_responses/[sensorname]_Spectral_Responses.mat

where [sensorname] is the name of the sensor (such as ALI, QB, WV3, etc.)

They shold contain the variables:

wavelength_nm: the wavelengths in nanometers (size: Nw x 1)
Spectral_Responses_Matrix: The spectral responses for each band (size: Nb x Nw)

where Nw is the number of the samples for which we discretize the wavelenghts axis
      Nb is the number of bands associated to the sensor
