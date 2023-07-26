# Module to calculate noise using the Pyxel tool

import pyxel
from pyxel.util import fit_into_array
import numpy as np
import time
import numpy.typing as npt

def calc_noise_pyxel( 
        input_photon: npt.NDArray, 
        yaml_file: str,
    ) -> np.ndarray:
    """ Takes photon array and yaml parameter file name and processes it to add noise

    Parameters
    ----------
    input_photon: ndarray (pass photon array)
    yaml_file: yaml file containing detector properties

    Returns
    -------
    output_image: ndarray

    """
    t1=time.time()
    
    # Add pyxel noise calculation stuff here
    # Make sure that the photon only has positive elements
    input_photon[input_photon<0.0] = 0.0
    row,col=np.shape(input_photon)

    # Get Pyxel configs from the yaml file

    config = pyxel.load(yaml_file=yaml_file)

    # Set pyxel stuff from the config file
    exposure = config.exposure  # class Single
    pipeline = config.pipeline  # class DetectionPipeline

    # Currently this only supports ccd and apd detector, 
    # but it's easy to add support for other detectors
    
    # set detector to ccd or apd
    if(config.ccd_detector is not None): 
        detector=config.ccd_detector
        # Need to remove the following later on - we multiply by 
        # 1e3 only because of a bug in the current script 
        input_photon=input_photon
        # print('using ccd detector')
    elif(config.apd_detector is not None): 
        detector=config.apd_detector
        # print('using apd detector')
    else: 
        raise ValueError("Failed to set detector. Must be either apd or ccd.")

    #Get detector size
    detector_row=detector.geometry.row
    detector_col=detector.geometry.col

    #Pad zeros to array to match detector size
    det_shape = (detector_row, detector_col)
    position_y=0 
    position_x=0
    padded_pyxel_photon_array = fit_into_array(array=input_photon,
    output_shape=det_shape, align='center')

    photon_array = padded_pyxel_photon_array
    try:
        detector.photon.array += photon_array
    except ValueError as ex:
        raise ValueError("Shapes of arrays do not match") from ex

    #run pyxel in exposure mode
    pyxel_result = pyxel.exposure_mode(
        exposure=exposure, detector=detector, pipeline=pipeline)
    #print(pyxel_result)

    #Get the image from Pyxel
    pyxel_readout=config.exposure.readout.times[0]
    output_image_padded=np.array(pyxel_result.image.sel(readout_time=pyxel_readout))

    #Remove padding
    output_image = fit_into_array(
        array=output_image_padded,
        output_shape=(row,col), align='center'
    )

    t2=time.time()
    # print(t2-t1, 'seconds taken for Pyxel computing')

    return np.array(output_image)