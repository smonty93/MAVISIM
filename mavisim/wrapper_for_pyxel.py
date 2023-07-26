# Module to calculate noise using the Pyxel tool

import pyxel
from pyxel.util import fit_into_array
import numpy as np
import numpy.typing as npt
from pyxel.detectors import CCD, APD


def calc_noise_pyxel(
    input_photon: npt.NDArray,
    yaml_file: str,
) -> np.ndarray | None:
    """ Takes photon array and yaml parameter file name and processes it to add noise

    Parameters
    ----------
    input_photon: ndarray (pass photon array)
    yaml_file: yaml file containing detector properties

    Returns
    -------
    output_image: ndarray

    """
    # Add pyxel noise calculation stuff here
    # Make sure that the photon only has positive elements
    input_photon[input_photon < 0.0] = 0.0
    row, col = np.shape(input_photon)

    # Get Pyxel configs from the yaml file

    config = pyxel.load(yaml_file=yaml_file)

    if config.exposure is None:
        return None

    # Set pyxel stuff from the config file
    exposure = config.exposure  # class Single
    pipeline = config.pipeline  # class DetectionPipeline

    # Currently this only supports ccd and apd detector,
    # but it's easy to add support for other detectors

    detector: APD | CCD
    # set detector to ccd or apd
    if (config.ccd_detector is not None):
        detector = config.ccd_detector
    elif (config.apd_detector is not None):
        detector = config.apd_detector
    else:
        raise ValueError("Failed to set detector. Must be either apd or ccd.")

    # Get detector size
    detector_row = detector.geometry.row
    detector_col = detector.geometry.col

    # Pad zeros to array to match detector size
    det_shape = (detector_row, detector_col)
    padded_pyxel_photon_array = fit_into_array(array=input_photon,
                                               output_shape=det_shape, align='center')

    photon_array = padded_pyxel_photon_array
    try:
        detector.photon.array += photon_array
    except ValueError as ex:
        raise ValueError("Shapes of arrays do not match") from ex

    # run pyxel in exposure mode
    if exposure:
        pyxel_result = pyxel.exposure_mode(
            exposure=exposure, detector=detector, pipeline=pipeline)
    else:
        raise NotImplementedError

    # Get the image from Pyxel
    if config.exposure:
        pyxel_readout = config.exposure.readout.times[0]
        output_image_padded = np.array(pyxel_result.image.sel(readout_time=pyxel_readout))

    # Remove padding
    output_image = fit_into_array(
        array=output_image_padded,
        output_shape=(row, col), align='center'
    )

    return np.array(output_image)
