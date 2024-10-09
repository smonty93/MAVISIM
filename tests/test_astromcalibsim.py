#!/usr/bin/env python
# coding: utf-8

from mavisim.astromsim import AstromCalibSimAna, AstromCalibSimE2E, AstromCalibSimGeneric
from astropy.io import ascii
import numpy as np
import pytest


def eval_astromsim(astrom_sim):
    # Values to evaluate distortion functions with. Can be any arangement of
    # coordinates within the science field.
    nsamp = 31
    xx, yy = np.meshgrid(np.linspace(-15, 15, nsamp), np.linspace(-15, 15, nsamp))
    xx = xx.flatten()
    yy = yy.flatten()
    rr = (xx**2+yy**2)**0.5
    xx = xx[rr <= 15.0]
    yy = yy[rr <= 15.0]

    # interpolate input distortions at xx/yy coordinates
    input_dist_xx, input_dist_yy = astrom_sim.input_dist(xx, yy)

    # interpolate recovered distortions at xx/yy coordinates
    recovered_dist_xx, recovered_dist_yy = astrom_sim.recovered_dist(xx, yy)

    # interpolate residual distortions at xx/yy coordinates
    residual_dist_xx, residual_dist_yy = astrom_sim.residual_dist(xx, yy)
    # remove tt and plate scale since not in the requirements
    X = np.array([xx, yy, xx*0 + 1]).T
    coeffs = np.linalg.solve(X.T @ X, X.T) @ np.c_[residual_dist_xx, residual_dist_yy]
    filt = np.eye(xx.shape[0]) - (X @ np.linalg.solve(X.T @ X, X.T))

    input_dist_xx = filt @ input_dist_xx
    input_dist_yy = filt @ input_dist_yy
    recovered_dist_xx = filt @ recovered_dist_xx
    recovered_dist_yy = filt @ recovered_dist_yy
    residual_dist_xx = filt @ residual_dist_xx
    residual_dist_yy = filt @ residual_dist_yy

    # Evaluate astrometric error between objects
    n_objects = 4000
    objects = np.random.random([n_objects, 2])*30-15
    filt_objects = ((objects**2).sum(axis=1)**0.5) <= 15.0
    objects = objects[filt_objects, :]
    n_objects = objects.shape[0]
    dists = []
    errs = []
    true_pos = objects.copy()
    distort = np.array([astrom_sim.residual_dist(ob[0], ob[1]) for ob in objects]) - objects @ coeffs[:2, :]
    distance = ((true_pos[:, None, :]-true_pos[None, :, :])**2).sum(axis=-1)**0.5
    errors = ((distort[:, None, :]-distort[None, :, :])**2).sum(axis=-1)**0.5
    indices = np.tril_indices(distance.shape[0], 1)
    dists = distance[indices]
    errs = errors[indices]
    rel_err = errs[dists < 1]
    rel_err_n = (dists < 1).sum()
    print(f"  rel astrometric error: {rel_err.std()*1e6:7.3f} uas rms")
    print(f"           err < 150uas: {(rel_err < 150e-6).sum() * 100 / rel_err_n:6.2f}%")
    print(f"           err <  50uas: {(rel_err < 50e-6).sum() * 100 / rel_err_n:6.2f}%")
    abs_err = errs[:]
    abs_err_n = errs.shape[0]
    print(f"  abs astrometric error: {abs_err.std()*1e6:7.3f} uas rms")
    print(f"          err < 2000uas: {(abs_err < 2000e-6).sum() * 100 / abs_err_n:6.2f}%")
    print(f"          err <  400uas: {(abs_err < 400e-6).sum() * 100 / abs_err_n:6.2f}%")
    return errs


ana_testdata = [
    (
        0.01,  # shift
        5,  # npoly
        0.0,  # noise
        20e-6  # rms performance threshold
    ),
    (
        0.3,
        5,
        10e-6,
        50e-6
    ),
    (
        0.3,
        5,
        20e-6,
        200e-6
    ),
]


@pytest.mark.parametrize("dx,n_poly,noise,thresh", ana_testdata)
def test_astromsim_analytical(dx, n_poly, noise, thresh):
    # input true distortion field:
    static_distort = ascii.read("tests/test_StaticDistort")

    # create the astrometry simulator object:
    astrom_sim = AstromCalibSimAna(
        static_distort,            # astropy parsed input file with appropriate headers
        centroid_noise_std=noise,  # centroiding noise in arcsec.
        dx=dx, dy=dx,              # Shift applied to mask for differential method.
        n_poly=n_poly)             # max order of polynomial to fit with.
    errs = eval_astromsim(astrom_sim)
    assert errs.std() < thresh


def noise_func(image, flux):
    image = np.random.poisson(image*flux)*1.0
    return image


e2e_testdata = [
    (
        0.3,  # shift
        5,  # npoly
        lambda x: x,  # noise func
        20e-6  # rms performance threshold
    ),
    (
        0.3,
        5,
        lambda x: noise_func(x, 6e5),
        50e-6
    ),
    (
        0.3,
        5,
        lambda x: noise_func(x, 1e5),
        200e-6
    ),
]


@pytest.mark.parametrize("dx,n_poly,noise_func,thresh", e2e_testdata)
def test_astromsim_e2e(dx, n_poly, noise_func, thresh):
    # input true distortion field:
    static_distort = ascii.read("tests/test_StaticDistort")
    # create the astrometry simulator object:
    astrom_sim = AstromCalibSimE2E(
        static_distort, num_pin_x=40, dx=dx, dy=dx, n_poly=n_poly,
        noise_fun=noise_func
        )
    errs = eval_astromsim(astrom_sim)
    assert errs.std() < thresh


def build_reconstructor(*, dx, dy,  # offset applied (in pixels)
                        n_poly,     # number of polynomials to use in reconstruction
                        p0):        # home positions of pinholes (in pixels)
    # create the astrometry simulator object:
    astrom_sim = AstromCalibSimGeneric(
        None,
        mask_scale=1.0,
        n_poly=n_poly,
        dx=dx, dy=dy)
    valid = np.ones(p0.shape[0], dtype=bool)
    astrom_sim._p0_nom = p0
    astrom_sim._valid = valid
    try:
        astrom_sim._fit_poly()
    except RuntimeError:
        pass
    return astrom_sim


def test_reconstructor():
    """This test also serves as an example for building the astrometric
    calibration reconstructor.
    """
    # replace this with array of actual pinhole positions:
    pinhole_positions = np.array([
        [0.0, 0.0],
        [30.0, 0.0],
        [0.0, 30.0],
        [30.0, 30.0],
        [0.0, 60.0],
        [60.0, 60.0],
        ])

    parameters = {
        "n_poly": 5,  # number of polynomials in each dimensions
        "dx": 1.0,    # pixels of shift in x
        "dy": 1.0,    # pixels of shift in y
        "p0": pinhole_positions,  # in pixels, [n_pinholes,2] numpy array
    }

    simulator = build_reconstructor(**parameters)
    recon = simulator._d_inv

    # now that you have the reconstructor, you need the positions after shifts
    # have been applied.
    dx_vector = np.r_[parameters["dx"], 0]
    dy_vector = np.r_[0, parameters["dy"]]

    # As an example, if there are no distortions present:
    positions_after_dx = (pinhole_positions + dx_vector)
    positions_after_dy = (pinhole_positions + dy_vector)

    # The sampled gradient is the measured distortion divided by the size of the
    # displacement induced:
    dpdx = (positions_after_dx) - (pinhole_positions) - dx_vector
    dpdy = (positions_after_dy) - (pinhole_positions) - dy_vector
    dpdx /= parameters["dx"]
    dpdy /= parameters["dy"]

    # estimated gradient samples in a vector:
    z_hat = np.c_[dpdx, dpdy].flatten()

    # estimated polynomial coefficients:
    distortion_modes = recon @ z_hat

    # evaluating the distortions at arbitrary co-ordinates:
    def evaluate(x, y):
        return np.c_[
            simulator._hbvpoly(np.c_[x, y], distortion_modes[:simulator._n_tot_poly], simulator._n_poly),
            simulator._hbvpoly(np.c_[x, y], distortion_modes[simulator._n_tot_poly:], simulator._n_poly)
        ]

    # here are the distortions (in pixels) evaluated at the pinhole positions
    dists = evaluate(pinhole_positions[:, 0], pinhole_positions[:, 1])
    assert (np.sum(dists**2) < 1e-5)


if __name__ == "__main__":
    test_astromsim_analytical(*ana_testdata[0])
    test_astromsim_e2e(*e2e_testdata[0])
    test_reconstructor()
