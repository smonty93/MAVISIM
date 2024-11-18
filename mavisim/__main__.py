#!/usr/bin/env python
import argparse
from mavisim.generate_image import ImageGenerator
from pydantic import BaseModel, ConfigDict
import numpy as np
from astropy.io import fits  # type: ignore
from scipy.interpolate import LinearNDInterpolator  # type: ignore
import json
from typing import Callable

imshape = (4000, 4000)


class Source(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)
    exp_time: float  # exposure time in seconds to simulate.
    star: np.ndarray  # unique ID of each star.
    flux: np.ndarray  # flux of each star.
    gauss_pos: np.ndarray  # X/Y-position of each star (in arcsec).
    gauss_cov: np.ndarray  # covariance of Gaussian kernel
    static_dist: None = None  # legacy
    cov_mat: None = None  # legacy

    @property
    def pos(self) -> np.ndarray:
        return self.gauss_pos

    @pos.setter
    def pos(self, value):
        self.gauss_pos = value

    def trim(self):
        valid = np.ones(self.pos.shape[0], dtype=bool)
        valid &= self.pos[:, 0] >= -20
        valid &= self.pos[:, 0] <= 20
        valid &= self.pos[:, 1] >= -20
        valid &= self.pos[:, 1] <= 20
        self.star = self.star[valid, ...]
        self.flux = self.flux[valid, ...]
        self.gauss_pos = self.gauss_pos[valid, ...]
        self.gauss_cov = self.gauss_cov[valid, ...]


def read_source(filename: str) -> Source:
    # read source file:
    with open(filename) as f:
        lines = f.readlines()
    try:
        float(lines[0].split()[0])
    except ValueError:
        # assume first row is headers
        lines = lines[1:]

    # create clean source data
    table = [line.split() for line in lines if len(line) > 0]
    xx = np.r_[[float(line[3]) for line in table]]
    yy = np.r_[[float(line[5]) for line in table]]
    flux = np.r_[[float(line[7]) for line in table]]
    nstars = xx.shape[0]

    # build mavisim-compatible Source object
    source = Source(
        exp_time=1.0,
        star=np.arange(nstars),
        flux=flux,
        gauss_pos=np.concatenate([xx[:, None], yy[:, None]], axis=1),
        gauss_cov=np.tile((np.eye(2)*1e-9)[None, :, :], reps=[nstars, 1, 1])
    )
    return source


def read_distortions(filename: str) -> Callable:
    """read distortion input file and return bivariate vector valued
    distortion function: (x_as,y_as) -> (dx_as,dy_as)"""
    platescale = 0.736  # "/mm

    with open(filename) as f:
        lines_raw = f.readlines()
    # remove comments
    lines_data = [line for line in lines_raw if line[0] != "#"]
    try:
        float(lines_data[0].split()[0])
    except ValueError:
        # assume first row is headers
        lines_data = lines_data[1:]
    lines = [[float(ell) for ell in line.split()] for line in lines_data]
    pos_x_mm = np.r_[[line[4] for line in lines]]
    pos_y_mm = np.r_[[line[5] for line in lines]]
    dist_x_mm = np.r_[[line[6] - line[4] for line in lines]]
    dist_y_mm = np.r_[[line[7] - line[5] for line in lines]]
    pos_x_as = pos_x_mm * platescale
    pos_y_as = pos_y_mm * platescale
    dist_x_as = dist_x_mm * platescale
    dist_y_as = dist_y_mm * platescale
    interp = LinearNDInterpolator(
        np.array([pos_x_as, pos_y_as]).T,
        np.array([dist_x_as, dist_y_as]).T,
        fill_value=0.0,
    )
    return interp


def main():
    parser = argparse.ArgumentParser(
        "mavisim",
        description="mavis image simulator",
    )
    parser.add_argument(
        "field", type=str, help="file containing source field to simulate",
    )
    parser.add_argument(
        "psfs", type=str,
        help="specify PSFs fits file on disk to use in simulation",
    )
    parser.add_argument(
        "-o", "--output", type=str, default="./out.fits",
        help="specify name for output fits file",
    )
    parser.add_argument(
        "-d", "--dists", type=str,
        help="specify distortion file on disk to pre-distort sources with",
    )
    parser.add_argument(
        "--static", "-s", type=int,
        help="use only 1 psf (not field varying), at the given psf fits hdu",
    )
    parser.add_argument(
        "--header", type=str, default=r'{}',
        help="dict to be parsed as header for saved fits file",
    )
    parser.add_argument(
        "--posfile", type=str, help="file to store distorted positions",
    )

    args = parser.parse_args()

    source = read_source(args.field)

    try:
        header_entries = json.loads(args.header)
    except json.decoder.JSONDecodeError as e:
        print(argparse.ArgumentError(parser._actions[-1], str(e)))
        exit(1)

    # read and apply distortions (if specified)
    if args.dists is not None:
        interp = read_distortions(args.dists)
        dist = interp(source.pos)
        source.pos += dist

    source.trim()

    if args.posfile:
        with open(args.posfile, "w") as f:
            for sr, fl in zip(source.pos, source.flux):
                f.write(f"{sr[0]},{sr[1]},{fl}\n")

    imgen = ImageGenerator(
        12_000, source, psfs_file=args.psfs, which_psf=args.static,
    )

    imgen.main()
    im = imgen.get_rebinned_cropped(2, 30)

    header = fits.Header()
    header.update(header_entries)
    hdu = fits.hdu.PrimaryHDU(data=im, header=header)
    hdu.writeto(args.output, overwrite=True)
    print(f"saved to: {args.output}")


if __name__ == "__main__":
    main()
