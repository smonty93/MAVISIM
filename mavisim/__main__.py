#!/usr/bin/env python
import argparse
from mavisim.generate_image import ImageGenerator
from pydantic import BaseModel, ConfigDict
import numpy as np
from astropy.io import fits
from scipy.interpolate import LinearNDInterpolator
import json


class Source(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)
    exp_time: float = None  # exposure time in seconds to simulate.
    star: np.ndarray = None  # unique ID of each star.
    flux: np.ndarray = None  # flux of each star.
    gauss_pos: np.ndarray = None  # X/Y-position of each star (in arcsec).
    gauss_cov: np.ndarray = None  # covariance of Gaussian kernel
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
    with open(args.field) as f:
        lines = f.readlines()
    try:
        float(lines[0].split()[0])
    except ValueError:
        # assume first row is headers
        lines = lines[1:]

    # create clean source data
    lines = [line.split() for line in lines if len(line) > 0]
    xx = np.r_[[float(line[3]) for line in lines]]
    yy = np.r_[[float(line[5]) for line in lines]]
    flux = np.r_[[float(line[7]) for line in lines]]
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


def read_distortions(filename: str) -> callable:
    with open(filename) as f:
        lines = f.readlines()
    lines = [line for line in lines if line[0] != "#"]  # remove comments
    try:
        float(lines[0].split()[0])
    except ValueError:
        # assume first row is headers
        lines = lines[1:]
    lines = [line.split() for line in lines]
    pos_x_deg = np.r_[[float(line[0]) for line in lines]]
    pos_y_deg = np.r_[[float(line[1]) for line in lines]]
    dist_x_mm = np.r_[[
        float(line[6]) - float(line[4])
        for line in lines
    ]]
    dist_y_mm = np.r_[[
        float(line[7]) - float(line[5])
        for line in lines
    ]]
    pos_x_as = pos_x_deg * 3600
    pos_y_as = pos_y_deg * 3600
    dist_x_as = dist_x_mm * 0.736  # 0.736 as/mm
    dist_y_as = dist_y_mm * 0.736  # 0.736 as/mm
    # spacing = pos_y_as[1] - pos_y_as[0]
    # x_vec = np.arange(pos_x_as.min(), pos_x_as.max()+spacing, spacing)
    # y_vec = np.arange(pos_y_as.min(), pos_y_as.max()+spacing, spacing)
    # dist_x_array = np.reshape(dist_x_as, [x_vec.shape[0], y_vec.shape[0]])
    # dist_y_array = np.reshape(dist_y_as, [x_vec.shape[0], y_vec.shape[0]])
    # dist_xx = RectBivariateSpline(x_vec, y_vec, dist_x_array).ev(xx, yy)
    # dist_yy = RectBivariateSpline(x_vec, y_vec, dist_y_array).ev(xx, yy)
    interp = LinearNDInterpolator(
        np.array([pos_x_as, pos_y_as]).T,
        np.array([dist_x_as, dist_y_as]).T,
        fill_value=0.0,
    )
    return interp


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
    help="specify distortion file on disk to pre-distort sources with"
)
parser.add_argument(
    "--static", "-s", type=int,
    help="use only 1 psf (not field varying), at the given psf fits hdu index"
)
parser.add_argument(
    "--header", type=str, default=r'{}',
    help="dict to be parsed as header for saved fits file"
)

args = parser.parse_args()

source = read_source(args.field)

try:
    header_entries = json.loads(args.header)
except json.decoder.JSONDecodeError as e:
    print(argparse.ArgumentError(parser._actions[-1], e))
    exit(1)

# read and apply distortions (if specified)
if args.dists is not None:
    interp = read_distortions(args.dists)
    dist = interp(source.pos)
    source.pos += dist

source.trim()

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
