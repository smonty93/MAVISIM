import numpy as np
import math
from scipy import interpolate
from mavisim import generate_image
import astropy.io.fits as fits
from tqdm import tqdm


class AstromCalibSimGeneric():
    """Generic class for astrometric calibration simulation.
    See: AstromCalibSimAna for the analytical simulation
         AstromCalibSimE2E for the end-to-end simulation
    """

    def __init__(self, static_distort, *, pixel_size_as=0.00736,
                 pixel_size_mm=10e-3, dist_amp=1.0, mask_scale=1e3 / 582,
                 hole_position_std=0.0, dx=0.2, dy=0.2,
                 dx_meas=None, dy_meas=None, n_poly=6, pin_pitch=0.5,
                 num_pin_x=40):
        self._true_cam_samp = pixel_size_as                    # arcsec/pixel at output plane
        self._plate_scale = self._true_cam_samp / pixel_size_mm  # arcsec/mm at output plane
        self._static_distort = static_distort                  # static distortion data
        self._dist_amp = dist_amp                              # distortion amplification factor
        self._mask_scale = mask_scale                          # arcsec/mm at input plane
        self._pin_pitch = pin_pitch                            # mm at input plane
        self._num_pins = num_pin_x**2                          # 30x30 pinholes
        self._num_pin_x = num_pin_x                            # number of pins across the x-dimension (sqrt(num_pins_total))
        self._hole_position_std = hole_position_std
        self._n_poly = n_poly
        self._dx = dx
        self._dy = dy
        if dx_meas is not None:
            self._dx_meas = dx_meas
        else:
            self._dx_meas = dx
        if dy_meas is not None:
            self._dy_meas = dy_meas
        else:
            self._dy_meas = dy
        self._n_tot_poly = ((self._n_poly + 1) * (self._n_poly + 2)) // 2 - 1

        # Process input distortions into interpolated evaluatable function
        if static_distort is None:
            self._input_distortions_func = lambda x, y: np.c_[x*y*0, x*y*0]
        else:
            self._input_distortions_func = self._input_distortions()
        self._recovered_distortions_func = None
        self._p0_meas = None
        self._ppx_meas = None
        self._ppy_meas = None

    def input_dist(self, x, y):
        """Evaluate interpolated input distortions at arbitrary coordinates.

        `x` and `y` (in arcsec) can be anywhere in the science field, but must
        be array-like and the same size.

        Args:
            x : array-like float : field x-coordinates (arcsec)
            y : array-like float : field y-coordinates (arcsec)

        Returns
            out_x : array-like float : x-component of distortion at each coord
            out_y : array-like float : y-component of distortion at each coord
        """
        x = np.array(x).copy()
        y = np.array(y).copy()
        xx = x.flatten()
        yy = y.flatten()
        out_xx = xx * 0
        out_yy = xx * 0
        for i in range(xx.shape[0]):
            out_xx[i], out_yy[i] = self._input_distortions_func(xx[i], yy[i])[0]
        out_xx = out_xx.reshape(x.shape)
        out_yy = out_yy.reshape(x.shape)
        return out_xx, out_yy

    def _input_distortions_degmm(self):
        field_x = self._static_distort["Field_x(deg)"]
        field_y = self._static_distort["Field_y(deg)"]

        dist_x = self._static_distort["Predicted_x(mm)"] - self._static_distort["Real_x(mm)"]
        dist_y = self._static_distort["Predicted_y(mm)"] - self._static_distort["Real_y(mm)"]

        grid_vals_x = np.unique(field_x)
        grid_vals_y = np.unique(field_y)
        # Create grids to save the distortion (difference) in x and y at each point in the grid
        # This is necessary for the interpolation of the distortion (grid-wise interpolation)
        dist_x_grid = np.zeros([len(grid_vals_x), len(grid_vals_y)])
        dist_y_grid = np.zeros([len(grid_vals_x), len(grid_vals_y)])

        for i in range(len(field_x)):
            row = np.where(grid_vals_x == field_x[i])
            col = np.where(grid_vals_y == field_y[i])
            dist_x_grid[row, col] = dist_x[i]
            dist_y_grid[row, col] = dist_y[i]

        dist_x_func_degmm = interpolate.RectBivariateSpline(grid_vals_x, grid_vals_y, dist_x_grid)
        dist_y_func_degmm = interpolate.RectBivariateSpline(grid_vals_x, grid_vals_y, dist_y_grid)

        return lambda x, y: np.c_[dist_x_func_degmm(x, y), dist_y_func_degmm(x, y)] * self._dist_amp

    def _input_distortions(self):
        """Distortion function generator wrapper for arcsec to arcsec.
        Call this function to generate
        """
        out_func = self._input_distortions_degmm()
        return lambda x, y: out_func(x / 3600, y / 3600) * self._plate_scale

    def recovered_dist(self, x, y):
        """Evaluate recovered/estimated input distortions at arbitrary coordinates.
        This is the estimated distortion via the differential calibration method
        based on the input static distortion.

        `x` and `y` (in arcsec) can be anywhere in the science field, but must
        be array-like and the same size.

        Args:
            x : array-like float : field x-coordinates (arcsec)
            y : array-like float : field y-coordinates (arcsec)

        Returns
            out_x : array-like float : x-component of distortion at each coord
            out_y : array-like float : y-component of distortion at each coord
        """
        if self._recovered_distortions_func is None:
            raise ValueError("Generic astrometric simulator does not implement distortion recovery method")
        x = np.array(x).copy()
        y = np.array(y).copy()
        xx = x.flatten()
        yy = y.flatten()
        out_xx = xx * 0
        out_yy = xx * 0
        for i in range(xx.shape[0]):
            out_xx[i], out_yy[i] = self._recovered_distortions_func(xx[i], yy[i])[0]
        out_xx = out_xx.reshape(x.shape)
        out_yy = out_yy.reshape(x.shape)
        return out_xx, out_yy

    def residual_dist(self, x, y):
        """Evaluate residual distortions at arbitrary coordinates.

        `x` and `y` (in arcsec) can be anywhere in the science field, but must
        be array-like and the same size.

        Args:
            x : array-like float : field x-coordinates (arcsec)
            y : array-like float : field y-coordinates (arcsec)

        Returns
            out_x : array-like float : x-component of distortion at each coord
            out_y : array-like float : y-component of distortion at each coord
        """
        if self._recovered_distortions_func is None:
            raise ValueError("Generic astrometric simulator does not implement distortion recovery method")

        x = np.array(x).copy()
        y = np.array(y).copy()
        xx = x.flatten()
        yy = y.flatten()
        out_xx = xx * 0
        out_yy = xx * 0
        for i in range(xx.shape[0]):
            out_xx[i], out_yy[i] = self._input_distortions_func(xx[i], yy[i])[0] - \
                self._recovered_distortions_func(xx[i], yy[i])[0]
        out_xx = out_xx.reshape(x.shape)
        out_yy = out_yy.reshape(x.shape)
        return out_xx, out_yy

    @staticmethod
    def _hbvpoly(p, a, n_poly):
        """ Evaluate the homogenous bi-variate polynomial defined by
        coefficients in a at position p.

        Arguments:
            p: np.ndarray : position to evaluate polynomial at, (M,2)
            a: np.ndarray : coefficients defining polynomial, (((N+1)(N+2))//2-1,)
            N: int: maximum homogenous polynomial order to go to.

        Returns:
            out: np.ndarray : evaluated polynomial, scalar or (M,1)
        """
        if len(p.shape) != 2:
            raise ValueError("p must be 2D, i.e., p.shape=(M,2)")
        out = np.zeros_like(p[:, 0])
        counter = 0
        for n in range(1, n_poly + 1):  # skip tip-tilt
            for j in range(n + 1):
                out[:] += a[counter] * p[:, 0]**j * p[:, 1]**(n - j) / math.factorial(n)
                counter += 1
        return out

    @staticmethod
    def _hbvpoly_grad(p, n_poly):
        """ Evaluate the gradient of the homogenous bi-variate polynomial
        defined by coefficients in a at position p.

        Arguments:
            p: np.ndarray : position to evaluate polynomial gradient at, (2,) or (M,2)
            n_poly: int: maximum homogenous polynomial order to go to.

        Returns:
            out: np.ndarray : evaluated polynomial gradient,

        """
        dx = np.zeros([p.shape[0], ((n_poly + 1) * (n_poly + 2)) // 2 - 1])
        dy = np.zeros([p.shape[0], ((n_poly + 1) * (n_poly + 2)) // 2 - 1])
        counter = -1
        for n in range(1, n_poly + 1):  # skip tip-tilt
            for j in range(n + 1):
                counter += 1
                if j == 0:
                    continue
                dx[:, counter] += j * p[:, 0]**(j - 1) * p[:, 1]**(n - j) / math.factorial(n)

        counter = -1
        for n in range(1, n_poly + 1):  # skip tip-tilt
            for j in range(n + 1):
                counter += 1
                if j == n:
                    continue
                dy[:, counter] += (n - j) * p[:, 0]**j * p[:, 1]**(n - j - 1) / math.factorial(n)
        return dx, dy

    @staticmethod
    def _make_pinhole_grid(xshift=0., yshift=0., sigma=0., grid="square", incl_dist=True,
                           pins_per_side=30, mask_scale=1e3 / 582, pin_pitch=0.571,
                           plate_scale=7.36e-3 / 10e-3, dist_func_degmm=None, seed=1234):
        """
        Generate arrays of x-y pinhole positions in pixels and arcsec.
        Optionally pass a global shift in x/y, in mm, to shift pinhole grid
        relative to distortion field. Can also provide an uncertainty in the hole
        positions (also in mm), which is treated as Gaussian. Distortions are
        included by default, but can be turned off to get "nominal" pinhole grid.


        Args:
            xshift (float, optional): Shift amount in x-axis (arcsec). Defaults to 0..
            yshift (float, optional): Shift amount in y-axis (arcsec). Defaults to 0..
            sigma (float, optional): Standard deviation on pinhole positions. Defaults to 0..
            incl_dist (bool, optional): Flag to include distortions in returned coordinates. Defaults to True.
            pins_per_side (int, optional): Number of pins per side of the grid. Defaults to 30.
            mask_scale (float, optional): arcsec/mm at the mask. Defaults to (1e3/582).
            pin_pitch (float, optional): spacing between the pinholes (in mm). Defaults to 0.582.
            plate_scale (float, optional): arcsec/mm at the sensor. Defaults to 7.36e-3/10e-3.
            dist_func_degmm (function, optional): Function to obtain distortion.
                Takes argument in degrees in field, returns distortion in mm at sensor. Defaults to None.

        Returns:
            numpy.ndarray: 2d array of x-y pinhole positions in arcsec.
        """
        if grid == "square":
            # individual hole positions in mm (relative to centre)
            x_mm, y_mm = np.meshgrid(
                (np.arange(pins_per_side) - (pins_per_side - 1) / 2.) * pin_pitch,
                (np.arange(pins_per_side) - (pins_per_side - 1) / 2.) * pin_pitch
            )
        elif grid == "hex":
            n_init = pins_per_side
            spacing = pin_pitch
            x_mm, y_mm = np.meshgrid(np.linspace(-n_init / 2, n_init / 2, n_init + 1),
                                     np.linspace(-n_init / 2, n_init / 2, n_init + 1))
            y_mm[:, ::2] += 1 / 2
            x_mm *= (np.sqrt(3) / 2)
            y_mm -= 0.5
            x_mm *= spacing
            y_mm *= spacing
        else:
            raise ValueError("grid type: " + grid + " is unsupported.")
        mask = (np.abs(x_mm) <= 15.0/mask_scale) * (np.abs(y_mm) <= 15.0 / mask_scale)
        x_mm = x_mm[mask]
        y_mm = y_mm[mask]
        x_mm = x_mm.flatten()
        y_mm = y_mm.flatten()

        pin_loc_mm_in = np.c_[x_mm, y_mm]

        # apply any shift terms
        pin_loc_mm_in[:, 0] += xshift
        pin_loc_mm_in[:, 1] += yshift

        valid = (pin_loc_mm_in**2).sum(axis=1)**0.5 <= (14.5 / mask_scale)

        # include uncertainty in hole position?
        np.random.seed(seed)
        pin_loc_mm_in += sigma * np.random.randn(*pin_loc_mm_in.shape)

        # convert from mm to deg. Use full grid to include individual point uncertainties
        pin_loc_as = pin_loc_mm_in * mask_scale

        if incl_dist:
            if dist_func_degmm is None:
                raise ValueError("Must provide distortion functions to include distortion")
            # get distortion terms at each location
            for ii, loc_deg in enumerate(pin_loc_as):
                pin_loc_as[ii, :] += plate_scale * dist_func_degmm(loc_deg[0] / 3600, loc_deg[1] / 3600).flatten()

        return pin_loc_as, valid

    def _fit_poly(self):
        """_summary_

        Returns:
            _type_: _description_
        """
        n_pos = self._valid.sum()

        d_mat = np.zeros([4 * n_pos, 2 * self._n_tot_poly])
        grad_tmp = self._hbvpoly_grad(self._p0_nom[self._valid], self._n_poly)

        d_mat[0::4, :self._n_tot_poly] = grad_tmp[0]
        d_mat[1::4, self._n_tot_poly:] = grad_tmp[0]
        d_mat[2::4, :self._n_tot_poly] = grad_tmp[1]
        d_mat[3::4, self._n_tot_poly:] = grad_tmp[1]

        d_inv = np.linalg.solve(d_mat.T @ d_mat, d_mat.T)
        self._d_inv = d_inv

        # compenent-wise gradients:
        if self._p0_meas is None:
            raise RuntimeError("measurements must be made before fitting polynomial")

        dx_arcsec = self._mask_scale * self._dx_meas
        dy_arcsec = self._mask_scale * self._dy_meas
        dpdx = (self._ppx_meas - self._p0_meas) - np.r_[dx_arcsec, 0]
        dpdx /= dx_arcsec
        dpdy = (self._ppy_meas - self._p0_meas) - np.r_[0, dy_arcsec]
        dpdy /= dy_arcsec

        # estimated gradients:
        z_hat = np.c_[dpdx[self._valid], dpdy[self._valid]].flatten()
        # estimated polynomial coefficients:
        return d_inv @ z_hat

    def _set_coords_square(self):
        kwargs_grid = {
            "pins_per_side": self._num_pin_x,
            "mask_scale": self._mask_scale,
            "pin_pitch": self._pin_pitch,
            "plate_scale": self._plate_scale,
            "dist_func_degmm": self._input_distortions_degmm(),
            "grid": "square"
        }
        self._p0_nom, v = self._make_pinhole_grid(incl_dist=False, **kwargs_grid)
        self._ppx_nom, v_ = self._make_pinhole_grid(xshift=self._dx, incl_dist=False, **kwargs_grid)
        v = v * v_
        self._ppy_nom, v_ = self._make_pinhole_grid(yshift=self._dy, incl_dist=False, **kwargs_grid)
        v = v * v_
        self._valid = v
        self._p0_true, _ = self._make_pinhole_grid(sigma=self._hole_position_std, **kwargs_grid)
        self._ppx_true, _ = self._make_pinhole_grid(xshift=self._dx, sigma=self._hole_position_std, **kwargs_grid)
        self._ppy_true, _ = self._make_pinhole_grid(yshift=self._dy, sigma=self._hole_position_std, **kwargs_grid)

    def _set_coords_hex(self):
        kwargs_grid = {
            "pins_per_side": self._num_pin_x,
            "mask_scale": self._mask_scale,
            "pin_pitch": self._pin_pitch,
            "plate_scale": self._plate_scale,
            "dist_func_degmm": self._input_distortions_degmm(),
            "grid": "hex"
        }
        self._p0_nom, v = self._make_pinhole_grid(incl_dist=False, **kwargs_grid)
        self._ppx_nom, v_ = self._make_pinhole_grid(xshift=self._dx, incl_dist=False, **kwargs_grid)
        v = v * v_
        self._ppy_nom, v_ = self._make_pinhole_grid(yshift=self._dy, incl_dist=False, **kwargs_grid)
        v = v * v_
        self._valid = v
        self._p0_true, _ = self._make_pinhole_grid(sigma=self._hole_position_std, **kwargs_grid)
        self._ppx_true, _ = self._make_pinhole_grid(xshift=self._dx, sigma=self._hole_position_std, **kwargs_grid)
        self._ppy_true, _ = self._make_pinhole_grid(yshift=self._dy, sigma=self._hole_position_std, **kwargs_grid)

###############################################################################
###############################################################################
###############################################################################


class AstromCalibSimAna(AstromCalibSimGeneric):
    """
    """

    def __init__(self, *args, centroid_noise_std=0.0, **kwargs):
        super().__init__(*args, **kwargs)
        self._centroid_noise_std = centroid_noise_std
        # Perform calibration process, and process recovered distortions into evaluatable function
        self._recovered_distortions_func = self._recovered_distortions_ana()

    def _do_measurements(self):
        rand_shape = self._p0_nom.shape
        self._p0_meas = self._p0_true + np.random.randn(*rand_shape) * self._centroid_noise_std
        self._ppx_meas = self._ppx_true + np.random.randn(*rand_shape) * self._centroid_noise_std
        self._ppy_meas = self._ppy_true + np.random.randn(*rand_shape) * self._centroid_noise_std

    def _recovered_distortions_ana(self):
        self._set_coords_square()
        self._do_measurements()
        u_hat = self._fit_poly()
        return lambda x, y: np.c_[self._hbvpoly(np.c_[x, y], u_hat[:self._n_tot_poly], self._n_poly),
                                  self._hbvpoly(np.c_[x, y], u_hat[self._n_tot_poly:], self._n_poly)]

###############################################################################
###############################################################################
###############################################################################


class AstromCalibSimE2E(AstromCalibSimGeneric):
    """
    """

    def __init__(self, *args, pin_size=10e-3, pinhole_os=4, pixel_os=2, wavelength=550e-9, pinhole_support_width=128,
                 noise_fun=None, centroid_win_rad=0.2, centroid_threshold=0.0, **kwargs):
        super().__init__(*args, **kwargs)
        self._pin_size = pin_size           # pinhole physical diameter in mm
        self._pinhole_os = pinhole_os       # oversampling factor
        self._pixel_os = pixel_os           # pixel oversampling factor
        self._wavelength = wavelength       # wavelength in metres
        self._fine_cam_samp = self._true_cam_samp / pixel_os
        self._noise_fun = noise_fun
        self._centroid_win_rad = centroid_win_rad
        self._centroid_threshold = centroid_threshold
        self._pinhole_support_width = pinhole_support_width

        # Perform calibration process, and process recovered distortions into evaluatable function
        self._recovered_distortions_func = self._recovered_distortions_e2e()

    def _recovered_distortions_e2e(self):
        self._set_coords_hex()
        self._do_measurements()
        u_hat = self._fit_poly()
        return lambda x, y: np.c_[self._hbvpoly(np.c_[x, y], u_hat[:self._n_tot_poly], self._n_poly),
                                  self._hbvpoly(np.c_[x, y], u_hat[self._n_tot_poly:], self._n_poly)]

    def _centroids(self, pos, im, origin):
        """_summary_

        Args:
            pos (_type_): _description_
            im (_type_): _description_
            win_rad (float, optional): _description_. Defaults to 0.2.

        Returns:
            _type_: _description_
        """
        cog_meas = []
        win_rad = self._centroid_win_rad
        # iterate over pinholes:
        for pin_x, pin_y in tqdm(pos, desc="centroids", leave=False):
            # integer-valued indices for desired pinhole-window
            win_idx = np.mgrid[int((pin_y - win_rad - origin[1]) / self._true_cam_samp):
                               int((pin_y + win_rad - origin[1]) / self._true_cam_samp) + 1,
                               int((pin_x - win_rad - origin[0]) / self._true_cam_samp):
                               int((pin_x + win_rad - origin[0]) / self._true_cam_samp) + 1]
            win_idx = np.clip(win_idx, 0, im.shape[0] - 1)
            # corresponding pixel coordinates within that window
            win_as = [((w+0.5) * self._true_cam_samp + origin[0]).flatten() for w in win_idx]
            # pixel intensities of window, flattened for centroiding
            window = im[win_idx[0], win_idx[1]].flatten()
            # centroid calculation
            cog_meas.append([(window @ p) / window.sum() for p in win_as][::-1])
        return np.r_[cog_meas]

    def _do_measurements(self):
        # Do this by making an image and measuring the positions:
        pin_image = self._pinhole_image()
        pin_image = pin_image / (self._pixel_os**2)
        pin_image_filename = "_pin_image.fits"
        fits.writeto(pin_image_filename, np.array([[]]), overwrite=True)
        fits.append(pin_image_filename, pin_image,
                    fits.Header({"YPOS": 0.0, "XPOS": 0.0, "LAMBDA": self._wavelength}))

        class SourceHack:
            def __init__(self, coords):
                self.flux = 1.0 * np.ones(len(coords))
                self.gauss_pos = coords.copy()
                self.cov_mat = None

        source_p0 = SourceHack(self._p0_true)
        source_ppx = SourceHack(self._ppx_true)
        source_ppy = SourceHack(self._ppy_true)

        image_gen_p0 = generate_image.ImageGenerator(6400 * self._pixel_os,
                                                     source_p0, pin_image_filename,
                                                     self._fine_cam_samp, which_psf=0, norm_psf=False)
        image_gen_ppx = generate_image.ImageGenerator(6400 * self._pixel_os,
                                                      source_ppx, pin_image_filename,
                                                      self._fine_cam_samp, which_psf=0, norm_psf=False)
        image_gen_ppy = generate_image.ImageGenerator(6400 * self._pixel_os,
                                                      source_ppy, pin_image_filename,
                                                      self._fine_cam_samp, which_psf=0, norm_psf=False)
        image_gen_p0.main()
        image_gen_ppx.main()
        image_gen_ppy.main()

        self._im0 = image_gen_p0.get_rebinned_cropped(self._pixel_os, self._true_cam_samp * 4000)
        self._impx = image_gen_ppx.get_rebinned_cropped(self._pixel_os, self._true_cam_samp * 4000)
        self._impy = image_gen_ppy.get_rebinned_cropped(self._pixel_os, self._true_cam_samp * 4000)

        if self._noise_fun is not None:
            self._im0 = self._noise_fun(self._im0)
            self._impx = self._noise_fun(self._impx)
            self._impy = self._noise_fun(self._impy)

        self._im0 -= self._centroid_threshold
        self._impx -= self._centroid_threshold
        self._impy -= self._centroid_threshold
        self._im0[self._im0 < 0] = 0
        self._impx[self._impx < 0] = 0
        self._impy[self._impy < 0] = 0

        self._p0_meas = np.zeros(self._p0_nom.shape)
        self._ppx_meas = np.zeros(self._ppx_nom.shape)
        self._ppy_meas = np.zeros(self._ppy_nom.shape)

        # TODO: this origin thing is ugly and maybe buggy, can be cleaned up with a convention change
        origin = np.array([-0.5, -0.5]) * (self._true_cam_samp * 4000)
        self._p0_meas[self._valid] = self._centroids(pos=self._p0_nom[self._valid], im=self._im0, origin=origin)
        self._ppx_meas[self._valid] = self._centroids(pos=self._ppx_nom[self._valid], im=self._impx, origin=origin)
        self._ppy_meas[self._valid] = self._centroids(pos=self._ppy_nom[self._valid], im=self._impy, origin=origin)

    @staticmethod
    def _pinhole(size, x, y, radius):
        """generate the sampled pinhole function

        Args:
            size (int): size of the output array (in pixels)
            x (float): x position of the pinhole (in pixels)
            y (float): y position of the pinhole (in pixels)
            radius (float): radius of the pinhole (in pixels)

        Returns:
            numpy.ndarray: 2d array of the pinhole function
        """
        # generate coordinates:
        pos = np.mgrid[:size, :size]
        # signed distance to pinhole edge:
        dr = np.sqrt((pos[1] - x)**2 + (pos[0] - y)**2) - radius
        # generate the pinhole function:
        weight = 0.5 - np.clip(dr, -0.5, 0.5)
        return weight

    def _pinhole_image(self):
        """generate the pinhole image

        Args:
            output_size (int): size of the output array (in pixels)

        Returns:
            numpy.ndarray: 2d array of the pinhole image
        """
        import poppy
        from astropy.convolution import convolve_fft
        # Pinhole Representation, fourier transform the pinhole directly
        pin_ap = (self._pin_size * self._mask_scale) * self._pinhole_os / self._fine_cam_samp
        output_size = self._pinhole_support_width
        # 2d point grid
        newsize = output_size * self._pinhole_os
        centre = newsize / 2

        # make a pinhole
        pin = self._pinhole(newsize, centre, centre, pin_ap)

        # generate the PSF with poppy
        osys = poppy.OpticalSystem()
        osys.add_pupil(poppy.CircularAperture(radius=4.))    # pupil radius in meters
        osys.add_detector(pixelscale=self._fine_cam_samp, fov_pixels=output_size, oversample=self._pinhole_os)
        psf = osys.calc_psf(self._wavelength)
        tel_psf = psf[0].data / psf[0].data.sum()

        # smoosh them together
        pin_image = convolve_fft(pin, tel_psf, allow_huge=True)

        # bin back down
        rb = np.copy(pin_image).reshape(
            output_size, self._pinhole_os, output_size, self._pinhole_os
        ).sum(axis=(-1, 1)) / self._pinhole_os**2

        return rb
