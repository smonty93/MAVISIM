import numpy as np
from scipy import interpolate

class AstromCalibSim():
    """
    Astrometric Calibration Simulator for MAVIS.

    MAVIS will feature an astrometric calibration mask characterise the static
    distortions present in the MAVIS system. This simulator allows one to quantify
    the performance of that mask for a given distortion field, using the nominal
    astrometric calibration identification algorithms intended for MAVIS.

    The constructor takes in an astropy-parsed distortion field with the appropriate 
    headers (see below), and performs a simulated astrometric calibration process. 
    The public methods allow the evaluation of input, recovered, and residual 
    distortion fields.

    Static distortion file should have at least the following columns:
    ```    
      Field_x(deg)    Field_y(deg)    Predicted_x(mm)  Predicted_y(mm)    Real_x(mm)      Real_y(mm)
    -4.16666667E-03 -4.16666667E-03   2.02613981E+01   2.02626354E+01   2.03494040E+01  2.03513749E+01
    -4.16666667E-03 -4.07407407E-03   2.02613981E+01   1.98123546E+01   2.03496805E+01  1.98994423E+01
    -4.16666667E-03 -3.98148148E-03   2.02613981E+01   1.93620738E+01   2.03499497E+01  1.94474891E+01
    ...
     4.07407407E-03  4.07407407E-03  -1.98111448E+01  -1.98123546E+01  -1.98997456E+01 -1.98980094E+01
    ```
    and it should be parsed by astropy first, like:
    ```python
    from astropy import ascii
    static_distort = ascii(distort_file)
    ```

    Args:
        static_distort : astropy table : table containing the distortions across the field
        
        centroid_noise_std : float : standard deviation of Gaussian noise applied to centroids (arcsec).
        
        hole_position_std  : float : standard deviation of manufacturing error on hole positions (mm).
        
        dx : float : shift applied in x direction for calibration process (mm).
        
        dy : float : shift applied in y direction for calibration process (mm).
        
        n_poly : int : maximum order of homogenous bivariate polynomial used to fit distortions.

        pin_pitch : float : distance between pins at input plane (mm).

        num_pin_x : int : number of pins across the x-dimension (sqrt(num_pins_total)).

        mask_scale : float : arcsec/mm at input plane.

        pixel_size_as : float : pixel size in arcsec of imager camera.

        pixel_size_mm : float : physical pixel size in mm of imager camera.

    Attributes:
        mask_scale : float : arcsec/mm at input plane.

        static_distort : astropy.Table : static distortion input table
    """
    def __init__(self, static_distort, centroid_noise_std=10e-6,
                hole_position_std=1e-2, dx=0.2, dy=0.2, n_poly=6, pin_pitch=0.5,
                num_pin_x=30, mask_scale=1e3/582, pixel_size_as=0.00736, 
                pixel_size_mm=10e-3):
        self.mask_scale = mask_scale                           # arcsec/mm at input plane
        self._init_cam_samp = pixel_size_as                    # arcsec/pixel at output plane
        self._pin_pitch = pin_pitch                            # mm at input plane
        self._num_pins = num_pin_x**2                          # 30x30 pinholes
        self._plate_scale = self._init_cam_samp/pixel_size_mm  # arcsec/mm at output plane
        self.static_distort = static_distort                   # static distortion data
        
        # Process input distortions into interpolated evaluatable function
        self._input_distortions_func = self._input_distortions()

        # Perform calibration process, and process recovered distortions into evaluatable function
        self._recovered_distortions_func = self._recovered_distortions(
                    centroid_noise_std=centroid_noise_std, 
                    hole_position_std=hole_position_std, 
                    dx=dx, dy=dy, n_poly=n_poly)

    def input_dist(self,x,y):
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
        out_xx = xx*0
        out_yy = xx*0
        for i in range(xx.shape[0]):
            out_xx[i],out_yy[i] = self._input_distortions_func(xx[i],yy[i])[0]
        out_xx = out_xx.reshape(x.shape)
        out_yy = out_yy.reshape(x.shape)
        return out_xx,out_yy
    
    def recovered_dist(self,x,y):
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
        x = np.array(x).copy()
        y = np.array(y).copy()
        xx = x.flatten()
        yy = y.flatten()
        out_xx = xx*0
        out_yy = xx*0
        for i in range(xx.shape[0]):
            out_xx[i],out_yy[i] = self._recovered_distortions_func(xx[i],yy[i])[0]
        out_xx = out_xx.reshape(x.shape)
        out_yy = out_yy.reshape(x.shape)
        return out_xx,out_yy

    def residual_dist(self,x,y):
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
        x = np.array(x).copy()
        y = np.array(y).copy()
        xx = x.flatten()
        yy = y.flatten()
        out_xx = xx*0
        out_yy = xx*0
        for i in range(xx.shape[0]):
            out_xx[i],out_yy[i] = self._input_distortions_func(xx[i],yy[i])[0] - \
                                self._recovered_distortions_func(xx[i],yy[i])[0]
        out_xx = out_xx.reshape(x.shape)
        out_yy = out_yy.reshape(x.shape)
        return out_xx,out_yy

    def _input_distortions_degmm(self):
        field_x = self.static_distort["Field_x(deg)"]
        field_y = self.static_distort["Field_y(deg)"]

        dist_x = self.static_distort["Predicted_x(mm)"] - self.static_distort["Real_x(mm)"]
        dist_y = self.static_distort["Predicted_y(mm)"] - self.static_distort["Real_y(mm)"]

        # Create an array of the field positions and the distortion at each pt (the difference)
        dist_all = np.empty([dist_x.shape[0], 4])
        dist_all[:, 0] = field_x
        dist_all[:, 1] = field_y
        dist_all[:, 2] = dist_x
        dist_all[:, 3] = dist_y

        grid_vals = np.unique(field_x)

        # Create grids to save the distortion (difference) in x and y at each point in the grid
        # This is necessary for the interpolation of the distortion (grid-wise interpolation)
        dist_x_grid = np.zeros([len(grid_vals),len(grid_vals)])
        dist_y_grid = np.zeros([len(grid_vals),len(grid_vals)])

        num = 0

        for x in grid_vals:
            sub_array = dist_all[np.where(dist_all[:, 0] == x), :][0]

            for row in np.arange(0, sub_array.shape[0]):

                dist_x_grid[num, row] = sub_array[row, 2]
                dist_y_grid[num, row] = sub_array[row, 3]

            num+=1

        dist_x_func_degmm = interpolate.RectBivariateSpline(grid_vals, grid_vals, dist_x_grid)
        dist_y_func_degmm = interpolate.RectBivariateSpline(grid_vals, grid_vals, dist_y_grid)
        
        return lambda x,y: np.c_[dist_x_func_degmm(x,y), dist_y_func_degmm(x,y)]

    def _input_distortions(self):
        """Distortion function generator wrapper for arcsec to arcsec.
        Call this function to generate 
        """
        out_func = self._input_distortions_degmm()
        return lambda x,y: out_func(x/3600,y/3600)*self._plate_scale

    def _recovered_distortions(self,centroid_noise_std=10e-6,
                        hole_position_std=1e-2,dx=0.2,dy=0.2,n_poly=6):

        # Static distortion map (good enough for now) & functions
        dist_func_degmm = self._input_distortions_degmm()

        def make_pinhole_grid(xshift=0., yshift=0., sigma=None):
            """
            Generate arrays of x-y pinhole positions in pixels and arcsec.
            Optionally pass a global shift in x/y, in mm, to shift pinhole grid
            relative to distortion field. Can also provide an uncertainty in the hole
            positions (also in mm), which is treated as Gaussian.
            
            """
            #holes per side assuming square grid
            pins_per_side = np.sqrt(self._num_pins)

            #individual hole positions in mm (relative to centre)
            pin_loc_x_mm = (np.arange(pins_per_side)-(pins_per_side-1)/2.)*self._pin_pitch
            pin_loc_y_mm = (np.arange(pins_per_side)-(pins_per_side-1)/2.)*self._pin_pitch

            # make full grid of pinhole positions
            pin_grid_x_mm, pin_grid_y_mm = np.meshgrid(pin_loc_x_mm, pin_loc_y_mm)

            #apply any shift terms
            pin_grid_x_mm += xshift
            pin_grid_y_mm += yshift
            
            #include uncertainty in hole position?
            if sigma is not None:
                rng = np.random.default_rng(1234)
                perturb = rng.multivariate_normal([0,0], [[1,0],[0,1]], pin_grid_x_mm.shape)
                
                pin_grid_x_mm += perturb[:,:,0]*sigma
                pin_grid_y_mm += perturb[:,:,1]*sigma
            
            # convert from mm to deg. Use full grid to include individual point uncertainties
            pin_grid_x_deg = pin_grid_x_mm * self.mask_scale / 3600.
            pin_grid_y_deg = pin_grid_y_mm * self.mask_scale / 3600.

            # get distortion terms at each location
            x_pos_dist_mm, y_pos_dist_mm = np.zeros(pin_grid_x_deg.size), np.zeros(pin_grid_x_deg.size)
            for ii, x, y in zip(np.arange(pin_grid_x_deg.size),pin_grid_x_deg.flatten(), pin_grid_y_deg.flatten()):
                x_pos_dist_mm[ii] = dist_func_degmm(x,y)[0,0]
                y_pos_dist_mm[ii] = dist_func_degmm(x,y)[0,1]
            x_pos_dist_mm = np.atleast_2d(x_pos_dist_mm).reshape(pin_grid_x_deg.shape)
            y_pos_dist_mm = np.atleast_2d(y_pos_dist_mm).reshape(pin_grid_x_deg.shape)

            # get "effective" position including optical distortion
            pin_eff_x_arcsec = pin_grid_x_deg*3600 + x_pos_dist_mm*self._plate_scale
            pin_eff_y_arcsec = pin_grid_y_deg*3600 + y_pos_dist_mm*self._plate_scale
            
            return np.c_[pin_eff_x_arcsec.flatten(),
                         pin_eff_y_arcsec.flatten()]

        def get_nominal_pos(xshift=0.0,yshift=0.0):
            """
            xshift,yshift in mm
            """
            #individual hole positions in arcsec (relative to centre)
            pins_per_side = np.sqrt(self._num_pins)
            #individual hole positions in mm (relative to centre)
            pin_loc_x_mm = (np.arange(pins_per_side)-(pins_per_side-1)/2.)*self._pin_pitch
            pin_loc_y_mm = (np.arange(pins_per_side)-(pins_per_side-1)/2.)*self._pin_pitch
            # make full grid of pinhole positions
            pin_grid_x_mm, pin_grid_y_mm = np.meshgrid(pin_loc_x_mm, pin_loc_y_mm)
            #apply any shift terms
            pin_grid_x_mm += xshift
            pin_grid_y_mm += yshift

            pin_eff_x_arcsec = pin_grid_x_mm * self.mask_scale
            pin_eff_y_arcsec = pin_grid_y_mm * self.mask_scale
            # get "effective" position including optical distortion
            return np.c_[pin_eff_x_arcsec.flatten(),
                         pin_eff_y_arcsec.flatten()]

        def hbvpoly(p,a,n_poly):
            """ Evaluate the homogenous bi-variate polynomial defined by 
            coefficients in a at position p.
            
            Arguments:
                p: np.ndarray : position to evaluate polynomial at, (M,2)
                a: np.ndarray : coefficients defining polynomial, (((N+1)(N+2))//2-1,)
                N: int: maximum homogenous polynomial order to go to.
            
            Returns:
                out: np.ndarray : evaluated polynomial, scalar or (M,1)
            """
            if len(p.shape)!=2:
                raise ValueError("p must be 2D, i.e., p.shape=(M,2)")
            out = np.zeros_like(p[:,0])
            counter = 0
            for n in range(1,n_poly+1): # skip tip-tilt
                for j in range(n+1):
                    out[:] += a[counter]*p[:,0]**j*p[:,1]**(n-j)/np.math.factorial(n)
                    counter += 1
            return out

        def hbvpoly_grad(p,n_poly):
            """ Evaluate the gradient of the homogenous bi-variate polynomial 
            defined by coefficients in a at position p.
            
            Arguments:
                p: np.ndarray : position to evaluate polynomial gradient at, (2,) or (M,2)
                n_poly: int: maximum homogenous polynomial order to go to.

            Returns:
                out: np.ndarray : evaluated polynomial gradient,
            
            """
            dx = np.zeros([p.shape[0],((n_poly+1)*(n_poly+2))//2-1])
            dy = np.zeros([p.shape[0],((n_poly+1)*(n_poly+2))//2-1])
            counter = -1
            for n in range(1,n_poly+1): # skip tip-tilt
                for j in range(n+1):
                    counter += 1
                    if j==0:
                        continue
                    dx[:,counter] += j*p[:,0]**(j-1)*p[:,1]**(n-j)/np.math.factorial(n)

            counter = -1
            for n in range(1,n_poly+1): # skip tip-tilt
                for j in range(n+1):
                    counter += 1
                    if j==n:
                        continue
                    dy[:,counter] += (n-j)*p[:,0]**j*p[:,1]**(n-j-1)/np.math.factorial(n)
            return dx,dy

        nominal_pos = get_nominal_pos()
        p0  = make_pinhole_grid(xshift=0,yshift=0,sigma=hole_position_std)
        ppx = make_pinhole_grid(xshift=dx,yshift=0,sigma=hole_position_std)
        ppy = make_pinhole_grid(xshift=0,yshift=dy,sigma=hole_position_std)
        
        dx_arcsec = self.mask_scale*dx
        dy_arcsec = self.mask_scale*dy
        
        n0 = np.random.randn(*p0.shape)*centroid_noise_std
        npx = np.random.randn(*ppx.shape)*centroid_noise_std
        npy = np.random.randn(*ppy.shape)*centroid_noise_std
        
        p0 += n0
        ppx += npx
        ppy += npy
        
        n_tot = ((n_poly+1)*(n_poly+2))//2-1
        n_pos = nominal_pos.shape[0]

        d_mat = np.zeros([4*n_pos,2*n_tot])
        grad_tmp = hbvpoly_grad(nominal_pos,n_poly)

        d_mat[0::4,:n_tot]   = grad_tmp[0]
        d_mat[1::4,n_tot:]   = grad_tmp[0]
        d_mat[2::4,:n_tot]   = grad_tmp[1]
        d_mat[3::4,n_tot:]   = grad_tmp[1]

        d_inv = np.linalg.solve(d_mat.T@d_mat,d_mat.T)

        # compenent-wise gradients:
        dpdx = (ppx-p0) - np.r_[dx_arcsec,0]
        dpdx /= dx_arcsec
        dpdy = (ppy-p0) - np.r_[0,dy_arcsec]
        dpdy /= dy_arcsec
        
        # estimated gradients:
        z_hat = np.c_[dpdx,dpdy].flatten()

        # estimated polynomial coefficients:
        u_hat = d_inv @ z_hat
        #return {    "d_true" : np.c_[make_pinhole_grid(xshift=0,yshift=0,sigma=None)] - nominal_pos,
        #        "d_estimate" : np.c_[hbvpoly(nominal_pos,u_hat[:N_tot],N),hbvpoly(nominal_pos,u_hat[N_tot:],N)]}
        
        return lambda x,y: np.c_[hbvpoly(np.c_[x,y],u_hat[:n_tot],n_poly),
                                 hbvpoly(np.c_[x,y],u_hat[n_tot:],n_poly)]