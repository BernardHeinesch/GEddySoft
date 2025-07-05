import numpy as np
import sys


def wind_rotations(ini, sonicdata, mean_uvw=False, R_tilt_PFM=False, mean_wind_dir=False):
    """
    Apply coordinate rotations to sonic anemometer wind measurements.

    This function implements three different tilt correction methods for sonic
    anemometer data to align the coordinate system with the mean streamlines:

    1. No rotation (TILT_CORRECTION_MODE = 0)
    2. Double rotation (TILT_CORRECTION_MODE = 1)
       - First rotation (yaw) around z-axis to nullify mean lateral wind
       - Second rotation (pitch) around y-axis to nullify mean vertical wind
    3. Angle-dependent planar fit (TILT_CORRECTION_MODE = 2)
       - Applies pre-calculated rotation matrices for different wind sectors
       - Followed by a yaw rotation

    Parameters
    ----------
    ini : dict
        Initialization dictionary containing processing parameters.
        Must include 'param' key with 'TILT_CORRECTION_MODE' setting:
        - 0: No rotation
        - 1: Double rotation
        - 2: Angle-dependent planar fit
    sonicdata : numpy.ndarray
        Raw sonic anemometer data array  (n_samples, x)
        Minimum columns are [timestamp, u, v, w] where:
        - u: streamwise wind component [m s⁻¹]
        - v: crosswind component [m s⁻¹]
        - w: vertical wind component [m s⁻¹]
    mean_uvw : numpy.ndarray, optional
        Mean wind components [u, v, w]. Required for double rotation.
        Default is False.
    R_tilt_PFM : dict, optional
        Dictionary of pre-calculated 3x3 rotation matrices for planar fit.
        Keys are wind direction sector bounds in degrees.
        Required for angle-dependent planar fit.
        Default is False.
    mean_wind_dir : float, optional
        Mean wind direction in degrees for sector selection.
        Required for angle-dependent planar fit.
        Default is False.

    Returns
    -------
    wind_rotated : numpy.ndarray
        Rotated wind components array with shape (n_samples, 3)
        Columns are [u, v, w] in the rotated coordinate system
    rot_phi : float
        Pitch angle in radians (only meaningful for double rotation)
        Returns NaN for planar fit method
    R_tilt : numpy.ndarray
        Final 3x3 rotation matrix that was applied to the data

    Notes
    -----
    For the double rotation method:
    1. Yaw rotation (θ): tan(θ) = v/u
    2. Pitch rotation (φ): tan(φ) = w/u

    The rotation matrix R combines both rotations:
    R = R_pitch * R_yaw

    For the planar fit method, pre-calculated rotation matrices for each wind
    sector are applied first, followed by a yaw rotation to align with the
    mean wind direction.

    References
    ----------
    .. [1] Wilczak, J.M., Oncley, S.P. and Stage, S.A. (2001). Sonic
           anemometer tilt correction algorithms. Boundary-Layer
           Meteorology, 99(1), 127-150.
    .. [2] Lee, X., Massman, W. and Law, B. (2004). Handbook of
           Micrometeorology: A Guide for Surface Flux Measurement and
           Analysis. Chapter 3. Kluwer Academic Publishers.

    See Also
    --------
    compute_wind_direction : Calculates wind direction from u and v components
    test_steady_state_FW96 : Tests stationarity of rotated wind components

    Examples
    --------
    >>> # Example with double rotation
    >>> import numpy as np
    >>> # Generate sample 10 Hz data for 30 minutes
    >>> n_samples = 18000
    >>> sonicdata = np.zeros((n_samples, 4))
    >>> sonicdata[:, 1:] = np.random.normal(0, 0.5, (n_samples, 3))
    >>> sonicdata[:, 1] += 2  # Add mean streamwise wind
    >>> mean_uvw = np.mean(sonicdata[:, 1:], axis=0)
    >>> ini = {'param': {'TILT_CORRECTION_MODE': 1}}
    >>> wind_rot, phi, R = wind_rotations(ini, sonicdata, mean_uvw)

    Author
    ------
    Written by Bernard Heinesch (March 2025)
    University of Liege, Gembloux Agro-Bio Tech
    """

    if (ini['param']['TILT_CORRECTION_MODE'] == 0):
        # no rotation
        wind_rotated = [sonicdata[:, 1],
                        sonicdata[:, 2],
                        sonicdata[:, 3]]  # untilted u, v and w

    elif ini['param']['TILT_CORRECTION_MODE'] == 1:
        # double rotation (yaw and pitch)
        rot_theta = np.arctan(mean_uvw[1]/mean_uvw[0])  #yaw
        u1 = mean_uvw[0]*np.cos(rot_theta) + mean_uvw[1]*np.sin(rot_theta)
        w1 = mean_uvw[2]
        rot_phi = np.arctan(w1/u1)  # pitch

        R_tilt = np.array([
            [np.cos(rot_theta)*np.cos(rot_phi) ,  np.sin(rot_theta)*np.cos(rot_phi), np.sin(rot_phi)],
            [-np.sin(rot_theta)                ,     np.cos(rot_theta)             ,       0        ],
            [-np.cos(rot_theta)*np.sin(rot_phi), -np.sin(rot_theta)*np.sin(rot_phi), np.cos(rot_phi)]
                            ])
        wind_rotated = np.transpose(np.dot(R_tilt, [
            np.transpose(sonicdata[:, 1]),
            np.transpose(sonicdata[:, 2]),
            np.transpose(sonicdata[:, 3])
            ]))  # untilted u, v and w

    elif (ini['param']['TILT_CORRECTION_MODE'] == 2):
        # angle-dependent planar fit

        # Get the planar fit matrix (pitch and roll) for the specific wind sector
        keys = list(R_tilt_PFM.keys())
        for i in range(len(keys)-1):
            lower_range = keys[i]
            upper_range = keys[i+1]
            if lower_range <= mean_wind_dir < upper_range:
                R_tilt = R_tilt_PFM[lower_range]
                break
        # Handle the case where mean_wind_dir is above the last key
        else:
            # If no match is found in the loop, assign the last key's value
            R_tilt = R_tilt_PFM[keys[-1]]

        # apply pitch and roll rotations
        wind_detilted_PFM = np.transpose(np.dot(R_tilt, [
            np.transpose(sonicdata[:, 1]),
            np.transpose(sonicdata[:, 2]),
            np.transpose(sonicdata[:, 3])
            ]))  # untilted u, v and w

        # apply yaw rotation
        mean_uvw_PFM = np.nanmean(wind_detilted_PFM[:, 0:3], 0)

        rot_theta = np.arctan(mean_uvw_PFM[1]/mean_uvw_PFM[0])

        R_tilt = np.array([
            [np.cos(rot_theta) , np.sin(rot_theta) , 0],
            [-np.sin(rot_theta), np.cos(rot_theta) , 0],
            [        0         ,         0         , 1]
                            ])

        wind_rotated = np.transpose(np.dot(R_tilt, [
            np.transpose(wind_detilted_PFM[:, 0]),
            np.transpose(wind_detilted_PFM[:, 1]),
            np.transpose(wind_detilted_PFM[:, 2])
            ]))

        # unused
        rot_phi = np.nan

    else:
        # invalid choice
        sys.exit('TILT_CORRECTION_MODE incorrectly set')

    return wind_rotated, rot_phi, R_tilt
