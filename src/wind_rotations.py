import numpy as np
import sys


def wind_rotations(ini, sonicdata, mean_uvw=False, R_tilt_PFM=False, mean_wind_dir=False):
    """
    Coordinate rotations for sonic anemometer data.

    parameters
    ----------
    ini: dict, initialisation information
    sonicdata: np array of floats, data from the sonic and the IRGA
    mean_uvw: np array of floats, mean wind components
    R_tilt_PFM: dict of 3x3 np arrays, needed only for the angle-dependent planar fit
    mean_wind_dir: np array of floats, mean wind direction, needed only for the angle-dependent planar fit
    
    returns
    -------
    wind_rotated: np array of floats, rotated mean wind components
    rot_phi: float, pitch angle (rad)
    R_tilt: 3x3 np array of floats, rotation matrix

    comments
    --------
    Written by B. Heinesch, March, 2025
    University of Liege, Gembloux Agro-Bio Tech.

    """

    if (ini['param']['TILT_CORRECTION_MODE'] == 0):
        # no rotation
        wind_rotated = [sonicdata[:, 1],
                        sonicdata[:, 2],
                        sonicdata[:, 3]]  # untilted u, v and w

    elif ini['param']['TILT_CORRECTION_MODE'] == 1:
        # double rotation (yaw and pitch)
        rot_theta = np.arctan2(mean_uvw[1], mean_uvw[0])  #yaw
        u1 = mean_uvw[0]*np.cos(rot_theta) + mean_uvw[1]*np.sin(rot_theta)
        w1 = mean_uvw[2]
        rot_phi = np.arctan2(w1, u1)  # pitch

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

        rot_theta = np.arctan2(mean_uvw_PFM[1], mean_uvw_PFM[0])

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
