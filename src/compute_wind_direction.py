"""Module for computing wind direction from sonic anemometer measurements.

This module provides functions to convert wind measurements from sonic anemometer
coordinate systems to meteorological wind direction conventions. It handles
different sonic anemometer types (Gill, Campbell, etc.) and accounts for their
specific mounting orientations.

See Also
--------
UCAR EOL Wind Direction Quick Reference:
    https://www.eol.ucar.edu/content/wind-direction-quick-reference
    (but there seems to be an error, with arguments of arctg2 eronously signed,
     and the routine is here vectorized for performance improvement)
Author
------
B. Heinesch
University of Liege, Gembloux Agro-Bio Tech

Created
-------
2020
"""

import numpy as np

def compute_wind_direction(Usonic, Vsonic, SonicType, SonicNorthAzimuth):
    """Convert sonic anemometer wind components to meteorological wind direction.

    This function converts wind velocity components measured in sonic anemometer
    coordinates to meteorological wind direction (direction FROM which the wind
    is blowing). It handles different sonic types and mounting orientations.
    The function is vectorized for improved performance.

    Coordinate Systems
    ==================
    Meteorological (target):
        - Umet: Positive = wind blowing to East
        - Vmet: Positive = wind blowing to North
        - Forms right-handed system with upward Wmet
        - Wind direction: 0° = North, 90° = East (direction wind is FROM)

    Sonic (source):
        - Coordinate system varies by manufacturer
        - Vsonic axis angle (Vaz) depends on type:
            * CSAT3: Vaz = array direction - 90°
            * Gill R2: Vaz = N arrow + 240° (U flipped)
            * Gill R3/HS50 AXIS: Vaz = N arrow + 240°
            * Gill R3/HS50 SPAR: Vaz = N arrow + 270°

    Coordinate Diagram::

                v>0
                ^
                |
                |
                +----> U>0

    Converting between Sonic and Meteorological Coordinates
    Determine the angle with respect to true north, (0=N,90=E) of the +Vsonic
    direction.
    Call this angle Vaz. Looking from above, the sonic coordinate system is
    therefore rotated Vaz degrees clockwise from meteorological coordinates.

    For ATI and Campbell CSAT3 sonics, Vaz is the direction relative to
    true north, straight into the array from the un-obstructed direction,
    minus 90 degrees.
    For Gill R2 sonics, if the sign of V is flipped, then Vaz is the angle
    of the N arrow + 60. If the sign of U is flipped, then Vaz is the N arrow
    direction + 240.
    For Gill R3s in the AXIS configuration, Vaz is the N arrow direction 
    + 240 degrees.

    From this Meteorological Coordinates, the direction the  wind is blowing 
    FROM is finaly computed 
    Parameters
    ----------
    Usonic : array_like
        U wind component in sonic coordinates [m/s]
    Vsonic : array_like
        V wind component in sonic coordinates [m/s]
    SonicType : str
        Sonic anemometer model/configuration:
        - 'R2': Gill R2
        - 'R3Spar': Gill R3 SPAR configuration
        - 'R3Axis': Gill R3 AXIS configuration
        - 'HS50Spar': Gill HS50 SPAR configuration
        - 'HS50Axis': Gill HS50 AXIS configuration
        - 'CSAT3': Campbell CSAT3
    SonicNorthAzimuth : float
        Azimuth angle of sonic anemometer reference direction
        relative to true North [degrees] (=0 if pointing to North)

    Returns
    -------
    Dirmet : ndarray
        Wind direction in meteorological coordinates [degrees]
        (0-360°, direction wind is coming FROM)

    Notes
    -----
    The conversion process:
    1. Convert U,V components to direction in sonic coordinates
    2. Apply sonic-specific rotation (Vaz)
    3. Correct for sonic mounting orientation (SonicNorthAzimuth)
    4. Convert from 'wind blowing TO' to 'wind blowing FROM'
    """

    # Ensure Usonic and Vsonic are arrays for vectorized operations
    Usonic = np.asarray(Usonic)
    Vsonic = np.asarray(Vsonic)

    if SonicType == 'R2':
        Usonic = -Usonic
        Vaz = 240
    if SonicType == 'R3Spar' or SonicType == 'HS50Spar':
        Vaz = 270
    if SonicType == 'R3Axis' or SonicType == 'HS50Axis':
        Vaz = 240
    if SonicType == 'CSAT3':
        Vaz = 180-90

    # compute wind direction in the sonic coordinate system
    Dirsonic = np.arctan2(Usonic, Vsonic) * 180./np.pi
    Dirsonic = np.where(Dirsonic < 0, Dirsonic + 360, Dirsonic)  # adjust negative angles

    # switch to the met coordinate system
    Dirmet = Dirsonic + Vaz
    Dirmet = np.where(Dirmet > 360, Dirmet - 360, Dirmet)  # adjust angles > 360

    # correct for SonicNorthAzimuth
    Dirmet = (Dirmet + SonicNorthAzimuth) % 360

    # switch from the "to" the "from" wind direction convention
    Dirmet = np.where(Dirmet < 180, Dirmet + 180, Dirmet - 180)

    return Dirmet
