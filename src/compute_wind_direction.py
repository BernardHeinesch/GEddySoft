import numpy as np

def compute_wind_direction(Usonic, Vsonic, SonicType, SonicNorthAzimuth):
    """
    see https://www.eol.ucar.edu/content/wind-direction-quick-reference
    (but there seems to be an error, with arguments of arctg2 eronously signed,
     and the routine is here vectorized for performance improvement)
    Meteorological wind coordinate system: Umet, Vmet
    A positive Umet component represents wind blowing to the East. +Vmet is
    wind to the North. This is right handed with respect to an upward +Wmet.
    ::
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
    Usonic : float or np array of float
        sonic u wind speed component
    Vsonic : float or np array of float
        sonic v wind speed component
    SonicType : string
        type of sonic coordinate system (see respective sonic manuals)
        (either 'R2Axis','R3Spar','R3Axis','HS50Spar','HS50Axis' or 'CSAT3')
    SonicNorthAzimuth : float
        direction to which the sonic is pointing to (=0 if pointing to North)

    Returns
    -------
    Dirmet : float or np array of float
        direction the wind is coming from (usual meteorological convention)

    comments
    --------
    Written by B. Heinesch.
    University of Liege, Gembloux Agro-Bio Tech.
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
