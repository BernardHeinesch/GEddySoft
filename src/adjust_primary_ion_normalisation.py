import numpy as np


def adjust_primary_ion_normalisation(tracerdata):
    """
    Adjustment of tracer concentration normalisation by primary ion counts.

    parameters
    ----------
    tracerdata: dict
        Dictionary holding the tracer input data for the current averaging interval.
        Required keys:
            - conc: np.ndarray, shape (n_time, n_tracers)
            - Xr0: np.ndarray, shape (n_tracers,)
            - FPH1: float (root attribute)
            - FPH2: float (root attribute)
            - conc_primary_21_022: np.ndarray, shape (n_time,)
            - conc_primary_38_033: np.ndarray, shape (n_time,)

    returns
    -------
    tracerdata: dict
        Same dict with tracerdata['conc'] normalised in-place.

    comments
    --------
    Implements the alternative primary-ion normalisation (Loubet et al., 2021, Biogeosciences).
    This function does not compute a new concentration from scratch; it adjusts the existing
    normalisation of tracerdata['conc'] by applying a time-dependent scaling factor based on
    primary ion counts.

    The goal is to convert a tracer signal that is normalised with instantaneous (high-frequency)
    primary ion counts into a signal normalised with averaged primary ion counts over the
    averaging interval.

    Practically, for each tracer j:
        - de-normalise using instantaneous primary ion counts (high-frequency cps)
        - re-normalise using the averaged primary ion counts over the interval

    The normalisation applied to each tracer j is:
        c_j(t) <- c_j(t) * [FPH1*c21(t) + Xr0_j*FPH2*c38(t)] / [FPH1*mean(c21) + Xr0_j*FPH2*mean(c38)]
    """

    conc_21 = tracerdata['conc_primary_21_022']
    conc_38 = tracerdata['conc_primary_38_033']
    conc_21_mean = np.nanmean(conc_21)
    conc_38_mean = np.nanmean(conc_38)

    norm_inst = (tracerdata['FPH1'] * conc_21[:, None] +
                 tracerdata['FPH2'] * conc_38[:, None] * tracerdata['Xr0'][None, :])
    norm_mean = (tracerdata['FPH1'] * conc_21_mean +
                 tracerdata['FPH2'] * conc_38_mean * tracerdata['Xr0'])

    tracerdata['conc'] = tracerdata['conc'] * (norm_inst / norm_mean[None, :])

    return tracerdata
