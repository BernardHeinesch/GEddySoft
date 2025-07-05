import numpy as np
from statsmodels.robust.scale import qn_scale

def robf_despiking(x, mfreq, file_length):
    """
    A parallel implementation of robust functionals for despiking high-frequency 
    eddy covariance raw data after Vitale et al. 2020, Biogeosciences, "Robust data cleaning procedure
    for eddy covariance flux measurements"" and associated inst_prob_test
    routine found on Gitlab

    parameters
    ----------
    x :  a vector of the observed time-series values.
    mfreq : the main frequency of the observed time series (24 or 48 for hourly 
            and halfhourly time series, respectively; 10 or 20 for raw EC data 
            acquired at 10Hz or 20Hz, respectively).
    
    file_length : the raw data file length in minutes (e.g. 30, 60 minutes).

    returns
    -------
    despiked_ts: A vector of the despiked time-series.
    spike_loc Integer: Location of the detected spikes.

    Comments
    --------
    Written by B. Heinesch, 2 November, 2022.
    University of Liege, Gembloux Agro-Bio Tech.

    """
    
	# Window width selection
    nL = mfreq*60*file_length
    x = x[1:nL]
    mod_rlm = rlm(as.vector(x)~poly(seq(1,nL,1), 5, raw = FALSE))
	ind = np.where((abs(as.vector(residuals(mod_rlm))) > 3*qn_scale(as.vector(residuals(mod_rlm))))[0]
	out_patch1 = statsNA(replace(x, ind, NA), print_only=FALSE)$longest_na_gap
	out_patch2 = max(rollapply(replace(x, ind, NA), width=mfreq*30+1, by=1, function(x) length(which(is.na(x)))))
	wl0 = max(c(mfreq*5+1, out_patch1*4, out_patch2*4), na.rm=TRUE)
	wl = wl0 if odd(wl0) else wl0+1
	
	if(wl > mfreq*60+1):
		dspk = despiking(x, mfreq=mfreq, variant="v3", wsignal=wl, wscale=wl, wby=1, zth=5)
		spike_loc = dspk$spike_loc
		ts_cleaned = dspk$ts_cleaned
		
	
	if(wl <= mfreq*60+1):
		mm = matrix(x, nrow=nL/6)
		plan(multicore)
		dspk = future_apply(mm, MARGIN=2, function(x) despiking(x, mfreq=mfreq, variant="v3", wsignal=wl, wscale=wl, wby=1, zth=5))
		spike_loc = na.omit(as.vector(c(dspk[[1]]$spike_loc, dspk[[2]]$spike_loc + mfreq*60*5, dspk[[3]]$spike_loc + mfreq*60*5*2, dspk[[4]]$spike_loc + mfreq*60*5*3, dspk[[5]]$spike_loc + mfreq*60*5*4, dspk[[6]]$spike_loc + mfreq*60*5*5)))
		ts_cleaned = as.vector(c(dspk[[1]]$ts_cleaned, dspk[[2]]$ts_cleaned,dspk[[3]]$ts_cleaned, dspk[[4]]$ts_cleaned, dspk[[5]]$ts_cleaned, dspk[[6]]$ts_cleaned))

        return [ts_cleaned, spike_loc]


