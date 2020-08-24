# roses_array.py
#
# Python module containing array and network methods for ROSES 2020
#
# Stephen Arrowsmith (sarrowsmith@smu.edu)

import numpy as np
from obspy.signal.array_analysis import *
from obspy.core import AttribDict
from pyproj import Geod
from numpy.linalg import inv
import matplotlib.pyplot as plt
import warnings, utm
warnings.filterwarnings("ignore")

def plot_sliding_window(st, T, B, V, C=None, v_max=5., element='BEAM', twin_plot=None):
    '''
    Plots the results of sliding-window array processing
    
    Inputs:
    st - ObsPy Stream object containing array data
    T - Times of array processing estimates (center of time windows) (s)
    B - Backazimuths
    V - Trace velocities (km/s)
    C - Optional color of points (e.g., Semblance, F-statistic, Correlation)
    v_max - Maximum trace velocity for y-axis on trace velocity
    element - Name of the element to plot the time series data for
    twin_plot - List containing the start and end times (in seconds) to plot
    '''
    
    tr = st.select(station=element)[0]
    
    ax1 = plt.subplot(3,1,1)
    t_tr = np.arange(0, tr.stats.npts*tr.stats.delta, tr.stats.delta)
    plt.plot(t_tr, tr.data/np.max(np.abs(tr.data)), 'k-')
    
    ax2 = plt.subplot(3,1,2, sharex=ax1)
    if C is not None:
        plt.scatter(T, B, s=4, c=C, cmap=plt.get_cmap('hot_r'))
    else:
        plt.plot(T, B, 'k.')
    ax2.set_ylim([0,360])
    ax2.set_ylabel('Backazimuth')
    if twin_plot is not None:
        plt.xlim(twin_plot)
    
    ax3 = plt.subplot(3,1,3, sharex=ax1)
    if C is not None:
        plt.scatter(T, V, s=4, c=C, cmap=plt.get_cmap('hot_r'))
    else:
        plt.plot(T, V, 'k.')
    ax3.set_ylim([0,v_max])
    ax3.set_ylabel('Phase vel.')
    ax3.set_xlabel('Time (s) after ' + str(tr.stats.starttime).split('.')[0].replace('T', ' '))
    
    ax1.get_yaxis().set_ticks([])

def sliding_time_array_fk(st, tstart, tend, win_len=20, win_frac=0.5, frqlow=0.5, frqhigh=4,
                          sll_x=-3.6, slm_x=3.6, sll_y=-3.6, slm_y=3.6, sl_s=0.18):
    '''
    Processes st_fk with sliding window FK analysis. Default parameters are suitable for most
    regional infrasound arrays.

    Returns a numpy array with timestamp, relative power, absolute power, backazimuth, slowness
    '''
     
    for st_i in st:
        st_i.stats.coordinates = AttribDict({
            'latitude': st_i.stats.sac.stla,
            'elevation': st_i.stats.sac.stel,
            'longitude': st_i.stats.sac.stlo})

    kwargs = dict(
            # slowness grid: X min, X max, Y min, Y max, Slow Step
            sll_x=sll_x, slm_x=slm_x, sll_y=sll_y, slm_y=slm_y, sl_s=sl_s,
            # sliding window properties
            win_len=win_len, win_frac=win_frac,
            # frequency properties
            frqlow=frqlow, frqhigh=frqhigh, prewhiten=0,
            # restrict output
            semb_thres=-1e9, vel_thres=-1e9, timestamp='mlabday',
            stime=st[0].stats.starttime+tstart, etime=st[0].stats.starttime+tend, verbose=False
        )

    out = array_processing(st, **kwargs)

    return out

def sliding_time_array_lsq(st, X, tstart, tend, twin, overlap):
    '''
    Performs sliding time-window array processing using the least-squares array processing
    method in array_lsq
    
    Inputs:
    st - ObsPy Stream object containing array data
    X - array coordinates
    tstart - Start time for processing (in seconds after Stream start-time)
    tend - End time for processing (in seconds after Stream start-time)
    twin - Time window for array processing (s)
    overlap - Overlap for array processing (s)
    
    Outputs:
    T - Times of array processing estimates (center of time windows) (s)
    V - Phase velocities
    B - Backazimuths
    
    Stephen Arrowsmith (sarrowsmith@smu.edu)
    '''

    time_start = st[0].stats.starttime + tstart
    time_end = time_start+tend

    time_start_i = time_start
    time_end_i = time_start_i+twin

    t = tstart; T = []; V = []; B = []
    while time_end_i < time_end:
        st_win = st.slice(time_start_i, time_end_i)
        vel, baz = array_lsq(st_win, X)
        T.append(t + twin/2); V.append(vel); B.append(baz)
        t = t + overlap
        time_start_i = time_start_i + overlap
        time_end_i = time_end_i + overlap
    T = np.array(T); V = np.array(V); B = np.array(B)
    
    return T, V, B

def array_lsq(st, X):
    '''
    Performs pairwise cross-correlation on each trace in st, and least-squares inversion
    for the slowness vector corresponding to the best-fitting plane wave
    
    Inputs:
    st - ObsPy Stream object containing array data (in a time window suitable for array processing)
    X - [Nx2] NumPy array of array coordinates (in km relative to a reference point)
    
    Outputs:
    baz - Backazimuth (in degrees from North)
    vel - Phase velocity (in km/s)
    
    Notes:
    This function requires that all traces in st begin at the same time (within the sampling interval)
    
    Stephen Arrowsmith (sarrowsmith@smu.edu)
    '''
    
    # Initializing empty arrays for array distances and delay times:
    N = len(st)           # Number of elements
    M = int(N*(N-1)/2)    # Number of pairs of elements
    R = np.zeros((M,2))   # Array to hold relative coordinates between elements
    tau = np.zeros((M,1)) # Array to hold delay times

    k = 0
    for i in range(0,N):
        for j in range(i+1,N):

            tr1 = st[i]; tr2 = st[j]
            C = np.correlate(tr1.data, tr2.data, mode='full')
            lags = np.arange(-np.floor(len(C)/2), np.floor(len(C)/2)+1, 1)*tr1.stats.delta

            # Computing lag corresponding to maximum correlation:
            ix = np.argmax(C); tau[k] = lags[ix]

            # Computing vector of distances between array coordinates:
            R[k,:] = X[i,:] - X[j,:]

            k = k + 1
    
    # Performing least squares inversion:
    R = np.matrix(R); tau = np.matrix(tau)
    u = (inv(np.transpose(R)*R)*np.transpose(R))*tau
    v = 1/np.sqrt(u[0]**2 + u[1]**2)
    azimut = 180 * math.atan2(u[0], u[1]) / math.pi
    baz = (azimut % -360 + 180) % 360
    
    return float(v), float(baz)

def gc_backzimuth(st, evlo, evla):
    '''
    Computes the Great-Circle backazimuth and epicentral distance for an array in st given a known event location
    
    Inputs:
    st - ObsPy Stream object containing array data
    evlo - A float containing the event longitude
    evla - A float containing the event latitude
    
    Outputs:
    baz - Backazimuth (degrees from North)
    dist - Great-circle distance (km)
    '''
    
    g = Geod(ellps='sphere')
    
    a12,a21,dist = g.inv(evlo,evla,st[0].stats.sac.stlo,st[0].stats.sac.stla); dist = dist/1000.
    
    return a21%360., dist

def plot_array_coords(X, stnm):
    '''
    Plots the array coordinates for a given array with coordinates X and element names stnm
    '''
    
    plt.plot(X[:,0], X[:,1], '.')
    for i in range(0, len(stnm)):
        plt.text(X[i,0], X[i,1], stnm[i])
    plt.xlabel('X (km)')
    plt.ylabel('Y (km)')

def get_array_coords(st, ref_station):
    '''
    Returns the array coordinates for an array, in km with respect to the reference array provided
    
    Inputs:
    st - ObsPy Stream object containing array data
    ref_station - A String containing the name of the reference station
    
    Outputs:
    X - [Nx2] NumPy array of array coordinates in km
    stnm - [Nx1] list of element names
    
    Stephen Arrowsmith (sarrowsmith@smu.edu)
    '''
    
    X = np.zeros((len(st), 2))
    stnm = []
    for i in range(0, len(st)):
        E, N, _, _ = utm.from_latlon(st[i].stats.sac.stla, st[i].stats.sac.stlo)
        X[i,0] = E; X[i,1] = N
        stnm.append(st[i].stats.station)

    # Adjusting to the reference station, and converting to km:
    ref_station_ix = np.where(np.array(stnm) == ref_station)[0][0]    # index of reference station
    X[:,0] = (X[:,0] - X[ref_station_ix,0])
    X[:,1] = (X[:,1] - X[ref_station_ix,1])
    X = X/1000.
    
    return X, stnm

def compute_tshifts(X, baz, vel):
    '''
    Computes the set of time shifts required to align the waveforms for an array and DOA
    
    Inputs:
    X - [Nx2] NumPy array of array coordinates in km
    baz - Backazimuth of plane wave
    vel - Phase velocity of plane wave
    
    Outputs:
    t_shifts - [Nx1] List of time shifts (in seconds)
    '''
    
    # Compute slowness vector:
    sl_y = np.abs(np.sqrt((1/vel**2)/((np.tan(np.deg2rad(baz)))**2+1)))
    sl_x = np.abs(sl_y * np.tan(np.deg2rad(baz)))
    if baz > 180:
        sl_x = -sl_x
    if (baz > 90) and (baz < 270):
        sl_y = -sl_y

    # Computes the time shifts for the slowness vector defined by sl_x, sl_y:
    t_shifts = []
    for i in range(0, X.shape[0]):
        t_shift = X[i,0]*sl_x + X[i,1]*sl_y
        t_shifts.append(t_shift)
        #st_copy[i].stats.starttime = st_copy[i].stats.starttime + t_shift
    
    return t_shifts

def plot_beam(st, int_shifts, ref_station, twin=None, trace_spacing=2., return_beam=False):
    '''
    Plots array beam and waveforms on each array element
    
    Inputs:
    st - ObsPy Stream object containing array data 
    int_shifts - List or NumPy array containing integer time shifts for each element
    ref_station - A String describing the reference station (which all shifts are relative to)
    twin - A list defining the time window for plotting (in seconds from start), e.g., twin=[0,100]
    trace_spacing - A float that controls vertical spacing between traces
    return_beam - Boolean to optionally return the array beam
    
    Notes:
    This function requires that all traces in st begin at the same time (within the sampling interval)
    
    Stephen Arrowsmith (sarrowsmith@smu.edu)
    '''
    
    tr_ref = st.select(station=ref_station)[0]
    
    ix = 0
    for tr in st:

        # Plotting data:
        if int_shifts[ix] < 0:    # Need to truncate data before plotting
            ix_start = np.abs(int_shifts[ix])
            data = tr.data[ix_start::]
        elif int_shifts[ix] > 0:  # Need to add zeros to the beginning before plotting:
            data = np.concatenate((np.zeros((np.abs(int_shifts[ix]))), tr.data))
        else:
            data = tr.data

        # Computing stack:
        if ix == 0:
            stack = data
        else:
            diff = len(stack)-len(data)
            if diff > 0:
                data = np.concatenate((data, np.zeros(diff)))
            elif diff < 0:
                stack = np.concatenate((stack, np.zeros(np.abs(diff))))
            stack = data + stack
        
        # Normalizing data prior to plotting
        data = data/(np.max(np.abs(data)))

        plt.plot(np.arange(0, len(data)*tr_ref.stats.delta, tr_ref.stats.delta), ix*trace_spacing + data, 'k')
        
        ix = ix + 1

    beam = stack/len(st)
    
    stack = stack/np.max(np.abs(stack))
    plt.plot(np.arange(0, len(data)*tr_ref.stats.delta, tr_ref.stats.delta), ix*trace_spacing + stack, 'r')
    plt.xlabel('Time (s) after ' + str(tr_ref.stats.starttime).split('.')[0].replace('T', ' '))
    plt.gca().get_yaxis().set_ticks([])
    
    if twin is not None:
        plt.gca().set_xlim([twin[0],twin[1]])
    
    if return_beam:
        return beam

def add_beam_to_stream(st, beam, ref_station):
    '''
    Adds the beam channel to an ObsPy Stream using the time of the reference station
    
    Inputs:
    st - ObsPy Stream object containing array data
    beam - NumPy array containing beam data
    ref_station - String containing the name of the reference station
    
    Outputs:
    st - ObsPy Stream object including beam channel with station name = 'BEAM'
    '''
    
    # Obtain trace for reference station:
    tr = st.select(station=ref_station)[0].copy()
    tr.data = beam[0:len(tr.data)]
    tr.stats.station = 'BEAM'
    st.append(tr)
    
    return st

def plotFK(st, startTime, endTime, frqlow, frqhigh,
           sll_x=-3.6, slm_x=3.6, sll_y=-3.6, slm_y=3.6, sl_s=0.18,
           plot=True, normalize=True, sl_corr=[0.,0.]):
    '''
    Computes and displays an FK plot for an ObsPy Stream object, st, given
    a start time and end time (as UTCDateTime objects) and a frequency band
    defined by frqlow and frqhigh. The slowness grid is defined as optional
    parameters (in s/km).

    This function implements code directly from ObsPy, which has been optimized,
    for simply plotting the FK spectrum

    It includes the option to normalize the data in the time window before running FK

    It also includes the option to apply a slowness correction, defined by sl_corr
    '''

    stream = st.copy()
    stream = stream.trim(startTime, endTime)

    if normalize:
        for st_i in stream:
            st_i.data = st_i.data/np.max(np.abs(st_i.data))
    
    for st_i in stream:
        st_i.stats.coordinates = AttribDict({
            'latitude': st_i.stats.sac.stla,
            'elevation': st_i.stats.sac.stel,
            'longitude': st_i.stats.sac.stlo})

    verbose = False
    coordsys = 'lonlat'
    method = 0

    prewhiten = 0

    grdpts_x = int(((slm_x - sll_x) / sl_s + 0.5) + 1)
    grdpts_y = int(((slm_y - sll_y) / sl_s + 0.5) + 1)

    geometry = get_geometry(stream, coordsys=coordsys, verbose=verbose)

    time_shift_table = get_timeshift(geometry, sll_x, sll_y,
                                     sl_s, grdpts_x, grdpts_y)
    nstat = len(stream)
    fs = stream[0].stats.sampling_rate
    nsamp = stream[0].stats.npts

    # generate plan for rfftr
    nfft = next_pow_2(nsamp)
    deltaf = fs / float(nfft)
    nlow = int(frqlow / float(deltaf) + 0.5)
    nhigh = int(frqhigh / float(deltaf) + 0.5)
    nlow = max(1, nlow)  # avoid using the offset
    nhigh = min(nfft // 2 - 1, nhigh)  # avoid using nyquist
    nf = nhigh - nlow + 1  # include upper and lower frequency

    # to speed up the routine a bit we estimate all steering vectors in advance
    steer = np.empty((nf, grdpts_x, grdpts_y, nstat), dtype=np.complex128)
    clibsignal.calcSteer(nstat, grdpts_x, grdpts_y, nf, nlow,
                         deltaf, time_shift_table, steer)
    _r = np.empty((nf, nstat, nstat), dtype=np.complex128)
    ft = np.empty((nstat, nf), dtype=np.complex128)

    # 0.22 matches 0.2 of historical C bbfk.c
    tap = cosine_taper(nsamp, p=0.22)
    relpow_map = np.empty((grdpts_x, grdpts_y), dtype=np.float64)
    abspow_map = np.empty((grdpts_x, grdpts_y), dtype=np.float64)

    for i, tr in enumerate(stream):
        dat = tr.data
        dat = (dat - dat.mean()) * tap
        ft[i, :] = np.fft.rfft(dat, nfft)[nlow:nlow + nf]

    ft = np.ascontiguousarray(ft, np.complex128)
    relpow_map.fill(0.)
    abspow_map.fill(0.)

    # computing the covariances of the signal at different receivers
    dpow = 0.
    for i in range(nstat):
        for j in range(i, nstat):
            _r[:, i, j] = ft[i, :] * ft[j, :].conj()
            if i != j:
                _r[:, j, i] = _r[:, i, j].conjugate()
            else:
                dpow += np.abs(_r[:, i, j].sum())
    dpow *= nstat

    clibsignal.generalizedBeamformer(
        relpow_map, abspow_map, steer, _r, nstat, prewhiten,
        grdpts_x, grdpts_y, nf, dpow, method)

    ix, iy = np.unravel_index(relpow_map.argmax(), relpow_map.shape)

    # here we compute baz, slow
    slow_x = sll_x + ix * sl_s
    slow_y = sll_y + iy * sl_s

    # ---------
    slow_x = slow_x - sl_corr[0]
    slow_y = slow_y - sl_corr[1]
    #print(slow_x, slow_y)
    # ---------

    slow = np.sqrt(slow_x ** 2 + slow_y ** 2)
    if slow < 1e-8:
        slow = 1e-8
    azimut = 180 * math.atan2(slow_x, slow_y) / math.pi
    baz = azimut % -360 + 180

    if plot:
        plt.pcolormesh(np.arange(sll_x, slm_x + sl_s, sl_s)+sl_corr[0],
                       np.arange(sll_x, slm_x + sl_s, sl_s)+sl_corr[1],
                       np.flipud(np.fliplr(relpow_map.transpose())))
        plt.xlim(sll_x,slm_x)
        plt.ylim(sll_y,slm_y)
        plt.plot(0, 0, 'w+')
        plt.xlabel('Slowness (x) s/km')
        plt.ylabel('Slowness (y) s/km')
        plt.title('Peak semblance at ' + str(round(baz % 360., 2)) + ' degrees ' + str(round(1/slow, 2)) + ' km/s')

    # only flipping left-right, when using imshow to plot the matrix is takes points top to bottom
    # points are now starting at top-left in row major
    return np.fliplr(relpow_map.transpose()), baz % 360, 1. / slow