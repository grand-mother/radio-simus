import glob
import numpy as np
from scipy.signal import butter, lfilter
from scipy.fftpack import rfft, irfft, rfftfreq

__all__ = ["add_noise", "digitization", "filter"]

##########################################################################
##########################################################################


def add_noise(vrms, voltages):
    """Add normal random noise on voltages

    inputs : (voltage noise rms, voltages)
    outputs : noisy voltages
    """
    voltages[:, 1:] = voltages[:, 1:] + \
        np.random.normal(0, vrms, size=np.shape(voltages[:, 1:]))
    return voltages


def digitization(voltages, tstep, tsampling):
    """Digitize the voltages at an specific sampling

    inputs : (voltages, time step, sampling time)
    outputs : digitized voltages
    """
    ratio = int(round(tsampling / tstep))
    ind = np.arange(0, len(voltages[:, 0]), ratio)
    vf = voltages[ind, 1:]
    tf = voltages[ind, 0]
    return np.array([tf, vf[:, 0], vf[:, 1], vf[:, 2]]).T

# Filter the voltages at a bandwidth
################################################################


def _butter_bandpass(lowcut, highcut, fs, order):
    """subfunction of filt
    """
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def _butter_bandpass_filter(data, lowcut, highcut, fs, order):
    """subfunction of filt
    """
    b, a = butter_bandpass(lowcut, highcut, fs, order)
    y = lfilter(b, a, data)
    return y


def _Filtering(x, fs, lowcut, highcut):
    """subfunction of filt
    """
    return butter_bandpass_filter(x, lowcut, highcut, fs, order=5)


def filter(f=None, d=None, FREQMIN=50.e6, FREQMAX=200.e6):
    """Filters signals at specific bandwith (by default 50-200MHz)

    inputs : (file of voltages or array of voltages, frequency min, frequency max)
    outputs : filtered voltages
    """

    if d is None:  # No array of data
        if f is None:  # No file is given
            raise RuntimeError("No data input! Abort.")
        else:  # Load file
            d = np.loadtxt(f)
    # Now filter
    t = d[:, 0]
    dt = np.mean(np.diff(t))  # Compute time step
    fs = 1 / dt
    v = np.array(d[:, 1:])  # Raw signal
    nCh = np.shape(v)[1]
    vout = np.zeros(shape=(len(t), nCh))
    res = []
    for i in range(nCh):
        vi = v[:, i]
        vout[:, i] = Filtering(vi, fs, FREQMIN, FREQMAX)
        imax = np.argmax(vout[:, i], axis=0)
    for i in range(nCh):
        vi = v[:, i]
        vout[:, i] = Filtering(vi, fs, FREQMIN, FREQMAX)
        imax = np.argmax(vout[:, i], axis=0)
        imin = np.argmin(vout[:, i], axis=0)
        res = res + [t[imax], vout[imax, i], t[imin], vout[imin, i]]

    return np.array([t, vout[:, 0], vout[:, 1], vout[:, 2]]).T
##########################################################################
##########################################################################
