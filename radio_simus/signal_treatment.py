import numpy as np
from scipy.signal import butter, lfilter, resample
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


def digitization(voltages, tsampling):
    """Digitize the voltages at an specific sampling

    inputs : (voltages, sampling rate)
    outputs : digitized voltages
    """
    dt = voltages[0, 1] - voltages[0, 0]
    num = len(voltages[0, :]) * int(tsampling / dt)
    if tsampling % dt != 0:
        raise ValueError("Sampling time not multiple of simulation time")
    t_sampled = resample(voltages[0, :], num)
    Vx_sampled = resample(voltages[1, :], num)
    Vy_sampled = resample(voltages[2, :], num)
    Vz_sampled = resample(voltages[3, :], num)
    return np.array([t_sampled, Vx_sampled, Vy_sampled, Vz_sampled])

# Filter the voltages at a bandwidth
################################################################


def _butter_bandpass_filter(data, lowcut, highcut, fs):
    """subfunction of filt
    """
    b, a = butter(5, [lowcut / (0.5 * fs), highcut / (0.5 * fs)],
                  btype='band')  # (order, [low, high], btype)
    return lfilter(b, a, data)


def filter(voltages, FREQMIN=50.e6, FREQMAX=200.e6):
    """Filters signals at specific bandwith (by default 50-200MHz)

    inputs : (array of voltages, frequency min, frequency max)
    outputs : filtered voltages
    """
    t = voltages[0, :]
    v = np.array(voltages[1:, :])  # Raw signal

    fs = 1 / np.mean(np.diff(t))  # Compute frequency step
    nCh = np.shape(v)[0]
    vout = np.zeros(shape=(nCh, len(voltages[0, :])))
    res = []
    for i in range(nCh):
        vi = v[:, i]
        vout[:, i] = _butter_bandpass_filter(vi, FREQMIN, FREQMAX, fs)
        imax = np.argmax(vout[:, i], axis=0)
    for i in range(nCh):
        vi = v[:, i]
        vout[:, i] = _butter_bandpass_filter(vi, FREQMIN, FREQMAX, fs)
        imax = np.argmax(vout[:, i], axis=0)
        imin = np.argmin(vout[:, i], axis=0)
        res = res + [t[imax], vout[imax, i], t[imin], vout[imin, i]]

    return np.array([t, vout[0, :], vout[1, :], vout[2, :]])
##########################################################################
##########################################################################
