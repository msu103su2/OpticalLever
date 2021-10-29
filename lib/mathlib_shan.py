import numpy as np
from scipy.signal import butter, lfilter

def Single_sided_fft(timeSignal, window = None):
    N = len(timeSignal)
    if window is None:
        fSignal = np.fft.fft(timeSignal, window)
    else:
        fSignal = np.fft.fft(np.multiply(timeSignal, window))
    if N % 2 == 0:
        fSignal = fSignal[0:int(N/2)+1]
    else:
        fSignal = fSignal[0:int((N+1)/2)]
    return fSignal

def Single_sided_PSD_arbU(timeSignal, window = None):
    N = len(timeSignal)
    fSignal = Single_sided_fft(timeSignal, window)
    PSD = np.square(np.absolute(fSignal))
    return PSD

def Single_sided_PSD(timeSignal, fs, window = None):
    N = len(timeSignal)
    fSignal = Single_sided_fft(timeSignal, window)
    PSD = 2*np.square(np.absolute(fSignal))/(fs*N) #PSD in unit of /Hz
    return PSD

def Single_sided_fft_chirp_avg(timeSignal, avg):
    N = len(timeSignal)
    assert N % avg == 0
    n = int(N/avg)
    ts_avg = np.zeros(n)
    for i in range(avg):
         ts_avg = ts_avg + timeSignal[i*n:(i+1)*n]
    ts_avg = ts_avg/avg
    return Single_sided_fft(timeSignal)

def HilbertTransform(timeSignal):
    N = len(timeSignal)
    fSignal = np.fft.fft(timeSignal)
    if N % 2 == 0:
        flo = fSignal[0:int(N/2)+1]
        fhi = fSignal[int(N/2)+1:]
    else:
        flo = fSignal[0:int((N+1)/2)]
        fhi = fSignal[int((N+1)/2):]
    fhi = -fhi
    transformed = np.concatenate((flo, fhi), axis = None)
    transformed  = transformed *(np.complex(0,-1))
    return np.real(np.fft.ifft(transformed))

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

def fitFunc(x, y, fm, offset, Noi, A):
    return offset + 10*np.log10(Noi + A/(np.power((fm**2 - np.power(x, 2)), 2)*(2*np.pi)**4+np.power((2*np.pi*x*y),2)))
