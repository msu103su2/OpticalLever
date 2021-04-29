import numpy as np

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

def Single_sided_fft_chirp_avg(timeSignal, avg):
    N = len(timeSignal)
    assert N % avg == 0
    n = int(N/avg)
    ts_avg = np.zeros(n)
    for i in range(avg):
         ts_avg = ts_avg + timeSignal[i*n:(i+1)*n]
    ts_avg = ts_avg/avg
    return Single_sided_fft(timeSignal)
