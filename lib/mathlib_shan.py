import numpy as np
import scipy
from scipy.signal import butter, lfilter
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

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

def fit_logLorentzian(x, y, p0, plot = 1):
    popt, pcov = curve_fit(fitFunc, x, y, p0 =p0)
    if plot:
        plt.plot(x, y, '.')
        fitx = np.linspace(x[0], x[-1], len(x)*10)
        plt.plot(fitx, fitFunc(np.array(fitx), *popt))
        plt.show()
    return popt, pcov

def fitFunc(x, y, fm, offset, Noi, A):
    x = np.array(x, dtype = np.float)
    y = np.float(y)
    fm = np.float(fm)
    offset = np.float(offset)
    Noi = np.float(Noi)
    A = np.float(A)
    return offset + 10*np.log10(Noi + A/(np.power((fm**2 - np.power(x, 2)), 2)*(2*np.pi)**4+np.power((2*np.pi*x*y),2)))

def Exclude_data(xdata, ydata, i0s, i1s):
    x = (np.array(xdata)).copy()
    y = (np.array(ydata)).copy()
    for i in range(len(i0s)):
        x[i0s[i]:i1s[i]] = np.zeros(i1s[i] - i0s[i])
    y = y[np.nonzero(x)]
    x = x[np.nonzero(x)]
    return x, y

def Approximate_index(x, value):
    x = np.array(x)
    return np.argmin(np.absolute(x - value))

def angle_of_two_vector(a, b):
    #a is the reference vector, angle is between -pi and pi
    theta = np.arccos(np.inner(a, b)/np.sqrt(np.inner(a,a)*np.inner(b,b)))
    diff = np.array(b) - np.array(a)
    if np.cross(b, a) > 0:
        return theta
    else:
        return -theta

def angle_diff(angleA, angleB):
    a = [np.cos(angleA), np.sin(angleA)]
    b = [np.cos(angleB), np.sin(angleB)]
    return angle_of_two_vector(a, b)

def HG(x, x0, w, P):
    return P*np.sqrt(2/np.pi)/w*np.exp(-2*(x-x0)**2/w**2)

def CumHG(x, x0, w, P):
    return scipy.integrate.quad(HG, -np.inf, x, args=(x0, w, P))
