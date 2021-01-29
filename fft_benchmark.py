from multiprocessing import Process
from multiprocessing.managers import BaseManager
import numpy as np
import time
from scipy.signal import periodogram
from scipy import signal
from scipy.optimize import minimize
#import pyfftw

def PDdualAC_balance(data, fs, freq):
    inspectF = 4e5
    span = 4e5
    N = len(data[:,0])
    T = N/fs

    fftC = np.fft.fft(data[:,2])
    if N % 2 == 0:
        f = np.arange(N/2+1)/T
        fftC = fftC[0:int(N/2)]
    else:
        f = np.arange((N+1)/2)/T
        fftC = fftC[0:int((N-1)/2)]
    PSD = 2*np.square(np.absolute(fftC))/(T*fs*fs*N/fs) #PSD in unit of /bin

    i1 = np.argmin(abs(f-inspectF+span))
    i2 = np.argmin(abs(f-inspectF-span))
    p = 10*np.log10(PSD[i1:i2]*20)
    loc, _ = signal.find_peaks(p, prominence=50)

    fftA = np.fft.fft(data[:,0])
    fftB = np.fft.fft(data[:,1])

    if N % 2 == 0:
        fftA = fftA[0:int(N/2)]
        fftB = fftB[0:int(N/2)]
    else:
        fftA = fftA[0:int((N-1)/2)]
        fftB = fftB[0:int((N-1)/2)]

    data = np.stack((fftA[loc+i1], fftB[loc+i1]), 1)
    res = minimize(meritFunc, 1, (data, ), bounds=[(0.9,1.1)])

    if res.success:
        return res.x
    else:
        return 0

def meritFunc(p, data):
    VA = data[:,0]*p
    VB = data[:,1]*(2-p)
    return np.sum(np.square(np.absolute(VA-VB)))

class test_class():
    def __init__(self):
        self.a = 1
    def func(self):
        return(self.a)

def func(test_inst):
    print(hex(id(test_inst.func())))

if __name__ == '__main__':
    file = r'Z:\data\optical lever project\NORCADA_NX53515C\20-SNR\TS_GS=12.300440899753948_count=999.bin'
    data = np.fromfile(file)
    data = np.reshape(data, [-1,3])
    VA = data[:,0]
    VB = data[:,1]
    VC = data[:,2]

    avg = 10
    start = time.time()
    for i in range(avg):
        np.fft.fft(VA-VB)
        #periodogram(VA-VC)
        #pyfftw.interfaces.numpy_fft.fft(VA-VB)
        #print(PDdualAC_balance(data, 1e8/6.4, range(len(VA))))
    print((time.time()-start)/10)
    '''
    BaseManager.register('test_class', test_class)
    mg = BaseManager()
    mg.start()
    test_inst = mg.test_class()
    test_inst2 = test_class()

    p1 = Process(target = func, args = (test_inst, ))
    p2 = Process(target = func, args = (test_inst, ))
    p1.start()
    p2.start()
    p1.join()
    p2.join()'''
