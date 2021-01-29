import pyvisa as py
import numpy as np
import lib.ps5000A as ps
import time
import os
import multiprocessing as mp
from multiprocessing.managers import BaseManager

def DN_avg(data, fs):
    N = len(data[:,0])
    T = N/fs
    fftBalance = np.fft.fft(data[:,0]-data[:,1])
    if N % 2 == 0:
        fftBalance = fftBalance[0:int(N/2)+1]
    else:
        fftBalance = fftBalance[0:int((N+1)/2)]
    balancedPSD = 2*np.square(np.absolute(fftBalance))/(T*fs*fs*N/fs) #*N/fs turns into unit of dB/bin
    return balancedPSD

def func(lock, pi, idleCores, count, avg, ps5000a, fs, data_cache):
    PSD = np.array([])
    p_count = 0
    while count.value < avg:
        lock.acquire()
        i = count.value
        ps5000a.getTimeSignal(check = False)
        data = ps5000a.getData()
        count.value += 1
        lock.release()

        if i < avg: #extra protection
            if len(PSD) == 0:
                PSD = DN_avg(data, fs)
            else:
                PSD += DN_avg(data, fs)
        p_count += 1
    PSD = PSD / p_count
    data_cache.set(pi, PSD)
    idleCores.put(pi)

class dataCollector():
    def __init__(self, num_processes):
        self.PSDs = [np.array([None])] * num_processes

    def get(self):
        PSD = self.PSDs[0]
        for i in range(1, len(self.PSDs)):
            PSD += self.PSDs[i]
        PSD = PSD/len(self.PSDs)
        return PSD

    def set(self, i, data):
        self.PSDs[i] = data

if __name__ == '__main__':
    avg = 100
    wd = r'Z:\data\optical lever project\NORCADA_NX53515C\31-test'
    BaseManager.register('picoscope', ps.picoscope)
    BaseManager.register('dataCollector', dataCollector)
    BM = BaseManager()
    BM.start()

    ps5000a = BM.picoscope(wd)
    ps5000a.defaultSetting()
    ps5000a.AC('A')
    ps5000a.AC('B')
    ps5000a.DC('C')

    ps5000a.ChangeRangeTo('A', 100)
    ps5000a.ChangeRangeTo('B', 100)
    ps5000a.AutoRange('C')
    ps5000a.configurePSD(8.515, 4e6)

    wd = wd+'\\'
    if not os.path.exists(wd):
        os.makedirs(wd)

    cores_touse = mp.cpu_count()
    lock = mp.Lock()
    idleCores = mp.Queue()
    count = mp.Value('i', 0)
    data_cache = BM.dataCollector(cores_touse)
    fs = ps5000a.getfs()

    start = time.time()
    processes = []
    for i in range(cores_touse):
        processes.append(mp.Process(target = func, args = (lock, i, idleCores, count, avg, ps5000a, fs, data_cache)))
        processes[i].start()
    for i in range(cores_touse):
        pi = idleCores.get()
        processes[pi].join()
        processes[pi].close()

    PSD = data_cache.get()
    filename = 'DN'
    PSD.tofile(wd+filename+'.bin', sep = '')
    print(time.time()-start)

    ps5000a.PSDfromTS(ps5000a.getTimeSignal('A'), ps5000a.getfs())
    nprows = np.array(ps5000a.getfreq())
    filename = 'PSDfreq'
    nprows.tofile(wd+filename+'.bin', sep = '')
