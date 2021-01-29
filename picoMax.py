import lib.redpitaya_scpi as scpi
import lib.Switch as Switch
import lib.AgilentNA as AG
import lib.ps5000A as ps
import lib.splitPair as sp
import lib.FuncGen33522b as FuncGen
import pyvisa as py
import numpy as np
from scipy.optimize import minimize
from scipy import signal
import time
import os
import math
import multiprocessing as mp
from multiprocessing.managers import BaseManager

def PDdualAC_balance(data, fs, dir_filename = None):
    inspectF = 4e5
    span = 4e5
    N = len(data[:,0])
    T = N/fs

    fftC = np.fft.fft(data[:,2])
    if N % 2 == 0:
        f = np.arange(N/2+1)/T
        fftC = fftC[0:int(N/2)+1]
    else:
        f = np.arange((N+1)/2)/T
        fftC = fftC[0:int((N+1)/2)]
    PSD = 2*np.square(np.absolute(fftC))/(T*fs*fs*N/fs) #PSD in unit of /bin

    i1 = np.argmin(abs(f-inspectF+span))
    i2 = np.argmin(abs(f-inspectF-span))
    p = 10*np.log10(PSD[i1:i2]*20)
    loc, _ = signal.find_peaks(p, prominence=50)

    fftA = np.fft.fft(data[:,0])
    fftB = np.fft.fft(data[:,1])

    if N % 2 == 0:
        fftA = fftA[0:int(N/2)+1]
        fftB = fftB[0:int(N/2)+1]
    else:
        fftA = fftA[0:int((N+1)/2)]
        fftB = fftB[0:int((N+1)/2)]

    data = np.stack((fftA[loc+i1], fftB[loc+i1]), 1)
    res = minimize(meritFunc, 1, (data, ), bounds=[(0.9,1.1)])

    if res.success:
        pbest = res.x[0]
        balancedPSD = 2*np.square(np.absolute(fftA*pbest-fftB*(2-pbest)))/(T*fs*fs*N/fs)
        return pbest, balancedPSD
    else:
        return 0, 0

def meritFunc(p, data):
    VA = data[:,0]*p
    VB = data[:,1]*(2-p)
    return np.sum(np.square(np.absolute(VA-VB)))

def func(lock, pi, idleCores, count, ps5000a, GS, fs, wd, data_cache, dataBlockIdx, dataBlockSize):
    PSDoverP = np.array([])
    PSDoverPP = np.array([])
    PSDoverP_sq = np.array([])
    PSDoverPP_sq = np.array([])
    p_count = 0
    cPowers = []
    indexs = []

    for idx in dataBlockIdx:
        data = ps5000a.getData([idx * dataBlockSize, (idx+1) * dataBlockSize])
        temp, newPSD = PDdualAC_balance(data, fs)
        if not temp == 0:
            P = np.mean(data[:,2])
            if p_count == 0:
                PSDoverP = newPSD/P
                PSDoverPP = newPSD/(P*P)
                PSDoverP_sq = np.square(newPSD/P)
                PSDoverPP_sq =  np.square(newPSD/(P*P))
                cPowers.append(P)
                indexs.append(idx)
            else:
                PSDoverP += newPSD/P
                PSDoverPP += newPSD/(P*P)
                PSDoverP_sq += np.square(newPSD/P)
                PSDoverPP_sq +=  np.square(newPSD/(P*P))
                cPowers.append(P)
                indexs.append(idx)
            p_count += 1
    data_cache.set(pi, p_count, PSDoverP, PSDoverPP, PSDoverP_sq, PSDoverPP_sq)
    data_cache.set_indexs(pi, np.array(indexs))
    data_cache.set_cPowers(pi, np.array(cPowers))
    print(indexs, cPowers)
    idleCores.put(pi)

class dataCollector():
    def __init__(self, num_processes, wd, GS):
        self.PSDoverPs = [np.array([None])] * num_processes
        self.PSDoverPPs = [np.array([None])] * num_processes
        self.PSDoverP_sqs = [np.array([None])] * num_processes
        self.PSDoverPP_sqs = [np.array([None])] * num_processes
        self.counts = [0] * num_processes
        self.wd = wd
        self.GS = GS
        self.cPowers = [np.array([None])] * num_processes
        self.indexs = [np.array([None])] * num_processes

    def tofile(self):
        count = sum(self.counts)

        PSD = self.PSDoverPs[0]
        for i in range(1, len(self.counts)):
            PSD += self.PSDoverPs[i]
        PSD = PSD/count
        filename = 'PSDoverP_GS='+str(self.GS)
        PSD.tofile(self.wd+filename+'.bin', sep = '')
        PSDoverP = PSD

        PSD = self.PSDoverPPs[0]
        for i in range(1, len(self.counts)):
            PSD += self.PSDoverPPs[i]
        PSD = PSD/count
        filename = 'PSDoverPP_GS='+str(self.GS)
        PSD.tofile(self.wd+filename+'.bin', sep = '')
        PSDoverPP = PSD

        PSD = self.PSDoverP_sqs[0]
        for i in range(1, len(self.counts)):
            PSD += self.PSDoverP_sqs[i]
        PSD = np.sqrt(PSD/count - np.square(PSDoverP))/np.sqrt(count)
        filename = 'PSDoverP_std_GS='+str(self.GS)
        PSD.tofile(self.wd+filename+'.bin', sep = '')

        PSD = self.PSDoverPP_sqs[0]
        for i in range(1, len(self.counts)):
            PSD += self.PSDoverPP_sqs[i]
        PSD = np.sqrt(PSD/count - np.square(PSDoverPP))/np.sqrt(count)
        filename = 'PSDoverPP_std_GS='+str(self.GS)
        PSD.tofile(self.wd+filename+'.bin', sep = '')

    def set(self, i, count, PSDoverP, PSDoverPP, PSDoverP_sq, PSDoverPP_sq):
        if self.counts[i] == 0:
            self.PSDoverPs[i] = PSDoverP
            self.PSDoverPPs[i] = PSDoverPP
            self.PSDoverP_sqs[i] = PSDoverP_sq
            self.PSDoverPP_sqs[i] = PSDoverPP_sq
            self.counts[i] = count
        else:
            self.PSDoverPs[i] += PSDoverP
            self.PSDoverPPs[i] += PSDoverPP
            self.PSDoverP_sqs[i] += PSDoverP_sq
            self.PSDoverPP_sqs[i] += PSDoverPP_sq
            self.counts[i] += count

    def setGS(self, GS):
        self.GS = GS

    def set_cPowers(self, i, data):
        if len(self.cPowers[i]) == 0:
            self.cPowers[i] = data
        else:
            self.cPowers[i] = np.append(self.cPowers[i], data)

    def set_indexs(self, i, indexs):
        if len(self.indexs[i]) == 0:
            self.indexs[i] = indexs
        else:
            self.indexs[i] = np.append(self.indexs[i], indexs)

    def get_cPowers(self):
        N = 0
        print(self.indexs[3].shape)
        for i in range(len(self.indexs)):
            print(self.indexs[i])
            N = max(N, max(self.indexs[i]))
        temp = np.zeros(N+1)
        print(self.cPowers[3].shape)
        for i in range(len(self.cPowers)):
            temp[self.indexs[i]] = self.cPowers[i]
        return temp[temp != 0]

if __name__ == '__main__':
    avg = 130
    steps = 2
    wd = r'Z:\data\optical lever project\NORCADA_NX53515C\46-test'
    rpip = '192.168.137.14'
    BaseManager.register('picoscope', ps.picoscope)
    BaseManager.register('dataCollector', dataCollector)
    BM = BaseManager()
    BM.start()

    rm = py.ResourceManager()
    FG = FuncGen.FuncGen()
    #NA = AG.AgilentNA(wd)
    #RP = scpi.scpi(rpip)
    DS = sp.splitPair()
    #zx80 = Switch.zx80(RP)

    Ax = DS.A.Position();
    Bx = DS.B.Position();

    ps5000a = BM.picoscope(wd)
    ps5000a.defaultSetting()
    ps5000a.AC('A')
    ps5000a.AC('B')
    ps5000a.DC('C')

    ps5000a.AutoRange('A')
    ps5000a.AutoRange('B')
    ps5000a.AutoRange('C')
    ps5000a.configurePSD(8.515, 4e6)

    wd = wd+'\\'
    if not os.path.exists(wd):
        os.makedirs(wd)

    for j in range(steps):
        # split mirror balancing process
        ps5000a.DC('A')
        ps5000a.DC('B')
        ps5000a.AutoRange('A')
        ps5000a.AutoRange('B')

        ps5000a.configurePSD(8.515, 4e6)
        time.sleep(0.5)
        ps5000a.getTimeSignal();
        data = ps5000a.getData()
        data = data[:,0] - data[:,1]
        unbalance = np.mean(data)
        std = np.std(data)

        flag = True
        if abs(unbalance)>0.5:
            step = 0.005
        else:
            step = 0.002

        while abs(unbalance)>0.3*std:
            if unbalance>0:
                if flag:
                    DS.B.MoveBy(-step)
                else:
                    DS.A.MoveBy(step)
            else:
                if not flag:
                    DS.B.MoveBy(step)
                else:
                    DS.A.MoveBy(-step)
            ps5000a.AutoRange('A')
            ps5000a.AutoRange('B')
            ps5000a.getTimeSignal();
            data = ps5000a.getData()
            diff = data[:,0] - data[:,1]
            newUnbalance = np.mean(diff)
            std = np.std(diff)
            print(DS.getQuasi_gapsize(),' ',newUnbalance)
            print(np.mean(data[:,0]),' ',np.mean(data[:,1]), ' ',np.mean(data[:,2]))
            if newUnbalance*unbalance<0:
                flag = not flag
            unbalance = newUnbalance

        ps5000a.ChangeRangeTo('A', 100)
        ps5000a.ChangeRangeTo('B', 100)
        ps5000a.AC('A')
        ps5000a.AC('B')
        nTimes, dataBlockSize = ps5000a.configurePSD(8.515, 4e6, 'Max')

        GS = DS.getQuasi_gapsize()
        fs = ps5000a.getfs()

        cores_touse = mp.cpu_count()
        lock = mp.Lock()
        idleCores = mp.Queue()
        count = mp.Value('i', 0)
        data_cache = BM.dataCollector(cores_touse, wd, GS)

        remainder = nTimes % cores_touse
        quotient = nTimes // cores_touse
        dataBlockIdx = [[]]*cores_touse
        for i in range(cores_touse):
            dataBlockIdx[i] = list(range(i*quotient, (i+1)*quotient))
            if i < remainder:
                dataBlockIdx[i].append(nTimes - i - 1)
        time.sleep(1)


        for nth in range(math.ceil(avg/nTimes)):
            start = time.time()
            ps5000a.getTimeSignal()
            processes = []
            for i in range(cores_touse):
                processes.append(mp.Process(target = func, args = (lock, i, idleCores, count, ps5000a, GS, fs, wd, data_cache, dataBlockIdx[i], dataBlockSize)))
                processes[i].start()
            for i in range(cores_touse):
                pi = idleCores.get()
                processes[pi].join()
                processes[pi].close()
            print(time.time()-start)

        data_cache.tofile()
        cPowers = data_cache.get_cPowers()
        filename = 'cPower_GS='+str(GS)
        cPowers.tofile(wd+filename+'.bin', sep = '')

        print((j+1)/steps*100)
        DS.OpenGapBy(0.05)

    ps5000a.PSDfromTS(ps5000a.getTimeSignal('A'), ps5000a.getfs())
    nprows = np.array(ps5000a.getfreq())
    filename = 'PSDfreq'
    nprows.tofile(wd+filename+'.bin', sep = '')
    DS.A.MoveTo(Ax)
    DS.B.MoveTo(Bx)
