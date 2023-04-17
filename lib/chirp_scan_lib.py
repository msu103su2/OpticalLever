import lib.mathlib_shan as ms
import numpy as np
import os
import re
import time
import matplotlib.pyplot as plt
from scipy import signal
from scipy.optimize import curve_fit, minimize, optimize
import json


def sweep(FG, ps5000a, UC2, Nsteps, AWGF):
    wd = r'Z:\data\optical lever project\NORCADA_NX53515C\90-LineScan'
    if not os.path.exists(wd+'\\'):
        os.makedirs(wd+'\\')
    err = ''
    oldpeakf = 182330
    peakf = oldpeakf
    fspan = 30
    f1 = peakf - fspan/2
    f2 = peakf + fspan/2
    dc = 0.365
    amp = 0.25
    df = 0.25
    avg = 5
    for i in range(Nsteps):
        [tsAs, tsBs, tsCs] = get_ts(df, avg, ps5000a, trigger = 'D', triggerThreshold = 8000)
        UC2.PR(int(1), -10, err)

        offset = -int(np.round(1e3/0.930*(1/df)/ps5000a.timeInternalns.value))
        ps5000a.PSDfromTS(ps5000a.chs['A'].timeSignal[0:offset], ps5000a.getfs())
        window = signal.windows.hann(len(ps5000a.chs['A'].timeSignal[0:offset]))
        check1 = np.argmin(np.absolute(ps5000a.f-182200))
        check2 = np.argmin(np.absolute(ps5000a.f-182600))
        f = ps5000a.f[check1:check2]

        dftAs = np.zeros((avg, len(f)), dtype = np.complex)
        dftBs = np.zeros((avg, len(f)), dtype = np.complex)
        dftCs = np.zeros((avg, len(f)), dtype = np.complex)
        dfts = np.zeros((avg, len(f)), dtype = np.complex)
        for j in range(avg):
            dftA = ms.Single_sided_fft(tsAs[j,:], window)
            if j == 0:
                maxI = np.argmax(np.absolute(dftA[check1: check2]))
                peakf = f[maxI]
                if np.absolute(peakf - oldpeakf)>3:
                    f1 = peakf - fspan/2
                    f2 = peakf + fspan/2
                    if AWGF == 'chirp':
                        FG.chirp('test', f1, f2, AWGfs, AWGN, dc+amp, dc-amp)
                    elif AWGF == 'sweep':
                        FG.sweep(f1, f2, 1/df, amp, dc)
                    oldpeakf = peakf
            dftB = ms.Single_sided_fft(tsBs[j,:], window)
            dftC = ms.Single_sided_fft(tsCs[j,:], window)
            dftAs[j,:] = dftA[check1: check2]
            dftBs[j,:] = dftB[check1: check2]
            dftCs[j,:] = dftC[check1: check2]
            dfts[j,:] = np.multiply(np.divide(dftAs[j,:], dftBs[j,:]), np.absolute(dftBs[j,:]))

        dftavg = np.zeros(len(f), dtype = np.complex)
        drivepsd = np.zeros(len(f))
        for j in range(avg):
            dftavg = dftavg + dfts[j,:]
            drivepsd = drivepsd + np.square(np.absolute(dftCs[j,:]))
        drivepsd = 10*np.log10(20 * drivepsd/avg)
        psdavg = 10*np.log10(20* np.square(np.absolute(dftavg/avg)))
        dftAs.tofile(wd+'\\'+'BalDFT_{step:d}'.format(step = i)+'.bin', sep = '')
        dftBs.tofile(wd+'\\'+'LaserPDFT_{step:d}'.format(step = i)+'.bin', sep = '')
        dftCs.tofile(wd+'\\'+'DriveDFT_{step:d}'.format(step = i)+'.bin', sep = '')
        psdavg.tofile(wd+'\\'+'PSD_{step:d}'.format(step = i)+'.bin', sep = '')
        drivepsd.tofile(wd+'\\'+'DrivePSD_{step:d}'.format(step = i)+'.bin', sep = '')
        print('{percent:%}'.format(percent = i/Nsteps))
    f.tofile(wd+'\\'+'freq'+'.bin', sep = '')
    print(peakf)

def get_chirp(df, avg, ps5000a, trigger = None, triggerThreshold = 0):
    factor = 0.25/df
    ps5000a.configurePSD(0.25/factor, 1953125) #gives sampling intervals of 256ns
    ps5000a.getTimeSignal(trigger = trigger, triggerThreshold = triggerThreshold)
    offset = -int(np.round(1e3/0.930*(1/df)/ps5000a.timeInternalns.value))
    tsAs = np.zeros((avg, len(ps5000a.chs['A'].timeSignal[0:offset])))
    tsBs = np.zeros((avg, len(ps5000a.chs['A'].timeSignal[0:offset])))
    tsCs = np.zeros((avg, len(ps5000a.chs['A'].timeSignal[0:offset])))
    for i in range(avg):
        ps5000a.getTimeSignal(trigger = trigger, triggerThreshold = triggerThreshold)
        tsAs[i,:] = ps5000a.chs['A'].timeSignal[0:offset]
        tsBs[i,:] = ps5000a.chs['B'].timeSignal[0:offset]
        tsCs[i,:] = ps5000a.chs['C'].timeSignal[0:offset]
    return [tsAs, tsBs, tsCs]

def get_ts(df, avg, ps5000a, trigger):
    factor = 0.25/df
    ps5000a.configurePSD(0.25/factor, 1953125) #gives sampling intervals of 256ns
    tsAs = np.zeros((avg, len(ps5000a.chs['A'].timeSignal[0:offset])))
    tsBs = np.zeros((avg, len(ps5000a.chs['A'].timeSignal[0:offset])))
    tsCs = np.zeros((avg, len(ps5000a.chs['A'].timeSignal[0:offset])))
    for i in range(avg):
        ps5000a.getTimeSignal(trigger = trigger)
        tsAs[i,:] = ps5000a.chs['A'].timeSignal[0:offset]
        tsBs[i,:] = ps5000a.chs['B'].timeSignal[0:offset]
        tsCs[i,:] = ps5000a.chs['C'].timeSignal[0:offset]
    return [tsAs, tsBs, tsCs]

def chirp_check(fcs, fspans, df, avg, ps5000a, FG, FGwrite, AWGF, filename = None, wd = None,  amps = 0.25, dc = 0.365, AWGfs = 3.2e6, AWGN = 12.8e6, offset = None):
    fcs = np.array(fcs)
    fspans = np.array(fspans)
    amps = np.array(amps)
    arbdata = np.zeros(int(AWGN))

    if FGwrite:
        if AWGF == 'chirp':
            arbdata = FG.chirps('test', fcs - fspans/2, fcs + fspans/2, AWGfs, AWGN, dc+amps, dc-amps)
        elif AWGF == 'sweep':
            FG.sweep(fcs[0]-fspans[0]/2, fcs[0]-fspan[0]/2, 1/df, amps[0], dc)
    factor = 0.25/df
    ps5000a.configurePSD(0.25/factor, 1953125) #gives sampling intervals of 256ns
    if offset is None:
        offset = -int(np.round(1e3/0.930*(1/df)/ps5000a.timeInternalns.value))
    ps5000a.getTimeSignal()
    temp = ps5000a.PSDfromTS(ps5000a.chs['A'].timeSignal[0:offset], ps5000a.getfs())
    window = signal.windows.tukey(len(ps5000a.chs['A'].timeSignal[0:offset]), 0)
    f = ps5000a.f[0:len(temp)]

    tsAs = np.zeros((avg, len(ps5000a.chs['A'].timeSignal[0:offset])))
    tsBs = np.zeros((avg, len(ps5000a.chs['A'].timeSignal[0:offset])))
    tsCs = np.zeros((avg, len(ps5000a.chs['A'].timeSignal[0:offset])))
    for i in range(avg):
        ps5000a.getTimeSignal(trigger = 'D', triggerThreshold = 8000)
        tsAs[i,:] = ps5000a.chs['A'].timeSignal[0:offset]
        tsBs[i,:] = ps5000a.chs['B'].timeSignal[0:offset]
        tsCs[i,:] = ps5000a.chs['C'].timeSignal[0:offset]
        print(i)

    powers, gamma, psdAs, psdCs, freqs, popts= Piecewise(fcs, fspans, tsAs[i], tsBs[i], f, ps5000a.getfs())
    '''
    if filename is not None:
        dftAs.tofile(wd+'\\'+'BalDFT_{filename:s}'.format(filename = filename)+'.bin', sep = '')
        dftBs.tofile(wd+'\\'+'LaserPDFT_{filename:s}'.format(filename = filename)+'.bin', sep = '')
        dftCs.tofile(wd+'\\'+'DriveDFT_{filename:s}'.format(filename = filename)+'.bin', sep = '')
        psdavg.tofile(wd+'\\'+'PSD_{filename:s}'.format(filename = filename)+'.bin', sep = '')
        drivepsd.tofile(wd+'\\'+'DrivePSD_{filename:s}'.format(filename = filename)+'.bin', sep = '')
        '''

    for psdavg, drivepsd, popt, freq in zip(psdAs, psdCs, popts, freqs):
        plt.plot(freq, psdavg, label='signal')
        plt.plot(freq, drivepsd, label='drive')
        plt.plot(freq, fitFunc(freq, *popt), label='fit')
        plt.plot(freq, psdavg - drivepsd, label='response')
        plt.legend(loc="upper left")
        plt.show()
        input("Press Enter to continue...")

    print(gamma)
    print(powers)


    for j in range(len(psdAs)):
        plt.plot(freqs[0], psdAs[j] + 50*j)
    plt.show()
    return np.array(arbdata)

def window_cp(window, ps5000a, drive, offset):
    s1 = ps5000a.chs['A'].timeSignal
    s2 = ps5000a.chs[drive].timeSignal
    s1 = s1[0:offset]
    s2 = s2[0:offset]
    f,p1 = signal.periodogram(s1, ps5000a.getfs(), window = window)
    f,p2 = signal.periodogram(s2, ps5000a.getfs(), window = window)
    check1 = np.argmin(np.absolute(f-182200+100))
    check2 = np.argmin(np.absolute(f-182400-100))
    f = f[check1:check2]
    p1 = p1[check1:check2]
    p2 = p2[check1:check2]
    plt.plot(f, np.log10(p1))
    plt.plot(f, np.log10(p2))
    plt.show()
    plt.plot(f, np.log10(p1)-np.log10(p2))
    plt.show()

def fitFunc(x, y, fm, offset):
    return offset - 10*np.log10(np.power((fm**2 - np.power(x, 2)), 2)*(2*np.pi)**4+np.power((2*np.pi*x*y),2))

def Piecewise(fcs, fspans, tsA, tsC, f, fs, gamma = None):
    fcs = np.array(fcs)
    fspans = np.array(fspans)
    fstarts = fcs - fspans/2
    fends = fcs + fspans/2
    Nj = len(fcs)
    psdA = 10*np.log10(20*ms.Single_sided_PSD(tsA, fs, window = None))
    psdC = 10*np.log10(20*ms.Single_sided_PSD(tsC, fs, window = None))
    powers = np.zeros(Nj)

    check1s = np.zeros(Nj, dtype = np.int64)
    check2s = np.zeros(Nj, dtype = np.int64)
    for f1, f2, j in zip(fstarts, fends, list(range(Nj))):
        check1s[j] = np.argmin(np.absolute(f-f1+100))
        check2s[j] = np.argmin(np.absolute(f-f2-100))

    psdAs = [[]]*Nj
    psdCs = [[]]*Nj
    freqs = [[]]*Nj

    for j in range(Nj):
        psdAs[j] = psdA[check1s[j]:check2s[j]]
        psdCs[j] = psdC[check1s[j]:check2s[j]]
        freqs[j] = f[check1s[j]:check2s[j]]

    maxIs = np.zeros(Nj, dtype = np.int64)
    popts = [[]]*Nj
    gamma_provided = 0
    if gamma is None:
        gamma = [10]*Nj
        gamma_provided = False
    else:
        gamma_provided = True
    for j in range(Nj):
        maxIs[j] = np.argmax(psdAs[j])
        fitlo = min(300, maxIs[j])
        fithi = min(100,len(psdAs[j] - maxIs[j]))
        if gamma_provided:
            try:
                popts[j], pcov = curve_fit(fitFunc, freqs[j][maxIs[j]-fitlo:maxIs[j]+fithi], psdAs[j][maxIs[j]-fitlo:maxIs[j]+fithi] - psdCs[j][maxIs[j]-fitlo:maxIs[j]+fithi] , p0 =[gamma[j], freqs[j][maxIs[j]],  180], \
                    bounds = ([0.99*gamma[j], freqs[j][maxIs[j]] - 50,  120],[1.01*gamma[j], freqs[j][maxIs[j]]+50,  370]))
            except RuntimeError:
                print('fit failed. 0 given as result')
                popts[j] = [0,0,0]
        else:
            try:
                popts[j], pcov = curve_fit(fitFunc, freqs[j][maxIs[j]-fitlo:maxIs[j]+fithi], psdAs[j][maxIs[j]-fitlo:maxIs[j]+fithi] - psdCs[j][maxIs[j]-fitlo:maxIs[j]+fithi] , p0 =[gamma[j], freqs[j][maxIs[j]],  180], \
                    bounds = ([0, freqs[j][maxIs[j]] - 50,  120],[50, freqs[j][maxIs[j]]+50,  370]))
            except RuntimeError:
                print('fit failed. 0 given as result')
                popts[j] = [0,0,0]
        gamma[j] = popts[j][0]
        powers[j] = popts[j][2]
    return powers, gamma, psdAs, psdCs, freqs, popts
def LineScan(fcs, fspans, UC2, ps5000a, DS, UC2Amp, UC2_relative_xs, UC2_controllerAddress, wd):
    if not os.path.exists(wd):
        os.makedirs(wd)
    err = ''
    UC2.SU_Set(UC2_controllerAddress, '+', UC2Amp, err)
    UC2.SU_Set(UC2_controllerAddress, '-', UC2Amp, err)
    avg = 1
    waitTime_after_moving = 30

    temp = np.zeros(len(UC2_relative_xs)+2)
    temp[1:-1] = UC2_relative_xs
    steps = np.zeros(len(UC2_relative_xs)+1)
    for i in range(len(steps)):
        steps[i] = temp[i+1] - temp[i]

    UC2.PR(UC2_controllerAddress, steps[0], err)
    time.sleep(30)
    ps5000a.configurePSD(8, 1953125)
    DS.AutoBalance(ps5000a, 1)

    ps5000a.configurePSD(1, 1953125)
    ps5000a.AC('A')
    ps5000a.ChangeRangeTo('A', 5000)
    ps5000a.AC('B')
    ps5000a.AutoRange('B')
    ps5000a.DC('C')
    ps5000a.AutoRange('C')

    ps5000a.configurePSD(0.25, 1953125)
    ps5000a.getTimeSignal()
    ps5000a.PSDfromTS(ps5000a.chs['A'].timeSignal, ps5000a.getfs())
    temp, _, _  = get_chirp(0.25, avg, ps5000a)
    temp = ms.Single_sided_fft(temp[0,:])
    f = ps5000a.f[0:len(temp)]
    range_A = 5000
    powers_on_step = [[]]*len(UC2_relative_xs)

    for i in range(len(UC2_relative_xs)):
        ps5000a.ChangeRangeTo('A', range_A)
        cache = np.zeros(len(temp))
        ps5000a.configurePSD(0.25, 1953125)
        [tsAs, tsBs, tsCs] = get_chirp(0.25, avg, ps5000a, trigger = 'D', triggerThreshold = 8000)
        window = signal.windows.tukey(len(tsAs[0,:]), 0.8)

        dftA = np.zeros(len(temp))
        dftC = np.zeros(len(temp))

        t0 = time.time()
        if i != len(UC2_relative_xs)-1:
            UC2.PR(UC2_controllerAddress, steps[i+1], err)

        for j in range(avg):
            dftA = ms.Single_sided_fft(tsAs[j,:], window)
            dftC = ms.Single_sided_fft(tsCs[j,:], window)
            cache = cache + np.multiply(np.divide(dftA, dftC), np.absolute(dftC))
        cache = cache/avg
        cache = np.square(np.absolute(cache))
        psdC = np.square(np.absolute(dftC))
        #cache.tofile(wd+'\\'+'signal_{i:d}.bin'.format(i = i), sep = '')
        #psdC.tofile(wd+'\\'+'drive_{i:d}.bin'.format(i = i), sep = '')
        print('Progress:{p:.2%}'.format(p = (i+1)/(len(UC2_relative_xs))))

        if i == 0:
            powers, gamma, psdAs, psdCs, freqs, popts = Piecewise(fcs, fspans, tsAs[0], tsBs[0], f, ps5000a.getfs(), gamma = None)
        else:
            powers, _, psdAs, psdCs, freqs, popts = Piecewise(fcs, fspans, tsAs[0], tsBs[0], f, ps5000a.getfs(), gamma)
        powers_on_step[i] = powers
        print(powers)

        t1 = time.time()
        time.sleep(max(waitTime_after_moving - (t1 - t0), 0))
        ps5000a.configurePSD(8, 1953125)
        DS.AutoBalance(ps5000a)

    powers_on_step = np.array(powers_on_step)
    for i in range(len(powers_on_step[0])):
        plt.plot(UC2_relative_xs, powers_on_step[:,i], label=str(fcs[i]))
    plt.legend(loc="upper left")
    plt.show()
    ps5000a.f.tofile(wd+'\\'+'freq.bin', sep = '')
    (np.array(UC2_relative_xs)).tofile(wd+'\\'+'UC2steps.bin', sep = '')
    (np.array(powers_on_step)).tofile(wd+'\\'+'powers.bin', sep = '')
    ps5000a.configurePSD(0.25, 1953125)
def LineScan_sp(fcs, fspans, ps5000a, DS, sp_relative_xs, wd, avg = 1):
    if not os.path.exists(wd):
        os.makedirs(wd)
    err = ''
    waitTime_after_moving = 30

    Nsteps = len(sp_relative_xs)

    temp = np.zeros(Nsteps+2)
    temp[1:-1] = sp_relative_xs
    dxs = np.zeros(Nsteps+1)
    for i in range(Nsteps):
        dxs[i] = temp[i+1] - temp[i]

    BalancedAx = DS.A.x
    BalancedBx = DS.B.x

    DS.MoveGapBy(dxs[0])

    ps5000a.configurePSD(1, 1953125)
    ps5000a.AC('A')
    ps5000a.ChangeRangeTo('A', 1000)
    ps5000a.DC('B')
    ps5000a.ChangeRangeTo('B', 2000)
    ps5000a.DC('C')
    ps5000a.AutoRange('C')

    ps5000a.configurePSD(0.25, 1953125)
    ps5000a.getTimeSignal()
    ps5000a.PSDfromTS(ps5000a.chs['A'].timeSignal, ps5000a.getfs())
    temp, _, _  = get_chirp(0.25, 1, ps5000a)
    temp = ms.Single_sided_fft(temp[0,:])
    f = ps5000a.f[0:len(temp)]
    range_A = 2000
    powers_on_step = [[]]*Nsteps
    drives_on_step = [[]]*Nsteps
    Axs = [0]*Nsteps
    Bxs = [0]*Nsteps
    spdxs = [0]*Nsteps
    spoffs = [0]*Nsteps
    DC_A = [0]*Nsteps
    DC_B = [0]*Nsteps

    for i in range(len(sp_relative_xs)):
        #ps5000a.ChangeRangeTo('A', range_A)
        ps5000a.AutoRange('A')
        cache = np.zeros(len(temp))
        ps5000a.configurePSD(0.25, 1953125)
        [tsAs, tsBs, tsCs] = get_chirp(0.25, avg, ps5000a, trigger = 'D', triggerThreshold = 8000)
        window = signal.windows.tukey(len(tsAs[0,:]), 0)

        cache = np.zeros(len(temp))
        cache_drive = np.zeros(len(temp))

        Axs[i] = DS.A.x
        Bxs[i] = DS.B.x
        DC_A[i] = np.mean(tsAs)
        DC_B[i] = np.mean(tsBs)
        spdxs[i] = ((Axs[i] - BalancedAx)*np.sin(DS.Aangle) + (Bxs[i] - BalancedBx)*np.sin(DS.Bangle))
        spoffs[i] = ((Axs[i] - BalancedAx)*np.sin(DS.Aangle) - (Bxs[i] - BalancedBx)*np.sin(DS.Bangle))/2

        DS.MoveGapBy(dxs[i+1])

        for j in range(avg):
            cache = cache + ms.Single_sided_PSD(tsAs[j,:], ps5000a.getfs(), window)
            cache_drive = cache_drive + ms.Single_sided_PSD(tsAs[j,:], ps5000a.getfs(), window)
        cache = cache/avg
        cache_drive = cache_drive/avg
        cache.tofile(wd+'\\'+'signal_{i:d}.bin'.format(i = i), sep = '')
        cache_drive.tofile(wd+'\\'+'drive_{i:d}.bin'.format(i = i), sep = '')
        print('Progress:{p:.2%}'.format(p = (i+1)/(len(sp_relative_xs))))

        powers = [0]*len(fcs)
        powers_drive = [0]*len(fcs)
        for j in range(len(fcs)):
            j1 = np.argmin(np.absolute(ps5000a.f-(fcs[j] - fspans[j]/2)))
            j2 = np.argmin(np.absolute(ps5000a.f-(fcs[j] + fspans[j]/2)))
            powers[j] = np.sum(cache[j1:j2])/(j2-j1+1)
            powers_drive[j] = np.sum(cache_drive[j1:j2])/(j2-j1+1)

        powers_on_step[i] = powers
        drives_on_step[i] = powers_drive
        print(10*np.log10(powers))
        print(10*np.log10(powers_drive))
        print(DC_A[i])

    powers_on_step = np.array(powers_on_step)
    drives_on_step = np.array(drives_on_step)
    for i in range(len(powers_on_step[0])):
        plt.plot(spoffs, powers_on_step[:,i],'*')
        plt.plot(spoffs, drives_on_step[:,i],'o')
    plt.plot(spoffs, DC_A, '+')
    plt.plot(spoffs, DC_B, 'D')
    plt.show()
    ps5000a.f.tofile(wd+'\\'+'freq.bin', sep = '')
    (np.array(Axs)).tofile(wd+'\\'+'Axs.bin', sep = '')
    (np.array(Bxs)).tofile(wd+'\\'+'Bxs.bin', sep = '')
    (np.array(DC_A)).tofile(wd+'\\'+'DC_A.bin', sep = '')
    (np.array(DC_B)).tofile(wd+'\\'+'DC_B.bin', sep = '')
    (np.array(spdxs)).tofile(wd+'\\'+'spdxs.bin', sep = '')
    (np.array(spoffs)).tofile(wd+'\\'+'spoffs.bin', sep = '')
    (np.array(powers_on_step)).tofile(wd+'\\'+'powers.bin', sep = '')
    (np.array(drives_on_step)).tofile(wd+'\\'+'drives.bin', sep = '')
    ps5000a.configurePSD(0.25, 1953125)

def chirp_check_cm(fcs, fspans, df, avg, ps5000a, FG, FGwrite, AWGF, filename = None, wd = None,  amps = 0.25, dc = 0.365, AWGfs = 3.2e6, AWGN = 12.8e6, offset = None, refOff = 0, phis = None):
    fcs = np.array(fcs)
    fspans = np.array(fspans)
    amps = np.array(amps)
    arbdata = np.zeros(int(AWGN))

    if FGwrite:
        if AWGF == 'chirp':
            arbdata = FG.chirps('test', fcs - fspans/2, fcs + fspans/2, AWGfs, AWGN, dc+amps, dc-amps, phis = phis)
        elif AWGF == 'sweep':
            FG.sweep(fcs[0]-fspans[0]/2, fcs[0]-fspan[0]/2, 1/df, amps[0], dc)
    factor = 0.25/df
    ps5000a.configurePSD(0.25/factor, 1953125) #gives sampling intervals of 256ns
    if offset is None:
        offset = -int(np.round(1e3/0.930*(1/df)/ps5000a.timeInternalns.value))
    ps5000a.getTimeSignal()
    temp = ps5000a.PSDfromTS(ps5000a.chs['A'].timeSignal[0:offset], ps5000a.getfs())
    window = signal.windows.tukey(len(ps5000a.chs['A'].timeSignal[0:offset]), 0)
    f = ps5000a.f[0:len(temp)]
    tsAs = np.zeros((avg, len(ps5000a.chs['A'].timeSignal[0:offset])))
    tsBs = np.zeros((avg, len(ps5000a.chs['A'].timeSignal[0:offset])))
    tsCs = np.zeros((avg, len(ps5000a.chs['A'].timeSignal[0:offset])))
    dftAs = np.zeros((avg, len(f)), dtype = np.complex)
    dftBs = np.zeros((avg, len(f)), dtype = np.complex)
    dftCs = np.zeros((avg, len(f)), dtype = np.complex)

    for i in range(avg):
        ps5000a.getTimeSignal(trigger = 'D', triggerThreshold = 8000)
        tsAs[i,:] = ps5000a.chs['A'].timeSignal[0:offset]
        tsBs[i,:] = ps5000a.chs['B'].timeSignal[0:offset]
        tsCs[i,:] = ps5000a.chs['C'].timeSignal[0:offset]
        print(i)

    psds_all = [[]]*avg
    for i in range(avg):
        powers, gamma, psds_all[i], psdCs, freqs, popts= Piecewise(fcs, fspans, tsAs[i,:], tsBs[i,:], f, ps5000a.getfs())
        psds_all[i][0] = np.roll(psds_all[i][0], int(refOff))
        if i == 0:
            psdavg = np.power(10, (psds_all[i][1]- psds_all[i][0])/10)
        else:
            psdavg = psdavg + np.power(10, (psds_all[i][1]- psds_all[i][0])/10)
    psdavg = 10*np.log10(psdavg/avg)
    '''
    plt.plot(psds_all[0][0], label='ref')
    plt.plot(psds_all[0][1], label='signal')
    plt.plot(psds_all[0][2], label='2f')
    plt.legend(loc="upper left")
    plt.show()
    plt.plot(freqs[1], psdavg)
    plt.show()
    print(powers)
    print(np.mean(tsCs))
    '''
    #return psds_all, arbdata
    return psds_all[0][0], psds_all[0][1]
def LineScan_sp_cm(fcs, fspans, ps5000a, DS, sp_relative_xs, wd, avg = 1):
    if not os.path.exists(wd):
        os.makedirs(wd)
    err = ''
    waitTime_after_moving = 10

    Nsteps = len(sp_relative_xs)

    temp = np.zeros(Nsteps+2)
    temp[1:-1] = sp_relative_xs
    dxs = np.zeros(Nsteps+1)
    for i in range(Nsteps):
        dxs[i] = temp[i+1] - temp[i]

    BalancedAx = DS.A.x
    BalancedBx = DS.B.x

    DS.MoveGapBy(dxs[0])

    ps5000a.configurePSD(1, 1953125)
    ps5000a.AC('A')
    ps5000a.ChangeRangeTo('A', 2000)
    ps5000a.DC('B')
    ps5000a.ChangeRangeTo('B', 2000)
    ps5000a.DC('C')
    ps5000a.ChangeRangeTo('C', 5000)

    ps5000a.configurePSD(0.25, 1953125)
    ps5000a.getTimeSignal()
    ps5000a.PSDfromTS(ps5000a.chs['A'].timeSignal, ps5000a.getfs())
    temp, _, _  = get_chirp(0.25, 1, ps5000a)
    temp = ms.Single_sided_fft(temp[0,:])
    f = ps5000a.f[0:len(temp)]

    range_A = 2000

    responses_on_step = [[]]*Nsteps
    DC_on_step = [[]]*Nsteps
    powers_on_step = [[]]*Nsteps

    Axs = [0]*Nsteps
    Bxs = [0]*Nsteps
    spdxs = [0]*Nsteps
    spoffs = [0]*Nsteps
    amps = np.array([0.05, 0.25])

    for i in range(len(sp_relative_xs)):
        #FG.chirps('test', fcs - fspans/2, fcs + fspans/2, 3.2e6, 12.8e6, dc+amps, dc-amps)
        #time.sleep(10)
        ps5000a.ChangeRangeTo('A', range_A)
        cache = np.zeros(len(temp))
        ps5000a.configurePSD(0.25, 1953125)
        [tsAs, tsBs, tsCs] = get_chirp(0.25, avg, ps5000a, trigger = 'D', triggerThreshold = 8000)
        window = signal.windows.tukey(len(tsAs[0,:]), 0)

        cache = np.zeros(len(temp))
        cache_drive = np.zeros(len(temp))

        Axs[i] = DS.A.x
        Bxs[i] = DS.B.x
        DC_on_step[i] = np.mean(tsCs)

        spdxs[i] = ((Axs[i] - BalancedAx)*np.sin(DS.Aangle) + (Bxs[i] - BalancedBx)*np.sin(DS.Bangle))
        spoffs[i] = ((Axs[i] - BalancedAx)*np.sin(DS.Aangle) - (Bxs[i] - BalancedBx)*np.sin(DS.Bangle))/2

        DS.MoveGapBy(dxs[i+1])

        for j in range(avg):
            _, _, psdAs, psdCs, freqs, _= Piecewise(fcs, fspans, tsAs[j,:], tsBs[j,:], f, ps5000a.getfs())
            if j == 0:
                responses_on_step[i] = np.power(10, (psdAs[1]- psdAs[0])/10)
                psdavg = np.power(10, np.array(psdAs)/10)
            else:
                responses_on_step[i] = responses_on_step[i] + np.power(10, (psdAs[1]- psdAs[0])/10)
                psdavg = psdavg + np.power(10, np.array(psdAs)/10)
        responses_on_step[i] = 10*np.log10(responses_on_step[i]/avg)
        psdavg = 10*np.log10(psdavg/avg)
        print('Progress:{p:.2%}'.format(p = (i+1)/(len(sp_relative_xs))))
        responses_on_step[i].tofile(wd+'\\'+'response_{i:d}.bin'.format(i = i), sep = '')

        powers = [0]*len(fcs)
        for j in range(len(fcs)):
            j1 = np.argmin(np.absolute(freqs[j]-(fcs[j] - fspans[j]/2)))
            j2 = np.argmin(np.absolute(freqs[j]-(fcs[j] + fspans[j]/2)))
            powers[j] = 10*np.log10(np.sum(np.power(10,psdavg[j,j1:j2]/10))/(j2-j1+1))

        powers_on_step[i] = powers
        print(powers_on_step[i])
        print(DC_on_step[i])

    (np.array(Axs)).tofile(wd+'\\'+'Axs.bin', sep = '')
    (np.array(Bxs)).tofile(wd+'\\'+'Bxs.bin', sep = '')
    (np.array(spdxs)).tofile(wd+'\\'+'spdxs.bin', sep = '')
    (np.array(spoffs)).tofile(wd+'\\'+'spoffs.bin', sep = '')
    (np.array(DC_on_step)).tofile(wd+'\\'+'DC_on_step.bin', sep = '')
    powers_on_step = np.array(powers_on_step)
    freqs[1].tofile(wd+'\\'+'freq.bin', sep = '')
    (np.array(powers_on_step))[:,0].tofile(wd+'\\'+'ref.bin', sep = '')
    (np.array(powers_on_step))[:,1].tofile(wd+'\\'+'powers.bin', sep = '')
    ps5000a.configurePSD(0.25, 1953125)

def tune_phi(ps5000a, FG, fc, fspan, ch = 'B'):
    samplef = 24
    ps5000a.configurePSD(samplef, 1953125)
    phis = np.linspace(0, 180, 181)
    re = np.linspace(0, 180 ,181)
    for i in range(181):
        FG.phi(phis[i])
        test = ps5000a.getPSD('B', offset = -1)
        re[i] = np.max(test[int((fc - fspan/2 - 10)/samplef):int((fc + fspan/2 + 10)/samplef)])
    FG.phi(phis[np.argmin(re)])
    ps5000a.configurePSD(0.25, 1953125)
    return phis, re

def fdrift(fcs, fspans, wd, monitorTime, ps5000a):
    if not os.path.exists(wd):
        os.makedirs(wd)
    I0 = ((np.array(fcs) - np.array(fspans)/2)*4).astype('i')
    I1 = ((np.array(fcs) + np.array(fspans)/2)*4).astype('i')
    data = [[]]*len(I0)
    times = [[]]
    t0 = time.time()
    t1 = time.time()
    while t1 - t0 < monitorTime:
        test = ps5000a.getPSD('A', avg = 1, offset = -18, trigger = 'D')
        t1 = time.time()
        print(t1-t0)
        times = times + [t1]
        for i in range(len(I0)):
            data[i] = data[i]+[(np.argmax(test[I0[i]:I1[i]])+I0[i])/4]
    del times[0]
    np.array(times).tofile(wd+'\\'+'times.bin', sep = '')
    for i in range(len(I0)):
        np.array(data[i]).tofile(wd+'\\'+'fm{i:d}.bin'.format(i = i), sep = '')
    for i in range(len(fcs)):
        plt.plot(times, (data[i] - np.mean(data[i]))/np.mean(data[i])*1e5)
    return times, data

def LineScan_QPD(fcs, fspans, UC8, ps5000a, UC8Amp, UC8_relative_xs, wd):
    if not os.path.exists(wd):
        os.makedirs(wd)

    axis_x = int(1)
    axis_y = int(2)
    err = ''
    UC8.CC_Set(1, err)
    time.sleep(1)
    UC8.SU_Set(axis_x, '+', UC8Amp, err)
    UC8.SU_Set(axis_x, '-', UC8Amp, err)
    UC8.SU_Set(axis_y, '+', UC8Amp, err)
    UC8.SU_Set(axis_y, '-', UC8Amp, err)
    UC8_controllerAddress = axis_y

    avg = 1
    waitTime_after_moving = 30

    temp = np.zeros(len(UC8_relative_xs)+2)
    temp[1:-1] = UC8_relative_xs
    steps = np.zeros(len(UC8_relative_xs)+1)
    for i in range(len(steps)):
        steps[i] = temp[i+1] - temp[i]

    UC8.PR(UC8_controllerAddress, steps[0], err)
    time.sleep(30)
    ps5000a.configurePSD(8, 1953125)
    UC8.CC_Set(2, err)
    time.sleep(1)
    balancer(ps5000a, UC8, axis_x, 'C', 1, 10)
    balancer(ps5000a, UC8, axis_y, 'A', -1, 50)
    UC8.CC_Set(1, err)
    time.sleep(1)

    ps5000a.configurePSD(1, 1953125)
    ps5000a.DC('A')
    ps5000a.ChangeRangeTo('A', 100)
    ps5000a.AC('B')
    ps5000a.ChangeRangeTo('B', 100)
    ps5000a.DC('C')
    ps5000a.ChangeRangeTo('C', 100)

    ps5000a.configurePSD(0.25, 1953125)
    ps5000a.getTimeSignal()
    ps5000a.PSDfromTS(ps5000a.chs['A'].timeSignal, ps5000a.getfs())
    temp, _, _  = get_chirp(0.25, avg, ps5000a)
    temp = ms.Single_sided_fft(temp[0,:])
    f = ps5000a.f[0:len(temp)]
    range_A = 100
    powers_on_step_A = [[]]*len(UC8_relative_xs)
    powers_on_step_C = [[]]*len(UC8_relative_xs)

    for i in range(len(UC8_relative_xs)):
        ps5000a.ChangeRangeTo('A', range_A)
        ps5000a.configurePSD(0.25, 1953125)
        [tsAs, tsBs, tsCs] = get_chirp(0.25, avg, ps5000a, trigger = 'D', triggerThreshold = 8000)
        window = signal.windows.tukey(len(tsAs[0,:]), 0.8)

        psdA = np.zeros(len(temp))
        psdB = np.zeros(len(temp))
        psdC = np.zeros(len(temp))

        t0 = time.time()
        UC8.CC_Set(1, err)
        time.sleep(1)
        if i != len(UC8_relative_xs)-1:
            UC8.PR(UC8_controllerAddress, steps[i+1], err)
        '''
        for j in range(avg):
            psdA = psdA + ms.Single_sided_PSD(tsAs[j,:], ps5000a.getfs(), window)
            psdB = psdB + ms.Single_sided_PSD(tsBs[j,:], ps5000a.getfs(), window)
            psdC = psdC + ms.Single_sided_PSD(tsCs[j,:], ps5000a.getfs(), window)
        psdA = 10*np.log10(20*np.square(psdA/avg))
        psdB = 10*np.log10(20*np.square(psdB/avg))
        psdC = 10*np.log10(20*np.square(psdC/avg))
        '''
        print('Progress:{p:.2%}'.format(p = (i+1)/(len(UC8_relative_xs))))

        if i == 0:
            powers_A, gamma, psdAs, psdBs, freqs, popts = Piecewise(fcs, fspans, tsAs[0], tsBs[0], f, ps5000a.getfs(), gamma = None)
            powers_C, gamma, psdCs, psdBs, freqs, popts = Piecewise(fcs, fspans, tsCs[0], tsBs[0], f, ps5000a.getfs(), gamma = None)
        else:
            powers_A, _, psdAs, psdBs, freqs, popts = Piecewise(fcs, fspans, tsAs[0], tsBs[0], f, ps5000a.getfs(), gamma)
            powers_C, _, psdCs, psdBs, freqs, popts = Piecewise(fcs, fspans, tsCs[0], tsBs[0], f, ps5000a.getfs(), gamma)
        powers_on_step_A[i] = powers_A
        powers_on_step_C[i] = powers_C
        print(powers_A)
        print(powers_C)

        t1 = time.time()
        time.sleep(max(waitTime_after_moving - (t1 - t0), 0))
        ps5000a.configurePSD(8, 1953125)
        UC8.CC_Set(2, err)
        time.sleep(1)
        balancer(ps5000a, UC8, axis_x, 'C', 1, 10)
        balancer(ps5000a, UC8, axis_y, 'A', -1, 50)

    powers_on_step_A = np.array(powers_on_step_A)
    powers_on_step_C = np.array(powers_on_step_C)
    for i in range(len(powers_on_step_A[0])):
        plt.plot(UC8_relative_xs, powers_on_step_A[:,i], label=str(fcs[i]))
        plt.plot(UC8_relative_xs, powers_on_step_C[:,i], label=str(fcs[i]))
    plt.legend(loc="upper left")
    plt.show()
    ps5000a.f.tofile(wd+'\\'+'freq.bin', sep = '')
    (np.array(UC8_relative_xs)).tofile(wd+'\\'+'UC2steps.bin', sep = '')
    (np.array(powers_on_step_A)).tofile(wd+'\\'+'powers_A.bin', sep = '')
    (np.array(powers_on_step_C)).tofile(wd+'\\'+'powers_C.bin', sep = '')
    ps5000a.configurePSD(0.25, 1953125)
    return powers_on_step_A, powers_on_step_C

def balancer(ps5000a, UC8, axis, ch, direction, amp, offset = 0, UC8_ch = 2, debug = True):
    err = ''
    imb = np.mean(ps5000a.getTimeSignal(ch, reportOverflow = debug)) - offset
    sigma = np.std(ps5000a.chs[ch].timeSignal)
    scale = 1
    while np.absolute(imb) > 0.5*sigma:
        UC8.PR(axis, np.sign(imb)*direction*np.max([1, int(amp*scale)]), err)
        WaitForReady(UC8, axis, UC8_ch)
        temp = np.mean(ps5000a.getTimeSignal(ch, reportOverflow = debug)) - offset
        scale = np.min([100, np.absolute(0.2*scale*temp/(imb-temp))])
        imb = temp
        sigma = np.std(ps5000a.chs[ch].timeSignal)
    return imb, ps5000a.chs['A'].overflow or ps5000a.chs['C'].overflow

def Piecewise_PSD(fcs, fspans, psdA, psdC, f, gamma = None):
    fcs = np.array(fcs)
    fspans = np.array(fspans)
    fstarts = fcs - fspans/2
    fends = fcs + fspans/2
    Nj = len(fcs)
    powers = np.zeros(Nj)

    check1s = np.zeros(Nj, dtype = np.int64)
    check2s = np.zeros(Nj, dtype = np.int64)
    for f1, f2, j in zip(fstarts, fends, list(range(Nj))):
        check1s[j] = np.argmin(np.absolute(f-f1+100))
        check2s[j] = np.argmin(np.absolute(f-f2-100))

    psdAs = [[]]*Nj
    psdCs = [[]]*Nj
    freqs = [[]]*Nj

    for j in range(Nj):
        psdAs[j] = psdA[check1s[j]:check2s[j]]
        psdCs[j] = psdC[check1s[j]:check2s[j]]
        freqs[j] = f[check1s[j]:check2s[j]]

    maxIs = np.zeros(Nj, dtype = np.int64)
    popts = [[]]*Nj
    gamma_provided = 0
    if gamma is None:
        gamma = [10]*Nj
        gamma_provided = False
    else:
        gamma_provided = True
    for j in range(Nj):
        maxIs[j] = np.argmax(psdAs[j])
        fitlo = min(300, maxIs[j])
        fithi = min(100,len(psdAs[j] - maxIs[j]))
        if gamma_provided:
            try:
                popts[j], pcov = curve_fit(fitFunc, freqs[j][maxIs[j]-fitlo:maxIs[j]+fithi], psdAs[j][maxIs[j]-fitlo:maxIs[j]+fithi] - psdCs[j][maxIs[j]-fitlo:maxIs[j]+fithi] , p0 =[gamma[j], freqs[j][maxIs[j]],  180], \
                    bounds = ([0.99*gamma[j], freqs[j][maxIs[j]] - 50,  50],[1.01*gamma[j], freqs[j][maxIs[j]]+50,  270]))
            except RuntimeError:
                print('fit failed. 0 given as result')
                popts[j] = [0,0,0]
        else:
            try:
                popts[j], pcov = curve_fit(fitFunc, freqs[j][maxIs[j]-fitlo:maxIs[j]+fithi], psdAs[j][maxIs[j]-fitlo:maxIs[j]+fithi] - psdCs[j][maxIs[j]-fitlo:maxIs[j]+fithi] , p0 =[gamma[j], freqs[j][maxIs[j]],  180], \
                    bounds = ([0, freqs[j][maxIs[j]] - 50,  50],[50, freqs[j][maxIs[j]]+50,  270]))
            except RuntimeError:
                print('fit failed. 0 given as result')
                popts[j] = [0,0,0]
        gamma[j] = popts[j][0]
        powers[j] = popts[j][2]
    return powers, gamma, freqs, popts
def z_sweep(axis_z, ps5000a, UC8, start, end, step, fc, offset_x = 0, offset_y = 0):
    HAs = []
    HBs = []
    HCs = []
    zs = []
    xs = []
    ys = []
    axis_x = int(1)
    axis_y = int(2)
    f = [fc - 200, fc + 200]
    ps5000a.getTimeSignal()
    ps5000a.PSDfromTS(ps5000a.chs['A'].timeSignal, ps5000a.getfs())
    idx = [0,0]
    idx[0] = np.argmin(np.absolute(ps5000a.f-f[0]))
    idx[1] = np.argmin(np.absolute(ps5000a.f-f[1]))
    avg = 20
    axis_z.Quasi_MoveTo(start)
    z = axis_z.Position()
    while np.sign(end-start)*z < np.sign(end-start)*end:
        axis_z.MoveBy(np.sign(end-start)*step)
        ps5000a.DC('A')
        ps5000a.DC('C')
        QPD_T_balance(ps5000a, UC8, offset_x, offset_y, 2)
        ps5000a.AC('A')
        ps5000a.AC('C')
        #T_search(216000, UC8, ps5000a, axis_x, axis_y, 1e10)

        x, y = GetOffsets(ps5000a)
        print('x = {:.4f}, y = {:.4f}'.format(x, y))
        xs = xs + [x]
        ys = ys + [y]

        HA, HB, HC = max_in_psd(avg, ps5000a, fc)
        HAs = HAs + [HA]
        HBs = HBs + [HB]
        HCs = HCs + [HC]
        z = axis_z.Position()
        print([z, HA, HB, HC])
        zs = zs + [z]
    return zs, HAs, HBs, HCs, xs, ys
def z_search(axis_z, ps5000a, UC8):
    HAs = []
    HBs = []
    HCs = []
    zs = []
    axis_x = int(1)
    axis_y = int(2)

    f = [215800, 216200]
    ps5000a.getTimeSignal()
    ps5000a.PSDfromTS(ps5000a.chs['A'].timeSignal, ps5000a.getfs())
    idx = [0,0]
    idx[0] = np.argmin(np.absolute(ps5000a.f-f[0]))
    idx[1] = np.argmin(np.absolute(ps5000a.f-f[1]))

    z = axis_z.Position()
    avg = 20
    search = True
    while search:
        QPD_T_balance(ps5000a, UC8)
        for i in range(avg):
            ps5000a.getTimeSignal()
            if i == 0:
                psdA = ms.Single_sided_PSD(ps5000a.chs['A'].timeSignal, ps5000a.getfs())
                psdB = ms.Single_sided_PSD(ps5000a.chs['B'].timeSignal, ps5000a.getfs())
                psdC = ms.Single_sided_PSD(ps5000a.chs['C'].timeSignal, ps5000a.getfs())
            else:
                psdA = psdA + ms.Single_sided_PSD(ps5000a.chs['A'].timeSignal, ps5000a.getfs())
                psdB = psdB + ms.Single_sided_PSD(ps5000a.chs['B'].timeSignal, ps5000a.getfs())
                psdC = psdC + ms.Single_sided_PSD(ps5000a.chs['C'].timeSignal, ps5000a.getfs())
        psdA = psdA / avg
        HA = np.max(psdA[idx[0]:idx[1]])
        HAs = HAs + [HA]
        psdB = psdB / avg
        HB = np.max(psdB[idx[0]:idx[1]])
        HBs = HBs + [HB]
        psdC = psdC / avg
        HC = np.max(psdC[idx[0]:idx[1]])
        HCs = HCs + [HC]
        z = axis_z.Position()
        zs = zs + [z]

        if len(zs) == 1:
            step = -0.5
            search = True
        else:
            step = - 6e6 * (HCs[-1] - HCs[-2])/(zs[-1]-zs[-2])
            step = np.sign(step)*np.min([2*np.absolute(zs[-1] - zs[-2]), np.absolute(step)])
            step = np.sign(step)*np.min([0.5, np.absolute(step)])
            step = np.sign(step)*np.max([0.01, np.absolute(step)])
            if len(zs) > 10:
                flag = True
                for index in range(1):
                    flag = flag and ApproxMin(HCs[-5:], index, np.min(HCs))
                if flag:
                    search = False
        if (search):
            axis_z.MoveBy(round(step,3))
        print([z, HA, HC, step])
    print(zs[-1], HCs[-1])
    return zs, HAs, HBs, HCs

def ApproxMin(samples, index, target_min):
    sample = samples[index]
    del samples[index]
    if sample - target_min < 2*np.std(samples):
        return True
    else:
        return False

def QPD_T_balance(ps5000a, UC8, offset_x = 0, offset_y = 0, UC8_ch = 2, debug = True):
    axis_x = int(1)
    axis_y = int(2)
    OF = True
    balanced = False
    while not balanced:
        temp, OF = balancer(ps5000a, UC8, axis_x, 'C', 1, 4, offset_x, UC8_ch, debug = debug)
        #print(np.mean(ps5000a.chs['A'].timeSignal), np.mean(ps5000a.chs['C'].timeSignal))
        if np.absolute(np.mean(ps5000a.chs['A'].timeSignal) - offset_y) < np.std(ps5000a.chs['A'].timeSignal):
            balanced = True
        if balanced:
            break
        temp, OF = balancer(ps5000a, UC8, axis_y, 'A', -1, 4, offset_y, UC8_ch, debug = debug)
        #print(np.mean(ps5000a.chs['A'].timeSignal), np.mean(ps5000a.chs['C'].timeSignal))
        if np.absolute(np.mean(ps5000a.chs['C'].timeSignal) - offset_x) < np.std(ps5000a.chs['C'].timeSignal):
            balanced = True

def WaitForReady(UC, axis, ch):
    status = 99
    err = ''
    count = 0
    while status is not 0:
        _, status, _ = UC.TS(axis, status, err)
        time.sleep(0.1)
        count = count+1
        if status == -2147483648:
            UC.RS(err)
            UC.MR(err)
            UC.CC_Set(ch, err)
            print(status)
            break

def max_in_psd(avg, ps5000a, fc):
    f = [fc-200, fc+200]
    ps5000a.getTimeSignal()
    ps5000a.PSDfromTS(ps5000a.chs['A'].timeSignal, ps5000a.getfs())
    idx = [0,0]
    idx[0] = np.argmin(np.absolute(ps5000a.f-f[0]))
    idx[1] = np.argmin(np.absolute(ps5000a.f-f[1]))
    for i in range(avg):
        ps5000a.getTimeSignal()
        if i == 0:
            psdA = ms.Single_sided_PSD(ps5000a.chs['A'].timeSignal, ps5000a.getfs())
            psdB = ms.Single_sided_PSD(ps5000a.chs['B'].timeSignal, ps5000a.getfs())
            psdC = ms.Single_sided_PSD(ps5000a.chs['C'].timeSignal, ps5000a.getfs())
        else:
            psdA = psdA + ms.Single_sided_PSD(ps5000a.chs['A'].timeSignal, ps5000a.getfs())
            psdB = psdB + ms.Single_sided_PSD(ps5000a.chs['B'].timeSignal, ps5000a.getfs())
            psdC = psdC + ms.Single_sided_PSD(ps5000a.chs['C'].timeSignal, ps5000a.getfs())
    psdA = psdA / avg
    HA = np.max(psdA[idx[0]:idx[1]])
    psdB = psdB / avg
    HB = np.max(psdB[idx[0]:idx[1]])
    psdC = psdC / avg
    HC = np.max(psdC[idx[0]:idx[1]])
    return HA, HB, HC

def maxs_in_psd(avg, ps5000a, fcs):
    fstarts = np.array(fcs) - 200
    fends = np.array(fcs) + 200
    ps5000a.getTimeSignal()
    ps5000a.PSDfromTS(ps5000a.chs['A'].timeSignal, ps5000a.getfs())
    HAs = np.zeros(len(fcs))
    HBs = np.zeros(len(fcs))
    HCs = np.zeros(len(fcs))
    idx0s = np.zeros(len(fcs), dtype = int)
    idx1s = np.zeros(len(fcs), dtype = int)
    for i in range(len(fcs)):
        idx0s[i] = np.argmin(np.absolute(ps5000a.f-fstarts[i]))
        idx1s[i] = np.argmin(np.absolute(ps5000a.f-fends[i]))
    for i in range(avg):
        ps5000a.getTimeSignal()
        if i == 0:
            psdA = ms.Single_sided_PSD(ps5000a.chs['A'].timeSignal, ps5000a.getfs())
            psdB = ms.Single_sided_PSD(ps5000a.chs['B'].timeSignal, ps5000a.getfs())
            psdC = ms.Single_sided_PSD(ps5000a.chs['C'].timeSignal, ps5000a.getfs())
        else:
            psdA = psdA + ms.Single_sided_PSD(ps5000a.chs['A'].timeSignal, ps5000a.getfs())
            psdB = psdB + ms.Single_sided_PSD(ps5000a.chs['B'].timeSignal, ps5000a.getfs())
            psdC = psdC + ms.Single_sided_PSD(ps5000a.chs['C'].timeSignal, ps5000a.getfs())
    psdA = psdA / avg
    psdB = psdB / avg
    psdC = psdC / avg
    for i in range(len(fcs)):
        HAs[i] = np.max(psdA[idx0s[i]:idx1s[i]])
        HBs[i] = np.max(psdB[idx0s[i]:idx1s[i]])
        HCs[i] = np.max(psdC[idx0s[i]:idx1s[i]])
    return HAs, HBs, HCs

def T_search(fc, UC8, ps5000a, axis_x, axis_y, alpha):
    p_last = T_search_cost(ps5000a, fc)
    dx = 20
    dy = 20
    err = ''
    powers = [p_last]
    steps = [20]

    converged = False
    while not converged:
        UC8.PR(axis_x, dx, err)
        dx_last = dx
        pdx = T_search_cost(ps5000a, fc)
        dpx = pdx - p_last
        dx = -alpha*(dpx/dx_last)
        dx = np.sign(dx)*InRange(1,int(np.absolute(2*dx_last)),np.absolute(dx))
        dx = np.sign(dx)*InRange(1,50,np.absolute(dx))
        dx = int(dx)
        if len(steps) > 4:
            if dx*dx_last < 0 and dx_last * steps[-2] < 0:
                dx = int(-dx *0.5)

        UC8.PR(axis_y, dy, err)
        dy_last = dy
        pdy = T_search_cost(ps5000a, fc)
        dpy = pdy - pdx
        dy = -alpha*(dpy/dy_last)
        dy = np.sign(dy)*InRange(1,int(np.absolute(2*dy_last)),np.absolute(dy))
        dy = np.sign(dy)*InRange(1,50,np.absolute(dy))
        dy = int(dy)
        if len(steps) > 4:
            if dy*dy_last < 0 and dy_last * steps[-1] < 0:
                dy = int(-dy*0.5)

        dp = pdy - p_last
        HA, HB, HC = max_in_psd(5, ps5000a, fc)
        #print(np.absolute(dp/p_last))
        #print('dpx = {:.4e}, dx = {:d}, dpy = {:.4e}, dy = {:d}'.format(pdx, dx_last, pdy, dy_last))
        print('p = {:.4e}, HA = {:.4e}, HB = {:.4e}, HC = {:.4e}'.format(pdy, HA, HB, HC))
        p_last = pdy
        powers = powers + [pdx, pdy]
        steps = steps +[dx_last, dy_last]

        if len(powers) > 7:
            powers = np.array(powers)
            steps = np.array(steps)
            r = np.divide(powers[-6:] - powers[-7:-1], np.multiply(powers[-7:-1], steps[-6:]))
            #print(r)
            if np.max(np.absolute(r)) < 3e-4 and (powers[-1] - np.min(powers))/min(powers)<0.05:
                converged = True
            if max(np.absolute(steps[-6:])) < 5:
                converged = True
            powers = powers.tolist()
            steps = steps.tolist()


def T_search_cost(ps5000a, fc):
    avg = 5
    HA, HB, HC = max_in_psd(avg, ps5000a, fc)
    return HA + HC

def InRange(lo, hi, x):
    x = min(hi, x)
    x = max(lo, x)
    return x

def GetOffsets(ps5000a):
    ps5000a.DC('A')
    ps5000a.DC('C')
    RangeA = ps5000a.chs['A'].getRange()
    RangeC = ps5000a.chs['C'].getRange()
    ps5000a.AutoRange('A')
    ps5000a.AutoRange('C')
    ps5000a.getTimeSignal()
    x = np.mean(ps5000a.chs['A'].timeSignal)
    y = np.mean(ps5000a.chs['C'].timeSignal)
    ps5000a.AC('A')
    ps5000a.AC('C')
    ps5000a.ChangeRangeTo('A', RangeA)
    ps5000a.ChangeRangeTo('C', RangeC)
    return x, y

def max_in_fft(avg, ps5000a, fc, refch):
    f = [fc-200, fc+200]
    ps5000a.getTimeSignal()
    ps5000a.PSDfromTS(ps5000a.chs['A'].timeSignal, ps5000a.getfs())
    idx = [0,0]
    idx[0] = np.argmin(np.absolute(ps5000a.f-f[0]))
    idx[1] = np.argmin(np.absolute(ps5000a.f-f[1]))
    for i in range(avg):
        ps5000a.getTimeSignal(trigger = 'D', triggerThreshold = 8000)
        fftRef = ms.Single_sided_fft(ps5000a.chs[refch].timeSignal)
        fftRef = fftRef[idx[0]:idx[1]]
        IRef = np.argmax(np.absolute(fftRef))
        angle = np.mean(np.angle(fftRef[IRef-20:IRef-5]))
        refNum = np.complex(np.cos(angle), np.sin(angle))
        if i == 0:
            fftA = ms.Single_sided_fft(ps5000a.chs['A'].timeSignal)/refNum
            fftB = ms.Single_sided_fft(ps5000a.chs['B'].timeSignal)/refNum
            fftC = ms.Single_sided_fft(ps5000a.chs['C'].timeSignal)/refNum
        else:
            fftA = fftA + ms.Single_sided_fft(ps5000a.chs['A'].timeSignal)/refNum
            fftB = fftB + ms.Single_sided_fft(ps5000a.chs['B'].timeSignal)/refNum
            fftC = fftC + ms.Single_sided_fft(ps5000a.chs['C'].timeSignal)/refNum
    N = len(ps5000a.chs['A'].timeSignal)
    fftA = fftA[idx[0]:idx[1]] / avg
    IA = np.argmax(np.absolute(fftA))
    HA = fftA[IA]/(np.sqrt(ps5000a.getfs()*N))
    fftB = fftB[idx[0]:idx[1]] / avg
    IB = np.argmax(np.absolute(fftB))
    HB = fftB[IB]/(np.sqrt(ps5000a.getfs()*N))
    fftC = fftC[idx[0]:idx[1]] / avg
    IC = np.argmax(np.absolute(fftC))
    HC = fftC[IC]/(np.sqrt(ps5000a.getfs()*N))
    angle1 = np.mean(fftC[IC-15:IC-5])
    angle2 = np.mean(fftC[IC+5:IC+15])
    return HA, HB, HC, fftB, fftC, angle1, angle2

def max_in_fft_debug(avg, ps5000a, fc, refch):
    f = [fc-200, fc+200]
    ps5000a.getTimeSignal()
    ps5000a.PSDfromTS(ps5000a.chs['A'].timeSignal, ps5000a.getfs())
    idx = [0,0]
    idx[0] = np.argmin(np.absolute(ps5000a.f-f[0]))
    idx[1] = np.argmin(np.absolute(ps5000a.f-f[1]))

    fftRef = ms.Single_sided_fft(ps5000a.chs['B'].timeSignal)
    maxI = np.argmax(np.absolute(fftRef[idx[0]:idx[1]]))
    N = len(ps5000a.chs['A'].timeSignal)

    complexA = []
    complexB = []
    complexC = []

    for i in range(avg):
        ps5000a.getTimeSignal(trigger = 'D', triggerThreshold = 8000)
        fftRef = ms.Single_sided_fft(ps5000a.chs[refch].timeSignal)
        fftRef = fftRef[idx[0]:idx[1]]
        angleRef = np.angle(fftRef[maxI])
        normRef = np.complex(np.cos(angleRef), np.sin(angleRef))
        fftA = ms.Single_sided_fft(ps5000a.chs['A'].timeSignal)
        fftB = ms.Single_sided_fft(ps5000a.chs['B'].timeSignal)
        fftC = ms.Single_sided_fft(ps5000a.chs['C'].timeSignal)
        fftA = fftA[idx[0]:idx[1]]
        fftB = fftB[idx[0]:idx[1]]
        fftC = fftC[idx[0]:idx[1]]
        complexA = complexA + [fftA[maxI]/(normRef*(np.sqrt(ps5000a.getfs()*N)))]
        complexB = complexB + [fftB[maxI]/(normRef*(np.sqrt(ps5000a.getfs()*N)))]
        complexC = complexC + [fftC[maxI]/(normRef*(np.sqrt(ps5000a.getfs()*N)))]

    return np.mean(complexA), np.mean(complexB), np.mean(complexC), fftC

def cost(x, *arg):
    A, B, z0, phi0, phi1 = x
    z, amp, angle = arg
    amp_n = (amp - (np.max(amp) + np.min(amp))/2)/((np.max(amp) - np.min(amp))/2)
    angle_n = (angle - (np.max(angle) + np.min(angle))/2)/((np.max(angle) - np.min(angle))/2)
    fit_amp, fit_angle = fit_amp_angle(z, A, B, z0, phi0, phi1)
    fit_amp_n = (fit_amp - (np.max(amp) + np.min(amp))/2)/((np.max(amp) - np.min(amp))/2)
    fit_angle_n = (fit_angle - (np.max(angle) + np.min(angle))/2)/((np.max(angle) - np.min(angle))/2)
    return (np.sum(np.square(fit_amp_n - amp_n)) + 10*np.sum(np.square(fit_angle_n - angle_n)))/len(z)

def fit_amp_angle(z, A, B, z0, phi0, phi1):
    amp = np.sqrt(np.square(-A*(z-z0)+B*np.cos(phi1))+np.square(B*np.sin(phi1)))
    phi = phi0 + np.arccos(np.divide(-A*(z-z0)+B*np.cos(phi1),amp))
    return amp, phi

def cost_complex_test(x, *arg):
    A, B, z0, phi0, phi1 = x
    z, complexC = arg
    data_complex = complexC
    norm_factor = np.absolute(np.mean(data_complex))
    fit_complex = fit_complex_test(z, A, B, z0, phi0, phi1)/norm_factor
    data_complex = data_complex/norm_factor
    return  np.std(fit_complex - data_complex)

def fit_complex_test(z, A, B, z0, phi0, phi1):
    return A*(z-z0)*np.exp(1j*phi0) + B*np.exp(1j*(phi0 + phi1))

def balancer_prism(ps5000a, axis, ch, direction, amp, offset = 0, debug = True, rel = 0.5):
    imb = np.mean(ps5000a.getTimeSignal(ch, reportOverflow = debug)) - offset
    sigma = np.std(ps5000a.chs[ch].timeSignal)
    scale = 1
    while np.absolute(imb) > rel*sigma:
        axis.MoveBy(np.sign(imb)*direction*amp*scale)
        temp = np.mean(ps5000a.getTimeSignal(ch, reportOverflow = debug)) - offset
        scale = np.min([scale*2, np.absolute(0.2*scale*temp/(imb-temp)), 10])
        scale = np.max([0.001, scale])
        if imb == temp:
            scale = 1
        imb = temp
        sigma = np.std(ps5000a.chs[ch].timeSignal)
        if debug:
            print(imb, axis.Position(), scale)
    return imb, ps5000a.chs['A'].overflow or ps5000a.chs['C'].overflow

def process_signal_fft(zs, fftBs, fftCs, idxs):
    cpy_fftCs = [[]]*len(fftCs)
    for i in range(len(fftCs)):
        cpy_fftCs[i] = fftCs[i].copy()

    cpy_fftBs = [[]]*len(fftBs)
    for i in range(len(fftBs)):
        cpy_fftBs[i] = fftBs[i].copy()

    index = [int(len(fftCs[0])/2 - 500), int(len(fftCs[0])/2 + 500)]
    complexC = np.zeros(len(cpy_fftBs), dtype = np.complex)
    for i in range(len(cpy_fftCs)):
        bavg = np.mean(cpy_fftBs[i][index[0]:index[1]])
        complexC[i] = bavg
        cpy_fftCs[i][index[0]:index[1]] = np.divide(cpy_fftCs[i][index[0]:index[1]], cpy_fftBs[i][index[0]:index[1]])*bavg

    for i in range(len(cpy_fftCs)):
        cpy_fftCs[i] = cpy_fftCs[i]/complexC[0]*np.absolute(complexC[0])

    complexC = np.zeros(len(cpy_fftCs), dtype = np.complex)
    for i in range(len(cpy_fftCs)):
        maxI = np.argmax(np.absolute(cpy_fftCs[i]))
        indexarra = np.array(list(range(index[0], maxI - 100)) + list(range(maxI + 100, index[1])))
        temp = cpy_fftCs[i]
        complexC[i] = np.mean(temp[indexarra])

    p = np.polyfit(np.real(complexC), np.imag(complexC), 1)
    b = p[1]
    k = p[0]
    phi0 = np.arctan(k)
    A = ((np.real(complexC[0])+k*np.imag(complexC[0]))/np.sqrt(1+k**2) + \
        (np.real(complexC[-1])+k*np.imag(complexC[-1]))/np.sqrt(1+k**2))/(zs[0] - zs[-1])
    outphase_noise = b*np.cos(phi0)

    temp = complexC.copy()
    complexC = complexC - 1j*outphase_noise*np.exp(1j*phi0)
    p = np.polyfit(zs, np.real(complexC), 1)
    z0 = -p[1]/p[0]
    res = {'b':b, 'k':k, 'phi0':phi0, 'A':A, 'z0':z0}

    for i in range(len(cpy_fftCs)):
        cpy_fftCs[i][index[0]:index[1]] = cpy_fftCs[i][index[0]:index[1]] - 1j*outphase_noise*np.exp(1j*phi0)

    dict_re = {'zs_measured':zs, 'fftBs':fftBs, 'fftCs':fftCs, 'processed_C':cpy_fftCs, 'fit':res, 'complexC':temp, 'processed_complexC':complexC, 'idxs':idxs}
    return dict_re

def lock_in_amplifier(signal, ref, fs, freq):
    lowcut = freq - 1000
    highcut = freq + 1000
    T = len(signal)/fs
    ref1 = ms.butter_bandpass_filter(ref, lowcut, highcut, fs, order=3)
    ref2 = ms.HilbertTransform(ref1)
    X = np.dot(signal, ref1)/T
    Y = np.dot(signal, ref2)/T
    amp = np.absolute(np.complex(X, Y))
    angle = np.angle(np.complex(X, Y))
    return amp, angle

def lock_in_amplifier_t(signal, ref, fs, freq):
    lowcut = freq - 1000
    highcut = freq + 1000
    T = len(signal)/fs
    ref1 = ms.butter_bandpass_filter(ref, lowcut, highcut, fs, order=3)
    ref2 = ms.HilbertTransform(ref1)
    X = ms.butter_lowpass_filter(np.multiply(signal, ref1)/T, 1e3, fs, order=3)
    Y = ms.butter_lowpass_filter(np.multiply(signal, ref2)/T, 1e3, fs, order=3)
    S = X+1j*Y
    amps = np.absolute(S)
    angles = np.angle(S)
    return amps, angles

def lock_in(freq, FG, ps5000a, num_of_pieces = 1):
    ch = 2
    freq_in_FG = float(FG.inst.query('SOUR{:d}:FREQ?'.format(ch)))
    if (np.absolute(freq_in_FG - freq) > 1e-4 ):
        FG.Sine(freq, 5e-3, ch = ch)
        time.sleep(1)
    ps5000a.getTimeSignal()
    n = int(len(ps5000a.chs['C'].timeSignal)/num_of_pieces)
    amps = []
    angles = []
    for i in range(num_of_pieces):
        signal = ps5000a.chs['C'].timeSignal[i*n:(i+1)*n]
        amp, angle = lock_in_amplifier(ps5000a.chs['C'].timeSignal[i*n:(i+1)*n], \
            ps5000a.chs['B'].timeSignal[i*n:(i+1)*n], ps5000a.getfs(), freq)
        angles = angles + [angle]
        amps = amps + [amp]
    return freq, amps, angles

def lock_in_check(freq, ps5000a, num_of_pieces = 1):
    ps5000a.getTimeSignal()
    n = int(len(ps5000a.chs['A'].timeSignal)/num_of_pieces)
    amps = []
    angles = []
    for i in range(num_of_pieces):
        signal = ps5000a.chs['A'].timeSignal[i*n:(i+1)*n]
        amp, angle = lock_in_amplifier(ps5000a.chs['A'].timeSignal[i*n:(i+1)*n], \
            ps5000a.chs['B'].timeSignal[i*n:(i+1)*n], ps5000a.getfs(), freq)
        angles = angles + [angle]
        amps = amps + [amp]
    return freq, amps, angles

def Quasi_MoveBy(step, upperlimit, ps5000a, axis_z, prism_x):
    initial_z = axis_z.Position()
    target_z = initial_z + step
    while np.absolute(axis_z.Position() - target_z) > upperlimit:
        axis_z.MoveBy(np.sign(step)*upperlimit)
        axis_z.WaitForReady()
        ps5000a.configurePSD(1, int(1953125/5))
        ps5000a.ChangeRangeTo('A', 1000)
        balancer_prism(ps5000a, prism_x, 'A', -1, 0.005, offset = 0, debug = False, rel = 1)
        ps5000a.ChangeRangeTo('A', 2000)
    axis_z.MoveBy(-(axis_z.Position() - target_z))

def sweep_back_forward(prism_x, axis_z, ps5000a, zs0, n_sweeps, wd, minzs, fcs, windowf = None):
    t0 = time.time()
    dict_to_save = {}
    for i in range(n_sweeps):
        fordata = [0]*len(zs0)
        zs = zs0
        if np.any(minzs):
            for j in range(len(minzs)):
                minz = minzs[j]
                sweep_range = ((np.random.rand(1)-0.5)*2*3 + 3.2)[0]
                target_zs = np.around(np.linspace(minz - sweep_range, minz + sweep_range, 16), decimals = 2)
                target_zs = target_zs.tolist()
                zs = zs + target_zs
                fordata = fordata + [j+1]*len(target_zs)

        zs = np.array(zs)
        fordata = np.array(fordata)
        inds = zs.argsort()
        zs = zs[inds]
        fordata = fordata[inds]

        zs_index = np.argwhere(np.logical_and(zs>0, zs<25))
        zs = zs[zs_index.flatten()]
        fordata = fordata[zs_index.flatten()]
        zs = zs.tolist()
        fordata = fordata.tolist()

        z_last = axis_z.Position()
        if np.absolute(z_last - zs[0]) > np.absolute(z_last - zs[-1]):
            zs.reverse()
            fordata.reverse()

        Nfs = int(len(fcs)/2)
        zs_measured = []
        refs_allfs = [[]]*Nfs
        signals_allfs = [[]]*Nfs
        idxs = []
        maxfs = []
        psdCs = []

        for z, j in zip(zs, range(len(zs))):
            axis_z.Quasi_MoveTo(z)
            ps5000a.configurePSD(1, int(1953125/5))
            ps5000a.ChangeRangeTo('A', 1000)
            balancer_prism(ps5000a, prism_x, 'A', -1, 0.005, offset = 0, debug = False, rel = 2)
            ps5000a.ChangeRangeTo('A', 2000)

            if fordata[j]:
                signal, ref, idx, psdC = getResult_fft_sep(5, ps5000a, fcs, prism_x, windowf = windowf)
            else:
                signal, ref, idx, psdC = getResult_fft_sep(10, ps5000a, fcs, prism_x, windowf = windowf)

            for ifs in range(Nfs):
                 refs_allfs[ifs] = refs_allfs[ifs] + [signal[2*ifs+1]]
                 signals_allfs[ifs] = signals_allfs[ifs] + [signal[2*ifs]]

            tempf = ps5000a.f[idx[2]:idx[3]]
            maxf = tempf[np.argmax(np.absolute(signal[1]))]

            maxfs = maxfs + [maxf]
            idxs = idxs + [idx.tolist()]
            zs_measured = zs_measured +[axis_z.Position()]
            psdCs = psdCs + [psdC.tolist()]

            print(zs_measured[-1])
            z_last = axis_z.Position()

        for ifs in range(Nfs):
            dict_to_save = process_signal_fft_ref(np.array(zs_measured), refs_allfs[ifs], signals_allfs[ifs], idxs)
            j = len(os.listdir(wd))
            dict_to_save['maxfs'] = maxfs
            dict_to_save['fordata'] = fordata
            dict_to_save['psdCs'] = psdCs
            dict_to_json(dict_to_save, wd, 'result_{:d}_ifs_{:d}.json'.format(j, ifs))

        datas_list = BKA_get_data(wd)
        N = len(datas_list)
        zs_data_list = [[]]*N
        relative_z_data_list = [[]]*N
        for i in range(N):
            relative_z_data_list[i], _, _, _, _, zs_data_list[i], _, _ = BKAdata(datas_list[i], [0])
            ind = np.argmin(relative_z_data_list[i])
            zs_data_index = np.argwhere(np.absolute(zs_data_list[i] - zs_data_list[ind]<5))
            zs_data = zs_data_list[i][zs_data_index.flatten()]
            relative_z_data = relative_z_data_list[i][zs_data_index.flatten()]
            p = np.polyfit(relative_z_data, zs_data, 1)
            minzs[i] = p[1]
        print(time.time() - t0)

def PID(setPoint, kp, ki, kd, fc, ps5000a, monitorTime = None):
    et = []
    t = []
    DCs = []
    ps = []
    t0 = time.time()
    InRange = True
    i = 0
    sum = 0
    Ti = 0.125
    lo = 0
    hi = 2
    initial = 1
    plt.axis([-np.pi, np.pi, 0, 12])
    while InRange < 10:
        freq, amps, angles = lock_in_check(fc, ps5000a, num_of_pieces = 1)
        t = t + [time.time() - t0]
        et = et + [ms.angle_diff(setPoint, np.mean(angles))]
        ps = ps + [np.mean(amps)]
        sum = sum + et[-1]
        if i == 0:
            freq, amps, angles = lock_in_check(fc, ps5000a, num_of_pieces = 1)
            t = t + [time.time() - t0]
            et = et + [ms.angle_diff(setPoint, np.mean(angles))]
            sum = sum + et[-1]
            ps = ps + [np.mean(amps)]
        DC = initial - kp*et[-1] - ki*Ti*sum - kd*(et[-1] - et[-2])/(t[-1] - t[-2])
        if DC > lo and DC < hi:
            ps5000a.SigGenDC(DC)
            DCs = DCs + [DC]
            InRange = 0
        else:
            InRange = InRange + 1
            del t[-1]
            del et[-1]
            del ps[-1]
        if monitorTime is not None:
            if t[-1] > monitorTime:
                InRange = 10
        i = i + 1
        print(et[-1], ps[-1], DC)
        plt.plot(et[-1], ps[-1], '.')
        plt.pause(0.001)
    plt.plot(t, et)
    plt.show()
    return t, et, DCs, ps

def PID_ver2(setPoint, kp, ki, kd, fc, ps5000a, wd, initial = 1, monitorTime = None):
    bufferN = 10000
    et = np.zeros(bufferN)
    t = np.zeros(bufferN)
    DCs = np.zeros(bufferN)
    ps = np.zeros(bufferN)
    t0 = time.time()
    InRange = True
    i = 0
    sum = 0
    Ti = 0.125
    filei = 0
    lo = 0
    hi = 2
    initial = initial
    while InRange < 10:
        freq, amps, angles = lock_in_check(fc, ps5000a, num_of_pieces = 1)
        t[i] = time.time() - t0
        et[i] = ms.angle_diff(setPoint, np.mean(angles))
        ps[i] = np.mean(amps)
        sum = sum + et[i]
        if i == 0:
            i = i + 1
            freq, amps, angles = lock_in_check(fc, ps5000a, num_of_pieces = 1)
            t[i] = time.time() - t0
            et[i] = ms.angle_diff(setPoint, np.mean(angles))
            ps[i] = np.mean(amps)
            sum = sum + et[i]
        DC = initial - kp*et[i] - ki*Ti*sum - kd*(et[i] - et[i-1])/(t[i] - t[i-1])
        if DC > lo and DC < hi:
            ps5000a.SigGenDC(DC)
            DCs[i] = DC
            InRange = 0
        else:
            InRange = InRange + 1
            i = i - 1
        if monitorTime is not None:
            if t[i] > monitorTime:
                InRange = 10
        print(i, et[i], ps[i], DC)
        i = i + 1
        newi = i%bufferN
        #write files
        if newi != i:
            PIDtofile(t, et, DCs, ps, filei, wd)
            i = newi
            filei = filei + 1
    filei = filei + 1
    PIDtofile(t[:i], et[:i], DCs[:i], ps[:i], filei, wd)

def dict_to_json(target_dict, wd, filename):
    if not os.path.exists(wd):
        os.makedirs(wd)
    fp = open(wd+'\\'+filename, 'w')
    temp = {}
    temp['zs'] = target_dict['zs_measured'].tolist()
    temp['idxs'] = target_dict['idxs']
    temp['fordata'] = target_dict['fordata']
    temp['maxfs'] = target_dict['maxfs']
    temp['complexC_Re'] = (np.real(target_dict['complexC'])).tolist()
    temp['complexC_Im'] = (np.imag(target_dict['complexC'])).tolist()
    temp['processed_complexC_Re'] = (np.real(target_dict['processed_complexC'])).tolist()
    temp['processed_complexC_Im'] = (np.imag(target_dict['processed_complexC'])).tolist()
    temp['fit'] = target_dict['fit']
    temp['psdCs'] = target_dict['psdCs']
    fftBs_Re = [[]]*len(target_dict['fftBs'])
    fftBs_Im = [[]]*len(target_dict['fftBs'])
    fftCs_Re = [[]]*len(target_dict['fftBs'])
    fftCs_Im = [[]]*len(target_dict['fftBs'])
    fftpCs_Re = [[]]*len(target_dict['fftBs'])
    fftpCs_Im = [[]]*len(target_dict['fftBs'])
    for j in range(len(target_dict['fftBs'])):
        fftBs_Re[j] = (np.real(target_dict['fftBs'][j])).tolist()
        fftBs_Im[j] = (np.imag(target_dict['fftBs'][j])).tolist()
        fftCs_Re[j] = (np.real(target_dict['fftCs'][j])).tolist()
        fftCs_Im[j] = (np.imag(target_dict['fftCs'][j])).tolist()
        fftpCs_Re[j] = (np.real(target_dict['processed_C'][j])).tolist()
        fftpCs_Im[j] = (np.imag(target_dict['processed_C'][j])).tolist()
    temp['fftBs_Re'] = fftBs_Re
    temp['fftBs_Im'] = fftBs_Im
    temp['fftCs_Re'] = fftCs_Re
    temp['fftCs_Im'] = fftCs_Im
    temp['fftpCs_Re'] = fftpCs_Re
    temp['fftpCs_Im'] = fftpCs_Im
    json.dump(temp, fp)
    fp.close()

def getResult_fft(avg, ps5000a, fcs, prism_x):
    #only normalize to channel B's phase, ie, fftB is pure real
    nT = int(39063230/5)
    nT = 7812489
    ps5000a.configurePSD(0.099/avg, int(1953125/5))
    ps5000a.getTimeSignal(trigger = 'D', triggerThreshold = 8000)
    while ps5000a.chs['A'].overflow:
        ps5000a.configurePSD(1, int(1953125/5))
        balancer_prism(ps5000a, prism_x, 'A', -1, 0.005, offset = 0, debug = False, rel = 2)
        ps5000a.configurePSD(0.099/avg, int(1953125/5))
        ps5000a.getTimeSignal(trigger = 'D', triggerThreshold = 8000)
    TSidx = [0]
    i = 0
    while i < len(ps5000a.chs['D'].timeSignal) - 6:
        if ps5000a.chs['D'].timeSignal[i+4] - ps5000a.chs['D'].timeSignal[i] > 3:
            temp = np.array(ps5000a.chs['D'].timeSignal[i:i+5])
            for j in range(len(temp)):
                if temp[j] > 3:
                    break
            TSidx = TSidx + [j+i]
            i = i + 5
        else:
            i = i + 1
    print(TSidx)
    idx = np.zeros(len(fcs)*2, dtype = np.int)
    ref = [0]*len(fcs)
    signal = [0]*len(fcs)
    for i in range(len(fcs)):
        idx[2*i] = np.argmin(np.absolute(ps5000a.f - fcs[i] +1000))
        idx[2*i+1] = np.argmin(np.absolute(ps5000a.f - fcs[i] -1000))
    for i in range(avg):
        fftB = fft_process(ps5000a.chs['B'].timeSignal[i*nT:(i+1)*nT], idx, ps5000a.getfs())
        fftC = fft_process(ps5000a.chs['C'].timeSignal[i*nT:(i+1)*nT], idx, ps5000a.getfs())
        #fftB = fft_process(ps5000a.chs['B'].timeSignal[TSidx[i]:TSidx[i+1]], idx)
        #fftC = fft_process(ps5000a.chs['C'].timeSignal[TSidx[i]:TSidx[i+1]], idx)
        for j in range(len(fcs)):
            norm = np.zeros(len(fftB[j]), dtype = np.complex)
            norm = np.cos(np.angle(fftB[j])) + 1j*np.sin(np.angle(fftB[j]))
            if i == 0:
                ref[j] = np.divide(fftB[j], norm)
                signal[j] = np.divide(fftC[j], norm)
            else:
                ref[j] = ref[j] + np.divide(fftB[j], norm)
                signal[j] = signal[j] + np.divide(fftC[j], norm)
    for j in range(len(fcs)):
        signal[j] = signal[j]/avg
        ref[j] = ref[j]/avg
    return signal, ref, idx

def fft_process(TS, idx, fs, window = None):
    fftch = ms.Single_sided_fft(TS, window = window)/np.sqrt(len(TS)*fs)
    ffts = []
    for i in range(int(len(idx)/2)):
        ffts = ffts +[fftch[idx[2*i]:idx[2*i+1]]]
    return ffts

def json_to_dict_171(wd, filename):
    with open(wd+'\\'+filename, 'r') as f:
        target_dict = json.load(f)
    temp = {}
    temp['zs'] = target_dict['zs']
    temp['idxs'] = target_dict['idxs']
    temp['maxfs'] = target_dict['maxfs']
    temp['complexC'] = np.array(target_dict['complexC_Re']) + 1j * np.array(target_dict['complexC_Im'])
    temp['processed_complexC'] = np.array(target_dict['processed_complexC_Re']) + 1j * np.array(target_dict['processed_complexC_Im'])
    temp['fit'] = optimize.OptimizeResult(target_dict['fit'])
    #temp['fit'].direc = np.array(temp['fit'].direc)
    #temp['fit'].x = np.array(temp['fit'].x)
    fftBs = [[]]*len(target_dict['fftBs_Re'])
    fftCs = [[]]*len(target_dict['fftBs_Re'])
    fftpCs = [[]]*len(target_dict['fftBs_Re'])
    for j in range(len(target_dict['fftBs_Re'])):
        fftBs[j] = np.array(target_dict['fftBs_Re'][j]) + 1j*np.array(target_dict['fftBs_Im'][j])
        fftCs[j] = np.array(target_dict['fftCs_Re'][j]) + 1j*np.array(target_dict['fftCs_Im'][j])
        fftpCs[j] = np.array(target_dict['fftpCs_Re'][j]) + 1j*np.array(target_dict['fftpCs_Im'][j])
    temp['fftBs'] = fftBs
    temp['fftCs'] = fftCs
    temp['processed_C'] = fftpCs
    return temp

def json_to_dict_172(wd, filename):
    with open(wd+'\\'+filename, 'r') as f:
        target_dict = json.load(f)
    temp = {}
    temp['zs'] = target_dict['zs']
    temp['idxs'] = target_dict['idxs']
    temp['maxfs'] = target_dict['maxfs']
    temp['fordata'] = target_dict['fordata']
    temp['complexC'] = np.array(target_dict['complexC_Re']) + 1j * np.array(target_dict['complexC_Im'])
    temp['processed_complexC'] = np.array(target_dict['processed_complexC_Re']) + 1j * np.array(target_dict['processed_complexC_Im'])
    temp['fit'] = target_dict['fit']
    #temp['fit'].direc = np.array(temp['fit'].direc)
    #temp['fit'].x = np.array(temp['fit'].x)
    fftBs = [[]]*len(target_dict['fftBs_Re'])
    fftCs = [[]]*len(target_dict['fftBs_Re'])
    fftpCs = [[]]*len(target_dict['fftBs_Re'])
    for j in range(len(target_dict['fftBs_Re'])):
        fftBs[j] = np.array(target_dict['fftBs_Re'][j]) + 1j*np.array(target_dict['fftBs_Im'][j])
        fftCs[j] = np.array(target_dict['fftCs_Re'][j]) + 1j*np.array(target_dict['fftCs_Im'][j])
        fftpCs[j] = np.array(target_dict['fftpCs_Re'][j]) + 1j*np.array(target_dict['fftpCs_Im'][j])
    temp['fftBs'] = fftBs
    temp['fftCs'] = fftCs
    temp['processed_C'] = fftpCs
    return temp

def PIDtofile(t, et, DCs, ps, index, wd):
    t.tofile(wd+'\\'+'time_{:d}'.format(index)+'.bin', sep = '')
    et.tofile(wd+'\\'+'et_{:d}'.format(index)+'.bin', sep = '')
    DCs.tofile(wd+'\\'+'DCs_{:d}'.format(index)+'.bin', sep = '')
    ps.tofile(wd+'\\'+'ps_{:d}'.format(index)+'.bin', sep = '')

def maxs_in_psd_trig(avg, ps5000a, ch, fcs, refch = None):
    psd = ps5000a.getPSD_trig(ch, avg, refch = refch)
    maxs = []
    for i in range(len(fcs)):
        maxs = maxs + [np.max(psd[int(fcs[i]-100):int(fcs[i]+100)])]
    return maxs

def maxs_in_psd_trig_withref(avg, ps5000a, ch, fcs, refch = None):
    psd = ps5000a.getPSD_trig(ch, avg, refch = refch)
    maxs = []
    for i in range(int(len(fcs)/2)):
        maxs = maxs + [np.max(psd[int(fcs[2*i]-100):int(fcs[2*i]+100)] - psd[int(fcs[2*i+1]-100):int(fcs[2*i+1]+100)])]
    return maxs

def maxs_in_psd_trig_withrefch(avg, ps5000a, ch, refch, fcs):
    psd = ps5000a.getPSD_trig(ch, avg)
    maxs = []
    for i in range(int(len(fcs)/2)):
        maxs = maxs + [np.max(psd[int(fcs[2*i]-100):int(fcs[2*i]+100)] - psd[int(fcs[2*i+1]-100):int(fcs[2*i+1]+100)])]
    return maxs

def process_signal_fft_ref(zs, fftBs, fftCs, idxs):
    cpy_fftCs = [[]]*len(fftCs)
    for i in range(len(fftCs)):
        cpy_fftCs[i] = fftCs[i].copy()

    cpy_fftBs = [[]]*len(fftBs)
    for i in range(len(fftBs)):
        cpy_fftBs[i] = fftBs[i].copy()

    index = [int(len(fftCs[0])/2 - 500), int(len(fftCs[0])/2 + 500)]
    complexC = np.zeros(len(cpy_fftCs), dtype = np.complex)
    for i in range(len(cpy_fftCs)):
        maxI = np.argmax(np.absolute(cpy_fftCs[i]))
        indexarra = np.array(list(range(index[0], maxI - 100)) + list(range(maxI + 100, index[1])))
        temp = cpy_fftCs[i]
        complexC[i] = np.mean(temp[indexarra])

    p = np.polyfit(np.real(complexC), np.imag(complexC), 1)
    b = p[1]
    k = p[0]
    phi0 = np.arctan(k)
    A = ((np.real(complexC[0])+k*np.imag(complexC[0]))/np.sqrt(1+k**2) + \
        (np.real(complexC[-1])+k*np.imag(complexC[-1]))/np.sqrt(1+k**2))/(zs[0] - zs[-1])
    outphase_noise = b*np.cos(phi0)

    temp = complexC.copy()
    complexC = complexC - 1j*outphase_noise*np.exp(1j*phi0)
    p = np.polyfit(zs, np.real(complexC), 1)
    z0 = -p[1]/p[0]
    res = {'b':b, 'k':k, 'phi0':phi0, 'A':A, 'z0':z0}

    for i in range(len(cpy_fftCs)):
        cpy_fftCs[i][index[0]:index[1]] = cpy_fftCs[i][index[0]:index[1]] - 1j*outphase_noise*np.exp(1j*phi0)

    dict_re = {'zs_measured':zs, 'fftBs':fftBs, 'fftCs':fftCs, 'processed_C':cpy_fftCs, 'fit':res, 'complexC':temp, 'processed_complexC':complexC, 'idxs':idxs}
    return dict_re

def LineScan_prism(prism_x, Nstep, direction, span, ps5000a):
    stepsize = direction*span/Nstep
    DCA = []
    xs = []
    for i in range(Nstep):
        prism_x.MoveBy(stepsize)
        time.sleep(5)
        ps5000a.configurePSD(10, int(1953125/5))
        ps5000a.AutoRange('A')
        ps5000a.configurePSD(1, int(1953125/5))
        DCA = DCA + [np.mean(ps5000a.getTimeSignal('A'))]
        xs = xs + [prism_x.Position()]
    return xs, DCA

def BKAdata(datas, clusterIs):
    signals_ffts = [[]]
    relative_z = []
    relative_z_data = []
    relative_z_ref = []
    zs = []
    zs_data = []
    fN = 4848
    for i in range(len(datas)):
        for j in range(len(datas[i]['fordata'])):
            relative_z = relative_z + [datas[i]['complexC'][j]]
            zs = zs + [datas[i]['zs'][j]]
            if datas[i]['fordata'][j] in clusterIs:
                zs_data = zs_data + [datas[i]['zs'][j]]
                relative_z_data = relative_z_data + [datas[i]['complexC'][j]]
                #relative_z_data = relative_z_data + [np.mean(datas[i]['fftBs'][j][fN-500:fN+500])]
                signals_ffts = signals_ffts + [datas[i]['fftCs'][j][fN-500:fN+500]]

    del signals_ffts[0]
    p = np.polyfit(np.real(relative_z), np.imag(relative_z), 1)
    b = p[1]
    k = p[0]
    phi0 = np.arctan(k)
    outphase_noise = b*np.cos(phi0)
    check_data = relative_z_data - 1j*outphase_noise*np.exp(1j*phi0)
    check = relative_z - 1j*outphase_noise*np.exp(1j*phi0)
    relative_z_data = (np.real(check_data)+k*np.imag(check_data))/np.sqrt(1+k**2)
    relative_z = (np.real(check)+k*np.imag(check)/np.sqrt(1+k**2))

    inds = np.array(relative_z_data).argsort()
    relative_z_data = np.array(relative_z_data)[inds]
    zs_data = np.array(zs_data)[inds]
    temp = [signals_ffts[i].copy() - 1j*outphase_noise*np.exp(1j*phi0) for i in inds]
    del signals_ffts
    signals = np.zeros([len(temp), len(temp[0])])
    for i in range(len(temp)):
        signals[i] = 10*np.log10(20) + 10*np.log10(np.square(np.absolute(temp[i])))
    return relative_z_data, signals, relative_z, p, check_data, zs_data, zs, check

def BKA_get_data(wd):
    os.chdir(wd)
    filenames = sorted(filter(os.path.isfile, os.listdir('.')), key=os.path.getmtime)
    max = 0
    for name in filenames:
        temp = re.findall('result_[0-9]+_ifs_([0-9]+)\.json',name)
        if len(temp) > 0:
            temp = int(temp[0])
            if temp > max:
                max = temp

    datas_list = [[]]*int(max+1)
    filenames_list = [[]]*int(max+1)
    lengths = [0]*int(max+1)
    for name in filenames:
        temp = re.findall('result_[0-9]+_ifs_([0-9]+)\.json',name)
        if len(temp) > 0:
            idx = int(temp[0])
            lengths[idx] = lengths[idx]+1
            filenames_list[idx] = filenames_list[idx] + [name]

    for i in range(len(datas_list)):
        datas_list[i] = [[]]*lengths[i]

    for j in range(len(filenames_list)):
        for i in range(len(filenames_list[j])):
            datas_list[j][i] = json_to_dict_172(wd, filenames_list[j][i])

    return datas_list

def getResult_fft_sep(avg, ps5000a, fcs, prism_x, windowf = None):
    #only normalize to channel B's phase, ie, fftB is pure real
    ps5000a.configurePSD(0.198, int(1953125/5*3))
    TSC = [[]]*avg
    TSB = [[]]*avg
    for i in range(avg):
        ps5000a.getTimeSignal(trigger = 'D', triggerThreshold = 8000)
        while ps5000a.chs['A'].overflow:
            ps5000a.configurePSD(1, int(1953125/5))
            balancer_prism(ps5000a, prism_x, 'A', -1, 0.005, offset = 0, debug = False, rel = 2)
            ps5000a.configurePSD(0.198, int(1953125/5*3))
            ps5000a.getTimeSignal(trigger = 'D', triggerThreshold = 8000)
        TSC[i] = ps5000a.chs['C'].timeSignal.copy()
        TSB[i] = ps5000a.chs['B'].timeSignal.copy()
    i = 0
    while i < len(ps5000a.chs['D'].timeSignal) - 6:
        if ps5000a.chs['D'].timeSignal[i+4] - ps5000a.chs['D'].timeSignal[i] > 3:
            temp = np.array(ps5000a.chs['D'].timeSignal[i:i+5])
            for j in range(len(temp)):
                if temp[j] > 3:
                    break
            TSidx = j+i
            i = i + 5
        else:
            i = i + 1
    print(TSidx)
    for i in range(avg):
        TSC[i] = TSC[i][:TSidx]
        TSB[i] = TSB[i][:TSidx]
    if windowf  is not None:
        window = windowf(TSidx)
    else:
        window = None
    idx = np.zeros(len(fcs)*2, dtype = np.int)
    ref = [0]*len(fcs)
    signal = [0]*len(fcs)
    psdC = np.zeros(1)
    for i in range(len(fcs)):
        idx[2*i] = np.argmin(np.absolute(ps5000a.f - fcs[i] +1000))
        idx[2*i+1] = np.argmin(np.absolute(ps5000a.f - fcs[i] -1000))
    for i in range(avg):
        fftB = fft_process(TSB[i], idx, ps5000a.getfs(), window = window)
        fftC = fft_process(TSC[i], idx, ps5000a.getfs(), window = window)
        if i == 0:
            psdC = ms.Single_sided_PSD(TSC[i], ps5000a.getfs())
        else:
            psdC = psdC + ms.Single_sided_PSD(TSC[i], ps5000a.getfs())
        #fftB = fft_process(ps5000a.chs['B'].timeSignal[TSidx[i]:TSidx[i+1]], idx)
        #fftC = fft_process(ps5000a.chs['C'].timeSignal[TSidx[i]:TSidx[i+1]], idx)
        for j in range(len(fcs)):
            norm = np.zeros(len(fftB[j]), dtype = np.complex)
            norm = np.cos(np.angle(fftB[j])) + 1j*np.sin(np.angle(fftB[j]))
            if i == 0:
                ref[j] = np.divide(fftB[j], norm)
                signal[j] = np.divide(fftC[j], norm)
            else:
                ref[j] = ref[j] + np.divide(fftB[j], norm)
                signal[j] = signal[j] + np.divide(fftC[j], norm)
    psdC = psdC/avg
    for j in range(len(fcs)):
        signal[j] = signal[j]/avg
        ref[j] = ref[j]/avg
    return signal, ref, idx, psdC

def step_sine(prism_x, ps5000a, FG, fcs, step_size):
    fms = np.linspace(fcs[0] - 10, fcs[0] + 10, int(20/step_size)+1)
    amps = np.zeros(len(fms))
    angles = np.zeros(len(fms))
    t0 = time.time()
    FG.Sine(fms[0], 2)
    time.sleep(30)
    for fm, i in zip(fms, range(len(fms))):
        FG.Sine(fm, 2)
        time.sleep(30)
        ps5000a.getTimeSignal()
        while ps5000a.chs['A'].overflow:
            ps5000a.configurePSD(1, int(1953125/5))
            balancer_prism(ps5000a, prism_x, 'A', -1, 0.005, offset = 0, debug = False, rel = 2)
            ps5000a.configurePSD(5, int(1953125*3))
            ps5000a.getTimeSignal()
        amps[i], angles[i] = lock_in_amplifier(ps5000a.chs['A'].timeSignal, ps5000a.chs['B'].timeSignal, ps5000a.getfs(), fm)
        print(time.time()-t0)
        print(amps[i], angles[i])
    return amps, angles, fms

def json_to_dict_cali(wd, filename):
    with open(wd+'\\'+filename, 'r') as f:
        target_dict = json.load(f)
    temp = {}
    temp['zs'] = target_dict['zs']
    temp['idxs'] = target_dict['idxs']
    temp['maxfs'] = target_dict['maxfs']
    temp['fordata'] = target_dict['fordata']
    temp['complexC'] = np.array(target_dict['complexC_Re']) + 1j * np.array(target_dict['complexC_Im'])
    temp['processed_complexC'] = np.array(target_dict['processed_complexC_Re']) + 1j * np.array(target_dict['processed_complexC_Im'])
    temp['fit'] = target_dict['fit']
    temp['psdCs'] = target_dict['psdCs']
    #temp['fit'].direc = np.array(temp['fit'].direc)
    #temp['fit'].x = np.array(temp['fit'].x)
    fftBs = [[]]*len(target_dict['fftBs_Re'])
    fftCs = [[]]*len(target_dict['fftBs_Re'])
    fftpCs = [[]]*len(target_dict['fftBs_Re'])
    for j in range(len(target_dict['fftBs_Re'])):
        fftBs[j] = np.array(target_dict['fftBs_Re'][j]) + 1j*np.array(target_dict['fftBs_Im'][j])
        fftCs[j] = np.array(target_dict['fftCs_Re'][j]) + 1j*np.array(target_dict['fftCs_Im'][j])
        fftpCs[j] = np.array(target_dict['fftpCs_Re'][j]) + 1j*np.array(target_dict['fftpCs_Im'][j])
    temp['fftBs'] = fftBs
    temp['fftCs'] = fftCs
    temp['processed_C'] = fftpCs
    return temp

def lock_in_amplifier_temp(signal, ref, fs, freq):
    lowcut = freq - 10
    highcut = freq + 10
    T = len(signal)/fs
    ref1 = ms.butter_bandpass_filter(ref, lowcut, highcut, fs, order=1)
    ref2 = ms.HilbertTransform(ref1)
    X = np.dot(signal, ref1)/T
    Y = np.dot(signal, ref2)/T
    amp = np.absolute(np.complex(X, Y))
    angle = np.angle(np.complex(X, Y))
    return amp, angle
