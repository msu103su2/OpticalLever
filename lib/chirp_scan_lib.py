import lib.mathlib_shan as ms
import numpy as np
import os
import time
import matplotlib.pyplot as plt
from scipy import signal
from scipy.optimize import curve_fit

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

def LineScan(f1, f2, UC2, ps5000a, DS, UC2Amp, UC2_relative_xs, UC2_controllerAddress, wd):
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

    #print('connect')
    #input("Press Enter to continue...")
    UC2.PR(UC2_controllerAddress, steps[0], err)
    #print('Disconnect')
    #input("Press Enter to continue...")
    time.sleep(30)
    ps5000a.configurePSD(8, 1953125)
    DS.AutoBalance(ps5000a, 1)

    ps5000a.configurePSD(1, 1953125)
    ps5000a.AC('A')
    ps5000a.ChangeRangeTo('A', 2000)
    ps5000a.AC('B')
    ps5000a.AutoRange('B')
    ps5000a.AC('C')
    ps5000a.AutoRange('C')

    ps5000a.configurePSD(0.25, 1953125)
    temp = ps5000a.getTimeSignal('A')
    ps5000a.PSDfromTS(temp, ps5000a.getfs())
    temp, _, _  = get_chirp(0.25, avg, ps5000a)
    temp = ms.Single_sided_fft(temp[0,:])
    range_A = 2000

    check1 = np.argmin(np.absolute(ps5000a.f-f1+100))
    check2 = np.argmin(np.absolute(ps5000a.f-f2-100))

    powers = np.zeros(len(UC2_relative_xs))
    gamma = 20

    for i in range(len(UC2_relative_xs)):
        ps5000a.ChangeRangeTo('A', range_A)
        cache = np.zeros(len(temp))
        ps5000a.configurePSD(0.25, 1953125)
        [tsAs, tsBs, tsCs] = get_chirp(0.25, avg, ps5000a, trigger = 'D', triggerThreshold = 8000)
        window = signal.windows.tukey(len(tsAs[0,:]), 0.8)

        dftA = np.zeros(len(temp))
        dftB = np.zeros(len(temp))

        t0 = time.time()
        #print('connect')
        #input("Press Enter to continue...")
        UC2.PR(UC2_controllerAddress, steps[i+1], err)
        #print('Disconnect')
        #input("Press Enter to continue...")

        for j in range(avg):
            dftA = ms.Single_sided_fft(tsAs[j,:], window)
            dftC = ms.Single_sided_fft(tsCs[j,:], window)
            cache = cache + np.multiply(np.divide(dftA, dftC), np.absolute(dftC))
        cache = cache/avg
        cache = np.square(np.absolute(cache))
        psdC = np.square(np.absolute(dftC))
        cache.tofile(wd+'\\'+'signal_{i:d}.bin'.format(i = i), sep = '')
        psdC.tofile(wd+'\\'+'drive_{i:d}.bin'.format(i = i), sep = '')
        print('Progress:{p:.2%}'.format(p = (i+1)/(len(UC2_relative_xs)+1)))

        dftA = dftA[check1:check2]
        dftC = dftC[check1:check2]
        f = ps5000a.f[check1:check2]

        psdA = np.square(np.absolute(dftA))
        I = np.argmax(psdA)
        psdA = 10*np.log10(psdA)
        psdC = 10*np.log10(np.square(np.absolute(dftC)))
        popt, pcov = curve_fit(fitFunc, f[I-300:I+100], psdA[I-300:I+100] - psdC[I-300:I+100] , p0 =[gamma, f[I],  180], \
            bounds = ([0.99*gamma, f[I] - 50,  120],[1.01*gamma, f[I]+50,  210]))
        if i == 0:
            gamma = popt[0]
        powers[i] = popt[2]
        print(popt[2])

        t1 = time.time()
        time.sleep(max(waitTime_after_moving - (t1 - t0), 0))
        ps5000a.configurePSD(8, 1953125)
        DS.AutoBalance(ps5000a)
    ps5000a.f.tofile(wd+'\\'+'freq.bin', sep = '')
    ps5000a.configurePSD(0.25, 1953125)
    plt.plot(UC2_relative_xs, powers)
    plt.show()

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

    if FGwrite:
        if AWGF == 'chirp':
            FG.chirps('test', fcs - fspans/2, fcs + fspans/2, AWGfs, AWGN, dc+amps, dc-amps)
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
    for i in range(avg):
        dftAs[i] = ms.Single_sided_fft(tsAs[i,:], window)
        dftBs[i] = ms.Single_sided_fft(tsBs[i,:], window)
        dftCs[i] = ms.Single_sided_fft(tsCs[i,:], window)
        print(i)
    powers, gamma, psdAs, psdCs, freqs, popts= Piecewise(fcs, fspans, dftAs[i], dftBs[i], f)

    if filename is not None:
        dftAs.tofile(wd+'\\'+'BalDFT_{filename:s}'.format(filename = filename)+'.bin', sep = '')
        dftBs.tofile(wd+'\\'+'LaserPDFT_{filename:s}'.format(filename = filename)+'.bin', sep = '')
        dftCs.tofile(wd+'\\'+'DriveDFT_{filename:s}'.format(filename = filename)+'.bin', sep = '')
        psdavg.tofile(wd+'\\'+'PSD_{filename:s}'.format(filename = filename)+'.bin', sep = '')
        drivepsd.tofile(wd+'\\'+'DrivePSD_{filename:s}'.format(filename = filename)+'.bin', sep = '')
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

def Piecewise(fcs, fspans, dftA, dftC, f, gamma = None):
    fcs = np.array(fcs)
    fspans = np.array(fspans)
    fstarts = fcs - fspans/2
    fends = fcs + fspans/2
    Nj = len(fcs)
    psdA = 10*np.log10(np.square(np.absolute(dftA)))
    psdC = 10*np.log10(np.square(np.absolute(dftC)))
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
        if gamma_provided:
            popts[j], pcov = curve_fit(fitFunc, freqs[j][maxIs[j]-300:maxIs[j]+100], psdAs[j][maxIs[j]-300:maxIs[j]+100] - psdCs[j][maxIs[j]-300:maxIs[j]+100] , p0 =[gamma[j], freqs[j][maxIs[j]],  180], \
                bounds = ([0.99*gamma[j], freqs[j][maxIs[j]] - 50,  120],[1.01*gamma[j], freqs[j][maxIs[j]]+50,  270]))
        else:
            popts[j], pcov = curve_fit(fitFunc, freqs[j][maxIs[j]-300:maxIs[j]+100], psdAs[j][maxIs[j]-300:maxIs[j]+100] - psdCs[j][maxIs[j]-300:maxIs[j]+100] , p0 =[gamma[j], freqs[j][maxIs[j]],  180], \
                bounds = ([0, freqs[j][maxIs[j]] - 50,  120],[300, freqs[j][maxIs[j]]+50,  270]))
        gamma[j] = popts[j][0]
        powers[j] = popts[j][2]
    return powers, gamma, psdAs, psdCs, freqs, popts
def LineScan_new(fcs, fspans, UC2, ps5000a, DS, UC2Amp, UC2_relative_xs, UC2_controllerAddress, wd):
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
    ps5000a.ChangeRangeTo('A', 2000)
    ps5000a.AC('B')
    ps5000a.AutoRange('B')
    ps5000a.AC('C')
    ps5000a.AutoRange('C')

    ps5000a.configurePSD(0.25, 1953125)
    ps5000a.getTimeSignal()
    ps5000a.PSDfromTS(ps5000a.chs['A'].timeSignal, ps5000a.getfs())
    temp, _, _  = get_chirp(0.25, avg, ps5000a)
    temp = ms.Single_sided_fft(temp[0,:])
    f = ps5000a.f[0:len(temp)]
    range_A = 2000
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
        UC2.PR(UC2_controllerAddress, steps[i+1], err)

        for j in range(avg):
            dftA = ms.Single_sided_fft(tsAs[j,:], window)
            dftC = ms.Single_sided_fft(tsCs[j,:], window)
            cache = cache + np.multiply(np.divide(dftA, dftC), np.absolute(dftC))
        cache = cache/avg
        cache = np.square(np.absolute(cache))
        psdC = np.square(np.absolute(dftC))
        cache.tofile(wd+'\\'+'signal_{i:d}.bin'.format(i = i), sep = '')
        psdC.tofile(wd+'\\'+'drive_{i:d}.bin'.format(i = i), sep = '')
        print('Progress:{p:.2%}'.format(p = (i+1)/(len(UC2_relative_xs))))

        if i == 0:
            powers, gamma, psdAs, psdCs, freqs, popts = Piecewise(fcs, fspans, dftA, dftC, f, gamma = None)
        else:
            powers, _, psdAs, psdCs, freqs, popts = Piecewise(fcs, fspans, dftA, dftC, f, gamma)
        powers_on_step[i] = powers
        print(powers)

        t1 = time.time()
        time.sleep(max(waitTime_after_moving - (t1 - t0), 0))
        ps5000a.configurePSD(8, 1953125)
        DS.AutoBalance(ps5000a)

    powers_on_step = np.array(powers_on_step)
    for i in range(len(powers_on_step[0])):
        plt.plot(UC2_relative_xs, powers_on_step[:,i])
    plt.show()
    ps5000a.f.tofile(wd+'\\'+'freq.bin', sep = '')
    (np.array(UC2_relative_xs)).tofile(wd+'\\'+'UC2steps.bin', sep = '')
    (np.array(powers_on_step)).tofile(wd+'\\'+'powers.bin', sep = '')
    ps5000a.configurePSD(0.25, 1953125)
