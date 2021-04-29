import lib.ps5000A as ps
import lib.FuncGen33522b as FuncGen
import lib.mathlib_shan as ms
import lib.splitPair as sp
import matplotlib.pyplot as plt
import numpy as np
import time
import sys
sys.path.append(r'C:\Program Files\Newport\Piezo Motion Control\Newport AG-UC2-UC8 Applet\Bin')
import clr
clr.AddReference('AgilisCmdLib')
from Newport.AgilisCmdLib import *
import lib.chirp_scan_lib as CH

wd = r'Z:\data\optical lever project\NORCADA_NX53515C\89-LineScan'
ps5000a = ps.picoscope(wd)
FG = FuncGen.FuncGen()
fstart = 182350
fend = 182380
amp = 0.25
dc = 0.365
AWGfs = 3.2e6
AWGN = 12.8e6
# initialize AG-UC2
UC2 = AgilisCmds()
DevicePorts = UC2.GetDevices()
UC2.OpenInstrument(DevicePorts[2])
controllerAddress = int(1);
err = ''
UC2.MR(err)
#
safeRangeA = [-1, 11]
safeRangeB = [-1, 11]
DS = sp.splitPair(safeRangeA , safeRangeB)
DS.Aangle = 70.1666/180*np.pi
DS.Bangle = 77.7557/180*np.pi
#
fstart = 182270
fend = 182290
#arbname = 'f={fstart:.0f}_{fend:.0f}'.format(fstart = fstart, fend = fend)
FG.chirp('test', fstart, fend, AWGfs, AWGN, dc+amp, dc-amp)

ps5000a.defaultSetting()
ps5000a.AC('A')
ps5000a.AutoRange('A')
ps5000a.AC('B')
ps5000a.AutoRange('B')
ps5000a.AC('C')
ps5000a.AutoRange('C')
ps5000a.DC('D')
ps5000a.ChangeRangeTo('D', 5000)
ps5000a.ChangeRangeTo('A', 2000)
factor = 1
ps5000a.configurePSD(0.25/factor, 1953125) #gives sampling intervals of 256ns

def chirp_check(fcs, fspans, df, avg, ps5000a, FG, FGwrite, AWGF, filename = None, wd = None,  amps = 0.25, dc = 0.365, AWGfs = 3.2e6, AWGN = 12.8e6, offset = None):
    if FGwrite:
        if AWGF == 'chirp':
            FG.chirps('test', fcs - fspans/2, fcs + fspans/2, AWGfs, AWGN, dc+amps, dc-amps)
        elif AWGF == 'sweep':
            FG.sweep(fcs[0]-fspans[0]/2, fcs[0]-fspan[0]/2, 1/df, amps[0], dc)
    factor = 0.25/df
    ps5000a.configurePSD(0.25/factor, 1953125) #gives sampling intervals of 256ns
    #offset = int(-15*factor)
    if offset is None:
        offset = -int(np.round(1e3/0.930*(1/df)/ps5000a.timeInternalns.value))
    ps5000a.getTimeSignal()
    psd = ps5000a.PSDfromTS(ps5000a.chs['A'].timeSignal[0:offset], ps5000a.getfs())

    f = ps5000a.f[0:len(psd)]
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
    powers, psdAs, psdCs= Piecewise(fcs, fspans, dftAs[i], dftCs[i])

    if filename is not None:
        dftAs.tofile(wd+'\\'+'BalDFT_{filename:s}'.format(filename = filename)+'.bin', sep = '')
        dftBs.tofile(wd+'\\'+'LaserPDFT_{filename:s}'.format(filename = filename)+'.bin', sep = '')
        dftCs.tofile(wd+'\\'+'DriveDFT_{filename:s}'.format(filename = filename)+'.bin', sep = '')
        psdavg.tofile(wd+'\\'+'PSD_{filename:s}'.format(filename = filename)+'.bin', sep = '')
        drivepsd.tofile(wd+'\\'+'DrivePSD_{filename:s}'.format(filename = filename)+'.bin', sep = '')
    for psdavg, drivepsd in zip(psdAs, psdCs)
        plt.plot(f, psdavg)
        plt.plot(f, drivepsd)
        plt.show()

def sweep(FG, ps5000a, UC2, Nsteps):
    wd = r'Z:\data\optical lever project\NORCADA_NX53515C\79-chirp_over_dy'
    err = ''
    oldpeakf = 0
    peakf = oldpeakf
    fspan = 20
    f1 = peakf - fspan/2
    f2 = peakf + fspan/2
    dc = 0.365
    amp = 0.25
    for i in range(Nsteps):
        [tsAs, tsBs, tsCs] = get_chirp(0.25, 5, ps5000a)
        UC2.PR(int(1), 1, err)
        avg = np.min(tsAs.shape)

        ps5000a.PSDfromTS(ps5000a.chs['A'].timeSignal[0:offset], ps5000a.getfs())
        check1 = np.argmin(np.absolute(ps5000a.f-182200))
        check2 = np.argmin(np.absolute(ps5000a.f-182600))
        f = ps5000a.f[check1:check2]

        dftAs = np.zeros((avg, len(f)), dtype = np.complex)
        dftBs = np.zeros((avg, len(f)), dtype = np.complex)
        dftCs = np.zeros((avg, len(f)), dtype = np.complex)
        for j in range(avg):
            dftA = ms.Single_sided_fft(tsAs[j,:])
            if j == 0:
                maxI = np.argmax(np.absolute(dftA))
                peakf = f[maxI]
                if np.absolute(peakf - oldpeakf)>3:
                    f1 = peakf - fspan/2
                    f2 = peakf + fspan/2
                    FG.chirp('test', f1, f2, 3.2e6, 12.8e6, dc+amp, dc-amp)
                    oldpeakf = peakf
            dftB = ms.Single_sided_fft(tsBs[j,:])
            dftC = ms.Single_sided_fft(tsCs[j,:])
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
    f.tofile(wd+'\\'+'freq'+'.bin', sep = '')


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

    fstarts = fcs - fspans/2
    fends = fcs + fspans/2

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

    check1s = np.zeros(len(fcs))
    check2s = np.zeros(len(fcs))
    for f1, f2, i in zip(fstarts, fends, list(range(len(fcs))))
        check1s[i] = np.argmin(np.absolute(ps5000a.f-f1+100))
        check2s[i] = np.argmin(np.absolute(ps5000a.f-f2-100))

    powers = np.zeros(len(UC2_relative_xs), len(fcs))
    gamma = 20

    dftAs = [[]]*len(fcs)
    dftCs = [[]]*len(fcs)
    freqs = [[]]*len(fcs)

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

        for j in range(len(fcs)):
            dftAs[i] = dftA[check1s[i]:check2s[i]]
            dftCs[i] = dftC[check1s[i]:check2s[i]]
            freqs[i] = ps5000a.f[check1s[i]:check2s[i]]

        for j in range(len(fcs)):
            psdAs[j] = 10*np.log10(np.square(np.absolute(dftAs[j])))
            psdCs[j] = 10*np.log10(np.square(np.absolute(dftCs[j])))
        maxIs = np.zeros(len(fcs))
        popts = [[]]*len(fcs)
        for j in range(len(fcs)):
            maxIs[j] = np.argmax(psdAs[j])
            popts[j], pcov = curve_fit(fitFunc, f[j][maxIs[j]-300:maxIs[j]+100], psdA[j][maxIs[j]-300:maxIs[j]+100] - psdC[j][maxIs[j]-300:maxIs[j]+100] , p0 =[gamma, f[j][maxIs[j]],  180], \
                bounds = ([0.99*gamma[j], f[j][maxIs[j]] - 50,  120],[1.01*gamma[j], f[j][maxIs[j]]+50,  210]))
            if i == 0:
                gamma[j] = popt[j][0]
            powers[i, j] = popts[j][2]
        print(powers[i])

        t1 = time.time()
        time.sleep(max(waitTime_after_moving - (t1 - t0), 0))
        ps5000a.configurePSD(8, 1953125)
        DS.AutoBalance(ps5000a)
    ps5000a.f.tofile(wd+'\\'+'freq.bin', sep = '')
    powers.tofile(wd+'\\'+'powers.bin', sep = '')
    (np.array(xs)).tofile(wd+'\\'+'UC2steps.bin', sep = '')
    ps5000a.configurePSD(0.25, 1953125)
    plt.plot(UC2_relative_xs, powers)
    plt.show()

def Piecewise(fcs, fspans, dftA, dftC):
    fstarts = fcs - fspans/2
    fends = fcs + fspans/2
    Nj = len(fcs)
    psdA = 10*np.log10(np.square(np.absolute(dftA)))
    psdC = 10*np.log10(np.square(np.absolute(dftC)))
    powers = np.zeros(Nj)

    check1s = np.zeros(Nj)
    check2s = np.zeros(Nj)
    for f1, f2, j in zip(fstarts, fends, list(range(Nj)))
        check1s[j] = np.argmin(np.absolute(ps5000a.f-f1+100))
        check2s[j] = np.argmin(np.absolute(ps5000a.f-f2-100))

    psdAs = [[]]*Nj
    psdCs = [[]]*Nj
    freqs = [[]]*Nj

    for j in range(Nj):
        psdAs[j] = psdA[check1s[j]:check2s[j]]
        psdCs[j] = psdC[check1s[j]:check2s[j]]
        freqs[j] = ps5000a.f[check1s[j]:check2s[j]]

    maxIs = np.zeros(Nj)
    popts = [[]]*Nj
    for j in range(Nj):
        maxIs[j] = np.argmax(psdAs[j])
        popts[j], pcov = curve_fit(fitFunc, f[j][maxIs[j]-300:maxIs[j]+100], psdA[j][maxIs[j]-300:maxIs[j]+100] - psdC[j][maxIs[j]-300:maxIs[j]+100] , p0 =[gamma, f[j][maxIs[j]],  180], \
            bounds = ([0.99*gamma[j], f[j][maxIs[j]] - 50,  120],[1.01*gamma[j], f[j][maxIs[j]]+50,  210]))
        if i == 0:
            gamma[j] = popt[j][0]
        powers[j] = popts[j][2]
    return powers, psdAs, psdCs
