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
nT = int(39063235/5)

wd = r'Z:\data\optical lever project\NORCADA_NX53515C\manuscript_data'
ps5000a = ps.picoscope()
FG = FuncGen.FuncGen()
# initialize AG-UC2
UC2 = AgilisCmds()
DevicePorts = UC2.GetDevices()
UC2.OpenInstrument('COM4')
controllerAddress = int(1);
err = ''
UC2.MR(err)
#
safeRangeA = [-1, 13]
safeRangeB = [-1, 13]
DS = sp.splitPair(safeRangeA , safeRangeB)
DS.Aangle = 70.1666/180*np.pi
DS.Bangle = 77.7557/180*np.pi
#

ps5000a.defaultSetting()
ps5000a.AC('A')
ps5000a.DC('B')
ps5000a.DC('C')
ps5000a.DC('D')
ps5000a.ChangeRangeTo('D', 5000)
ps5000a.ChangeRangeTo('A', 2000)
ps5000a.ChangeRangeTo('B', 2000)
ps5000a.ChangeRangeTo('C', 5000)
ps5000a.configurePSD(0.25, 1953125) #gives sampling intervals of 256ns
fcs = [182300, 260900, 315200]
amps = [0.005]*3
fspans = [1e3]*3
factor = 1

fcs = [146798, 279959, 602338]

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

monitorTime = 3*3600
wd = wd = r'Z:\data\optical lever project\NORCADA_NX53515C\145-fdrift'
fcs = [320000, 530680, 580680]
amps = [0.5, 1, 0]
fspans = [80, 80, 80]
sinef = 50000
I0 = [180000*4, 259000*4, 529500*4, 580400*4]
I1 = [184500*4, 263000*4, 532000*4, 581200*4]
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
    if np.absolute(data[3][-1] - fcs[2])>fspans[2]/2*0.75:
        fcs[1] = data[3][-1] - sinef
        fcs[2] = data[3][-1]
        ref, sig = CH.chirp_check_cm(fcs, fspans, 0.25, 1, ps5000a, FG, 1, 'chirp', amps = amps, refOff = 1, dc = 0)

del times[0]
np.array(times).tofile(wd+'\\'+'times.bin', sep = '')
for i in range(len(fcs)):
    np.array(data[i]).tofile(wd+'\\'+'fm{i:d}.bin'.format(i = i), sep = '')
t0 = time.time()
test = ps5000a.getPSD('A', offset = -18, trigger = 'D', avg = 10)
t1 = time.time()
print(t0-t1)

wd = r'Z:\data\optical lever project\NORCADA_NX53515C\155-90degree'
monitorTime = 1200
fcs = [132750, 250400, 343262]
amps = [0.005]*3
fspans = [1e3]*3
def fdrift(fcs, fspans, wd, monitorTime, ps5000a):
    if not os.path.exists(wd):
        os.makedirs(wd)
    I0 = ((np.array(fcs) - np.array(fspans)/2)*4).astype('i')
    I1 = ((np.array(fcs) + np.array(fspans)/2)*4).astype('i')
    data = [[]]*len(I0)
    times = [[]]
    LP = [[]]
    t0 = time.time()
    t1 = time.time()
    while t1 - t0 < monitorTime:
        test = ps5000a.getPSD('A', avg = 1, offset = -18, trigger = 'D', triggerThreshold = 8000)
        t1 = time.time()
        print(t1-t0)
        times = times + [t1]
        LP = LP + [np.mean(ps5000a.chs['B'].timeSignal)]
        for i in range(len(I0)):
            data[i] = data[i]+[(np.argmax(test[I0[i]:I1[i]])+I0[i])/4]
    del times[0]
    np.array(times).tofile(wd+'\\'+'times.bin', sep = '')
    for i in range(len(I0)):
        np.array(data[i]).tofile(wd+'\\'+'fm{i:d}.bin'.format(i = i), sep = '')
    for i in range(len(fcs)):
        plt.plot(times, (data[i] - np.mean(data[i]))/np.mean(data[i])*1e5)
    return times, data, LP

def thorTH(T):
    A = -1.547e1
    B = 5.602e3
    C = -3.789e5
    D = 2.497e7
    R0 = 1e4
    R = R0*np.exp(A+B/T+C/T**2+D/T**3)
    return R

def proTH(T):
    R0 = 1e4
    B = 3870.447
    T0 = 298.15
    R = R0*np.exp(B*(1/T - 1/T0))
    return R

popt, pcov = curve_fit(func, T, tR, p0 =[1e4, 1e3, 298.15]))
t0 = time.time()
for i in range(100):
    pro8000.getTa()


t1 = time.time()

def PIDtest(Ts, P, I, D, timespan, wd):
    pro8000.setP(P)
    pro8000.setI(I)
    pro8000.setD(D)
    pro8000.setTs(Ts)
    print("Ts = {Ts:f}, P = {P:f}, I = {I:f}, D = {D:f}".format(Ts = pro8000.getTs(), P = pro8000.getP(), I = pro8000.getI(), D = pro8000.getD()))
    T = []
    times = []
    t0 = time.time()
    t1 = time.time()
    while t1-t0 < timespan:
        T = T+[pro8000.getTa()]
        times = times + [float(t1) - float(t0)]
        t1 = time.time()
        print(T[len(T) - 1])
        time.sleep(1)
    T = np.array(T)
    times = np.array(times)
    times.tofile(wd+'\\'+'times_T.bin', sep = '')
    T.tofile(wd+'\\'+'T.bin', sep = '')
    Tset = np.zeros(len(times))+Ts
    return times, T, P, I, D, Tset

import lib.CONEX_CC as cc
import lib.ps5000A as ps
import sys
import numpy as np
sys.path.append(r'C:\Program Files\Newport\Piezo Motion Control\Newport AG-UC2-UC8 Applet\Bin')
import clr
clr.AddReference('AgilisCmdLib')
from Newport.AgilisCmdLib import *
import lib.mathlib_shan as ms
import lib.chirp_scan_lib as CH

axis_z = cc.CONEX_CC('ASRL3::INSTR', [-1,26], False)
prism_x = cc.CONEX_CC('ASRL6::INSTR', [-1,26], False)
UC8 = AgilisCmds()
DevicePorts = UC8.GetDevices()
UC8.OpenInstrument('COM9')
controllerAddress = int(1);
err = ''
UC8.MR(err)
UC8.CC_Set(2, err)
axis_x = int(1)
axis_y = int(2)

ps5000a = ps.picoscope()
ps5000a.defaultSetting()
ps5000a.DC('A')
ps5000a.DC('B')
ps5000a.DC('C')
ps5000a.DC('D')
ps5000a.ChangeRangeTo('A', 200)
ps5000a.ChangeRangeTo('B', 10000)
ps5000a.ChangeRangeTo('C', 200)
ps5000a.ChangeRangeTo('D', 5000)
ps5000a.configurePSD(10, 1953125)

HAs = []
HCs = []
zs = []
f = [49500, 50500]
ps5000a.getTimeSignal()
ps5000a.PSDfromTS(ps5000a.chs['A'].timeSignal, ps5000a.getfs())
idx = [0,0]
idx[0] = np.argmin(np.absolute(ps5000a.f-f[0]))
idx[1] = np.argmin(np.absolute(ps5000a.f-f[1]))
z = axis_z.Position()
avg = 20
while z < 22.6:
    axis_z.MoveBy(0.2)
    temp = CH.balancer(ps5000a, UC8, axis_x, 'C', 1, 10)
    temp = CH.balancer(ps5000a, UC8, axis_y, 'A', -1, 50)
    for i in range(avg):
        ps5000a.getTimeSignal()
        if i == 0:
            psdA = ms.Single_sided_PSD(ps5000a.chs['A'].timeSignal, ps5000a.getfs())
        else:
            psdA = psdA + ms.Single_sided_PSD(ps5000a.chs['A'].timeSignal, ps5000a.getfs())
    psdA = psdA / avg
    HA = np.max(psdA[idx[0]:idx[1]])
    HAs = HAs + [HA]
    for i in range(avg):
        ps5000a.getTimeSignal()
        if i == 0:
            psdC = ms.Single_sided_PSD(ps5000a.chs['C'].timeSignal, ps5000a.getfs())
        else:
            psdC = psdC + ms.Single_sided_PSD(ps5000a.chs['C'].timeSignal, ps5000a.getfs())
    psdC = psdC / avg
    HC = np.max(psdC[idx[0]:idx[1]])
    HCs = HCs + [HC]
    z = axis_z.Position()
    print(z)
    zs = zs + [z]

def balancer(ps5000a, UC2, axis, ch, direction, amp):
    imb = np.mean(ps5000a.getTimeSignal(ch))
    while np.absolute(imb) > 3e-3:
         UC2.PR(axis, np.sign(imb)*direction*50, err)
         imb = np.mean(ps5000a.getTimeSignal(ch))
    return imb

CH.balancer(ps5000a, UC8, axis_x, 'C', 1, 4)
CH.balancer(ps5000a, UC8, axis_y, 'A', -1, 4)

fcs = [135450, 144900, 145275, 160975, 169300, 179300, 184100, 189925, 194275, 198575, 199237]
fcs = [132070, 144529, 145039, 160605, 169075, 178935, 183493, 189350, 194025, 198046, 198811]
fspans = [400]*len(fcs)
amps = [0.005]*len(fcs)

[131675.0, 144171, 144796, 160224, 168839, 178549, 182899, 188783, 193744, 197519, 198391]
fcs = [135450, 144900, 145275, 160975, 169300, 179300, 184100, 189925, 194275, 198575, 199237]

ps5000a.AC('B')
ps5000a.ChangeRangeTo('B', 500)
ps5000a.ChangeRangeTo('A', 500)
ps5000a.ChangeRangeTo('C', 500)
RFf = np.linspace(79.9, 80, 2)
RFp = np.linspace(-5, -3, 9)
mZs = []
mBs = []
mCs = []
mAs = []
ops = []
z_start = 0.2
z_end = 5

m = len(RFp)
zs_all = [[]]*m
HAs_all = [[]]*m
HBs_all = [[]]*m
HCs_all = [[]]*m
i= 0
for f in RFf.tolist():
#for p in RFp.tolist():
    E4436B.write(':FREQ:CW {f:.2f}MHz'.format(f = f))
    #E4436B.write(':POW:AMPL {p:.2f}dBm'.format(p = p))
    ps5000a.AC('B')
    ps5000a.ChangeRangeTo('B', 500)
    if z_start > 25 or z_end < 0:
        z_start = 0.2
        z_end = 5
        continue
    zs_all[i], HAs_all[i], HBs_all[i], HCs_all[i] = CH.z_sweep(axis_z, ps5000a, UC8, z_start, z_end, 0.2)
    fit = fitparaballic(zs, HCs)
    min_z = -fit[1]/(2*fit[0])
    if min_z > axis_z.safeRange[0] and min_z < axis_z.safeRange[1]:
        axis_z.Quasi_MoveTo(min_z)
        HA, HB, HC = CH.max_in_psd(20, ps5000a)
    else:
        HA = 0
        HB = 0
        HC = fit[2]-fit[1]**2/(4*fit[0])
    mZs = mZs + [min_z]
    mAs = mAs + [HA]
    mBs = mBs + [HB]
    mCs = mCs + [HC]
    ps5000a.DC('B')
    ps5000a.ChangeRangeTo('B', 10000)
    ops = ops + [np.mean(ps5000a.getTimeSignal('B'))]
    z_start = max([0.2, min_z - 4])
    z_end = min([24.2, min_z + 4])
    i = i+1
mCs_old = mCs

ps5000a.AC('B')
ps5000a.ChangeRangeTo('B', 500)
ps5000a.ChangeRangeTo('A', 500)
ps5000a.ChangeRangeTo('C', 500)
RFp = np.linspace(5, 10, 11)
mZs = []
mBs = []
mCs = []
mAs = []
ops = []
z_start = 0.2
z_end = 8
for p in RFp.tolist():
    E4436B.write(':POW:AMPL {p:.2f}dBm'.format(p = p))
    ps5000a.AC('B')
    ps5000a.ChangeRangeTo('B', 500)
    zs, HAs, HBs, HCs = CH.z_sweep(axis_z, ps5000a, UC8, z_start, z_end, 0.2)
    fit = fitparaballic(zs, HCs)
    min_z = -fit[1]/(2*a)
    if min_z > axis_z.safeRange[0] and min_z < axis_z.safeRange[1]:
        axis_z.Quasi_MoveTo(min_z)
        HA, HB, HC = CH.max_in_psd(20, ps5000a)
    else:
        HA = 0
        HB = 0
        HC = fit[2]-fit[1]**2/(4*a)
    mZs = mZs + [min_z]
    mAs = mAs + [HA]
    mBs = mBs + [HB]
    mCs = mCs + [HC]
    ps5000a.DC('B')
    ps5000a.ChangeRangeTo('B', 10000)
    ops = ops + [np.mean(ps5000a.getTimeSignal('B'))]

UC8.PR(axis_y, 500, err)
while status is not 0:
    _, status,_ = UC8.TS(axis_y, status, err)
    print(status)

10.17
ps5000a.configurePSD(0.09, 1953125)
nT = 39063235
avg = 10
for i in range(avg)
ps5000a.getTimeSignal(trigger = 'D', triggerThreshold = 8000)
ref = ms.Single_sided_PSD(ps5000a.chs['B'].timeSignal[0:nT], ps5000a.getfs())
ref = 10*np.log10(20*ref)
signal = ms.Single_sided_PSD(ps5000a.chs['C'].timeSignal[0:nT], ps5000a.getfs())
signal = 10*np.log10(20*signal)
plt.plot(np.linspace(0,len(ref)/10, len(ref)), signal-ref)

def getResult(avg, ps5000a, fc):
    nT = int(39063235/5)
    for i in range(avg):
        ps5000a.getTimeSignal(trigger = 'D', triggerThreshold = 8000)
        if i == 0:
            ref = ms.Single_sided_PSD(ps5000a.chs['B'].timeSignal[0:nT], ps5000a.getfs())
            signal = ms.Single_sided_PSD(ps5000a.chs['C'].timeSignal[0:nT], ps5000a.getfs())
        else:
            ref = ref + ms.Single_sided_PSD(ps5000a.chs['B'].timeSignal[0:nT], ps5000a.getfs())
            signal = signal + ms.Single_sided_PSD(ps5000a.chs['C'].timeSignal[0:nT], ps5000a.getfs())
    ref = ref/avg
    signal = signal/avg
    ref = 10*np.log10(20*ref)
    signal = 10*np.log10(20*signal)
    result = signal-ref
    result = result[(fc-100)*10:(fc+100)*10]
    return result



plt.plot(np.linspace(0,len(ref)/10, len(ref)), signal-ref)
plt.show()

steps = [0.01]*100
psds = [[]]*len(steps)
for step,i in zip(steps, range(len(steps))):
    axis_z.MoveBy(step)
    ps5000a.DC('A')
    ps5000a.DC('C')
    ps5000a.configurePSD(3, int(1953125/5))
    CH.QPD_T_balance(ps5000a, UC8)
    ps5000a.AC('A')
    ps5000a.AC('C')
    ps5000a.configurePSD(0.09, int(1953125/5))
    psds[i] = getResult(10, ps5000a)
    print(i)

steps = [1]*23
psds = [[]]*len(steps)
t0 = time.time()
for step,i in zip(steps, range(len(steps))):
    axis_z.MoveBy(step)
    ps5000a.configurePSD(10, 1953125)
    CH.QPD_T_balance(ps5000a, UC8)
    ps5000a.configurePSD(0.09, 1953125)
    psds[i] = getResult(10, ps5000a)
    print(i)
    print(time.time() - t0)

for i,z in zip(range(len(psds)), zpositions):
    psds[i].tofile(wd+'\\'+'psds-z={z:.2f}.bin'.format(z = z))

rm = visa.ResourceManager('C:\\Program Files (x86)\\IVI Foundation\\VISA\\WinNT\\agvisa\\agbin\\visa32.dll')
E4436B = rm.open_resource('GPIB1::19::INSTR')

t0 = time.time()
test = ms.Single_sided_PSD(ps5000a.chs['A'].timeSignal[0:39063235], ps5000a.getfs())
print(time.time()-t0)

amp =
dc =
FG.Sine(0.01, amp, offset = dc)

dcs = np.linspace(0, 0.2, 21)
dcs = dcs.tolist()
zs = [[]]*len(dcs)
for dc, i in zip(dcs, range(len(dcs))):
    FG.DC(dc, ch = 2)
    temp, HAs, HBs, HCs = CH.z_search(axis_z, ps5000a, UC8)
    zs[i] = temp[-1]

t0 = time.time()
monitorTime = 12*3600
minz = [[]]
timestamps = [[]]
step = 0.5
start = 0.1
end = 8
while time.time() - t0 < monitorTime:
    while np.absolute(axis_z.Position() - start) > 0.5:
        axis_z.MoveBy(-np.sign(axis_z.Position() - start)*0.5)
    zs, HAs, HCs = CH.z_sweep(axis_z, ps5000a, UC8, start, end, step)
    popt = fitparaballic(zs, HCs)
    a,b,c = popt
    minz = minz + [-b/(2*a)]
    timestamps = timestamps + [time.time()]
    start = max(0.1, minz[-1] - 4)
    end = min(24.8, minz[-1] + 4)

del minz[0]


def fitparaballic(x, y):
    popt, pcov = curve_fit(parabalic, np.array(x), np.array(y))
    return popt

def parabalic(x, a, b, c):
    return a*x**2 + b*x +c

HAs = []
HBs = []
HCs = []
for i in range (8):
    UC8.PR(axis_y, 50, err)
    HA, HB, HC = CH.max_in_psd(20, ps5000a)
    HAs = HAs + [HA]
    HBs = HBs + [HB]
    HCs = HCs + [HC]

ps5000a.AC('B')
ps5000a.DC('A')
ps5000a.DC('C')
ps5000a.ChangeRangeTo('B', 2000)
ps5000a.ChangeRangeTo('A', 500)
ps5000a.ChangeRangeTo('C', 500)
RFf = np.linspace(79, 81, 11)
RFp = np.linspace(7, 10, 6)
mBs = []
mCs = []
for f in RFf.tolist():
#for p in RFp.tolist():
    E4436B.write(':FREQ:CW {f:.2f}MHz'.format(f = f))
    #E4436B.write(':POW:AMPL {p:.2f}dBm'.format(p = p))
    #CH.QPD_T_balance(ps5000a, UC8)
    HA, HB, HC = CH.max_in_psd(5, ps5000a, 216800)
    mBs = mBs + [HB]
    mCs = mCs + [HC]

xss = np.linspace(-0.2, 0.2, 11)
yss = np.linspace(-0.2. 0.2, 11)
for x in xss:
    for y in yss:
        CH.QPD_T_balance(ps5000a, UC8, x, y)
        HA, HB, HC = max_in_psd(10, ps5000a, 216000)

ps5000a = ps.picoscope()
ps5000a.defaultSetting()
ps5000a.AC('A')
ps5000a.DC('B')
ps5000a.AC('C')
ps5000a.DC('D')
ps5000a.ChangeRangeTo('A', 500)
ps5000a.ChangeRangeTo('B', 10000)
ps5000a.ChangeRangeTo('C', 500)
ps5000a.ChangeRangeTo('D', 5000)
ps5000a.configurePSD(10, int(1953125/5))
ps5000a.LPF20M('C', 1)

temp = np.zeros(100)
for i in range(100):
    temp[i] = np.argmax(psds[i]) - np.argmin(psds[i])

for i in range(len(psds)):
    plt.plot(psds[i])
    plt.show()


ps5000a.DC('A')
ps5000a.DC('C')
ps5000a.configurePSD(3, int(1953125/5))
CH.QPD_T_balance(ps5000a, UC8)
ps5000a.AC('A')
ps5000a.AC('C')
ps5000a.configurePSD(0.099, int(1953125/5))
psd = getResult(1, ps5000a)
plt.plot(psd)
plt.show()


t0 = time.time()
monitorTime = 600
minz = [[]]
timestamps = [[]]
HCs_all = [[]]
while time.time() - t0 < monitorTime:
    zs, HAs, HBs, HCs, xs, ys = CH.z_sweep(axis_z, ps5000a, UC8, start, end, step)
    popt = fitparaballic(zs, HCs)
    a,b,c = popt
    minz = minz + [-b/(2*a)]
    timestamps = timestamps + [time.time()]
    HCs_all = HCs_all + [HCs]

t0 = time.time()
while time.time() - t0 < monitorTime:
    zs, HAs, HBs, HCs, xs, ys = CH.z_sweep(axis_z, ps5000a, UC8, start, end, step, 0, 0)

for i in range(len(HCs_all)-1):
    plt.plot(zs - minz[i+1], HCs_all[i+1])
    plt.pause(0.05)

for i in range(10):
    HA, HB, HC, fftC = CH.max_in_fft(1, ps5000a, 216000)
    plt.plot(np.angle(fftC, deg = True))
    plt.pause(0.1)

for i in range(5):
    HA, HB, HC, fftB, fftC = CH.max_in_fft(10, ps5000a, 216000, 'B')
    fftB = np.array(fftB)
    fftC = np.array(fftC)
    #plt.plot(np.angle(np.divide(fftC,fftB), deg = True))
    plt.plot(np.angle(fftC, deg = True))
    plt.pause(0.1)

axis_z.Quasi_MoveTo(0.1)
steps = [2]*10
steps = [0.9, 4, 2, 1, 1, 1, 1, 1, 3, 4]
zs = np.zeros(len(steps))
complexA = np.zeros(len(steps), dtype = np.complex)
complexB = np.zeros(len(steps), dtype = np.complex)
complexC = np.zeros(len(steps), dtype = np.complex)
for step,i in zip(steps, range(len(steps))):
    axis_z.MoveBy(step)
    ps5000a.DC('A')
    ps5000a.DC('C')
    ps5000a.configurePSD(3, int(1953125/5))
    CH.QPD_T_balance(ps5000a, UC8, debug = False)
    ps5000a.AC('A')
    ps5000a.AC('C')
    ps5000a.configurePSD(1, int(1953125/5))
    complexA[i], complexB[i], complexC[i]= CH.max_in_fft_debug(20, ps5000a, 222000, 'B')
    zs[i] = axis_z.Position()

amp = np.log10(np.square(np.absolute(fftC)))
amp = (amp - np.max(amp)/2 - np.min(amp)/2)/(np.max(amp) - np.min(amp))
angle = np.angle(fftC)

amp = np.sqrt(np.square(A*(z-z0)+B*np.sin(phi1))+B**2*np.sin(phi1)**2)
phi = phi0 + np.arccos(B*np.sin(-phi1)/amp)

def cost(x, *arg):
    A, B, z0, phi0, phi1 = x
    z, amp, angle = arg
    amp_n = (amp - (np.max(amp) + np.min(amp))/2)/((np.max(amp) - np.min(amp))/2)
    angle_n = (angle - (np.max(angle) + np.min(angle))/2)/((np.max(angle) - np.min(angle))/2)
    fit_amp, fit_angle = fit_amp_angle(z, A, B, z0, phi0, phi1)
    fit_amp_n = (fit_amp - (np.max(amp) + np.min(amp))/2)/((np.max(amp) - np.min(amp))/2)
    fit_angle_n = (fit_angle - (np.max(angle) + np.min(angle))/2)/((np.max(angle) - np.min(angle))/2)
    return (np.sum(np.square(fit_amp_n - amp_n)) + np.sum(np.square(fit_angle_n - angle_n)))/len(z)

def fit_amp_angle(z, A, B, z0, phi0, phi1):
    amp = np.sqrt(np.square(-A*(z-z0)+B*np.cos(phi1))+np.square(B*np.sin(phi1)))
    phi = phi0 + np.arccos(np.divide(-A*(z-z0)+B*np.cos(phi1),amp))
    return amp, phi

ampsC = lin_p
anglesC = np.angle(complexC0s/complexC0s[0])
bnd = ((-0.01, 0.01), (-0.01, 0.01), (0, 25), (-6.283185307179586, 6.283185307179586), (-6.283185307179586, 6.283185307179586))
x0 = [5e-04, np.min(np.absolute(complexC)), 15, 0, 0.2]
res = CH.minimize(CH.cost, x0, args = (zs, ampsC, anglesC), bounds = bnd, method = 'Powell', tol = 1e-7)

test = np.zeros(len(ampsC), dtype = np.complex)
test = np.multiply(ampsC, np.cos(anglesC)+1j* np.sin(anglesC))
residule = B*np.complex(np.cos(phi0 + phi1), np.sin(phi0+phi1))
plt.plot(np.angle(test - residule))
plt.plot(np.absolute(test - residule))

def param_manual(z, complex, A, B, z0, phi0, phi1):
    amp, phi = fit_amp_angle(z, A, B, z0, phi0, phi1)
    plt.plot(z, np.absolute(complex))
    plt.plot(z, amp)
    plt.show()
    plt.plot(z, np.angle(complex/complex[0]))
    plt.plot(z, phi)
    plt.show()
    return [A, B, z0, phi0, phi1]

monitorTime = 12*3600
t0 = time.time()
res = []
while time.time() - t0< monitorTime:
    axis_z.Quasi_MoveTo(0.1)
    steps = [1]*20
    zs = np.zeros(len(steps))
    complexA = np.zeros(len(steps), dtype = np.complex)
    complexB = np.zeros(len(steps), dtype = np.complex)
    complexC = np.zeros(len(steps), dtype = np.complex)
    for step,i in zip(steps, range(len(steps))):
        axis_z.MoveBy(step)
        ps5000a.DC('A')
        ps5000a.DC('C')
        ps5000a.configurePSD(3, int(1953125/5))
        CH.QPD_T_balance(ps5000a, UC8)
        ps5000a.AC('A')
        ps5000a.AC('C')
        ps5000a.configurePSD(1, int(1953125/5))
        complexA[i], complexB[i], complexC[i]= CH.max_in_fft_debug(20, ps5000a, 216000, 'B')
        zs[i] = axis_z.Position()
    x0 = [5e-04, np.min(np.absolute(complexC)), 10, 0, 0.2]
    res = res + [minimize(cost, x0, args = (zs,  np.absolute(complexC), np.angle(complexC/complexC[0])), bounds = bnd, method = 'Powell', tol = 1e-7)]

zs = [0.1, 1, 5, 7, 8, 8.8, 8.9, 9, 9.1, 9.2, 10, 11, 12, 15, 20]
steps = np.array(zs[1:]) - np.array(zs[0:-1])

monitorTime = 1
t0 = time.time()
dicts = []
while time.time() - t0< monitorTime:
    axis_z.Quasi_MoveTo(0.1)
    zs = np.zeros(len(steps))
    complexA = np.zeros(len(steps), dtype = np.complex)
    complexB = np.zeros(len(steps), dtype = np.complex)
    complexC = np.zeros(len(steps), dtype = np.complex)
    for step,i in zip(steps, range(len(steps))):
        axis_z.MoveBy(step)
        ps5000a.DC('A')
        ps5000a.DC('C')
        ps5000a.configurePSD(3, int(1953125/5))
        CH.QPD_T_balance(ps5000a, UC8, debug = False)
        ps5000a.AC('A')
        ps5000a.AC('C')
        complexA[i], complexB[i], complexC[i]= CH.max_in_fft_debug(20, ps5000a, 216000, 'B')
        if step < 0.2:
            ps5000a.configurePSD(0.09, int(1953125/5))
            psds = psds + [np.array(getResult(10, ps5000a))]
            psds_z = psds_z + [axis_z.Position()]
        zs[i] = axis_z.Position()
    x0 = [5e-04, np.min(np.absolute(complexC)), 10, 0, 0.2]
    res = minimize(cost, x0, args = (zs,  np.absolute(complexC), np.angle(complexC/complexC[0])), bounds = bnd, method = 'Powell', tol = 1e-7)
    dicts = dicts + [{'fit':res, 'z':zs, 'complexA':complexA, 'complexC':complexC, 'psds':psds, 'psds_z':}]

def getResult_fft(avg, ps5000a, fc):
    #only normalize to channel B's phase, ie, fftB is pure real
    nT = int(39063230/5)
    ps5000a.configurePSD(0.099/avg, int(1953125/5))
    ps5000a.getTimeSignal(trigger = 'D', triggerThreshold = 8000)
    TSidx = [0]
    for i in range(len(ps5000a.chs['D'].timeSignal) - 2):
        if ps5000a.chs['D'].timeSignal[i+1] - ps5000a.chs['D'].timeSignal[i] > 2:
            TSidx = TSidx + [i+1]
    print(TSidx)
    idx = [0, 0]
    idx[0] = np.argmin(np.absolute(ps5000a.f - fc +100))
    idx[1] = np.argmin(np.absolute(ps5000a.f - fc -100))
    for i in range(avg):
        #fftB = fft_process(ps5000a.chs['B'].timeSignal[i*nT:(i+1)*nT], idx)
        #fftC = fft_process(ps5000a.chs['C'].timeSignal[i*nT:(i+1)*nT], idx)
        fftB = fft_process(ps5000a.chs['B'].timeSignal[TSidx[i]:TSidx[i+1]], idx)
        fftC = fft_process(ps5000a.chs['C'].timeSignal[TSidx[i]:TSidx[i+1]], idx)
        norm = np.zeros(len(fftB), dtype = np.complex)
        norm = np.cos(np.angle(fftB)) + 1j*np.sin(np.angle(fftB))
        if i == 0:
            ref = np.divide(fftB, norm)
            signal = np.divide(fftC, norm)
        else:
            ref = ref + np.divide(fftB, norm)
            signal = signal + np.divide(fftC, norm)
    signal = signal/avg
    ref = ref/avg
    return signal, ref, idx

def fft_process(TS, idx):
    fftch = ms.Single_sided_fft(TS)/np.sqrt(len(TS)*ps5000a.getfs())
    fftch = fftch[idx[0]:idx[1]]
    return fftch

zs0 = [0.1, 5, 7, 8, 9, 10, 11, 12, 15, 20]
monitorTime = 12*3600
t0 = time.time()
dicts = []
dict_to_save = {}
iter = 0
while time.time() - t0< monitorTime:
    axis_z.Quasi_MoveTo(0.1)
    if iter == 0:
        minz = 9
    else:
        minz = dict_to_save['fit'].x[2]
    sweep_range = ((np.random.rand(1)-0.5)*2*0.16 + 0.24)[0]
    target_zs = np.around(np.linspace(minz - sweep_range, minz + sweep_range, 5), decimals = 2)
    target_zs = target_zs.tolist()
    zs = zs0 + target_zs
    zs.sort()
    steps = np.array(zs[1:]) - np.array(zs[0:-1])
    zs = np.zeros(len(steps))
    fftCs = []
    fftBs = []
    for step,i in zip(steps, range(len(steps))):
        axis_z.MoveBy(step)
        ps5000a.DC('A')
        ps5000a.DC('C')
        ps5000a.configurePSD(3, int(1953125/5))
        CH.QPD_T_balance(ps5000a, UC8, debug = False)
        ps5000a.AC('A')
        ps5000a.AC('C')
        if step > 0.2:
            fftC, fftB = getResult_fft(1, ps5000a, fc)
        else:
            fftC, fftB = getResult_fft(5, ps5000a, fc)
        fftCs = fftCs + [fftC]
        fftBs = fftBs + [fftB]
        zs[i] = axis_z.Position()
    dict_to_save = CH.process_fft(zs, fftBs, fftCs)
    if iter == 0:
        i = 0
    else:
        i = len(os.listdir(wd))
    dict_to_json(dict_to_save, wd, 'result_{:d}.json'.format(i))
    print(time.time() - t0)
    iter = iter + 1

zs = [0.1, 5, 7, 8, 8.8, 8.9, 9, 9.1, 9.2, 9.3, 9.4, 10, 11, 12, 15, 20]
steps = np.array(zs[1:]) - np.array(zs[0:-1])
axis_z.Quasi_MoveTo(0.1)
zs = np.zeros(len(steps))
fftCs = []
fftBs = []
for step,i in zip(steps, range(len(steps))):
    axis_z.MoveBy(step)
    ps5000a.DC('A')
    ps5000a.DC('C')
    ps5000a.configurePSD(3, int(1953125/5))
    CH.QPD_T_balance(ps5000a, UC8, debug = False)
    ps5000a.AC('A')
    ps5000a.AC('C')
    if step > 0.2:
        fftC, fftB = getResult_fft(1, ps5000a, fc)
    else:
        fftC, fftB = getResult_fft(5, ps5000a, fc)
    fftCs = fftCs + [fftC]
    fftBs = fftBs + [fftB]
    zs[i] = axis_z.Position()
    print(i)

cpy_fftCs = [[]]*len(fftCs)
for i in range(len(fftCs)):
    cpy_fftCs[i] = fftCs[i].copy()

cpy_fftBs = [[]]*len(fftBs)
for i in range(len(fftBs)):
    cpy_fftBs[i] = fftBs[i].copy()


index = [500, 1500]
for i in range(len(cpy_fftCs)):
    bavg = np.mean(cpy_fftBs[i][index[0]:index[1]])
    cpy_fftCs[i][index[0]:index[1]] = np.divide(cpy_fftCs[i][index[0]:index[1]], cpy_fftBs[i][index[0]:index[1]])*bavg

complexC = np.zeros(len(cpy_fftCs), dtype = np.complex)
for i in range(len(cpy_fftCs)):
    maxI = np.argmax(np.absolute(cpy_fftCs[i]))
    indexarra = np.array(list(range(index[0], maxI - 100)) + list(range(maxI + 100, index[1])))
    temp = cpy_fftCs[i]
    complexC[i] = np.mean(temp[indexarra])

for i in range(len(cpy_fftCs)):
    cpy_fftCs[i] = cpy_fftCs[i]/complexC[0]*np.absolute(complexC[0])

complexC = np.zeros(len(cpy_fftCs), dtype = np.complex)
for i in range(len(cpy_fftCs)):
    maxI = np.argmax(np.absolute(cpy_fftCs[i]))
    indexarra = np.array(list(range(index[0], maxI - 100)) + list(range(maxI + 100, index[1])))
    temp = cpy_fftCs[i]
    complexC[i] = np.mean(temp[indexarra])


res = minimize(cost, x0, args = (zs,  np.absolute(complexC), np.angle(complexC)), bounds = bnd, method = 'Powell', tol = 1e-7)
A, B, z0, phi0, phi1 = res.x
for i in range(len(cpy_fftCs)):
    cpy_fftCs[i] = cpy_fftCs[i] - B*np.complex(np.cos(phi0+phi1), np.sin(phi0+phi1))

for i in range(len(cpy_fftCs)):
    complexC[i] = complexC[i] - B*np.complex(np.cos(phi0+phi1), np.sin(phi0+phi1))

plt.plot(zs, np.angle(complexC - B*np.complex(np.cos(phi0+phi1), np.sin(phi0+phi1))))
wd = 'Z:\\data\\optical lever project\\NORCADA_NX53515C\\166-chirp'
for i in range(len(dicts)):
    #dicts[i]['fit'].direc = np.array(dicts[i]['fit'].direc)
    #dicts[i]['fit'].x = np.array(dicts[i]['fit'].x)
    dict_to_json(dicts[i], wd, 'result_{:d}.json'.format(i))

zs0 = [0.1, 5, 7, 8, 9, 10, 11, 12, 15, 20]

axis_z.Quasi_MoveTo(0.1)
steps = [1]*7
zs = np.zeros(len(steps))
complexA = np.zeros(len(steps), dtype = np.complex)
complexB = np.zeros(len(steps), dtype = np.complex)
complexC = np.zeros(len(steps), dtype = np.complex)
HAs = np.zeros(len(steps))
HBs = np.zeros(len(steps))
HCs = np.zeros(len(steps))
ps5000a.configurePSD(1, int(1953125/5))
for step,i in zip(steps, range(len(steps))):
    axis_z.MoveBy(step)
    ps5000a.DC('C')
    CH.balancer_prism(ps5000a, prism_x, 'C', -1, 0.01, offset = 0, debug = False)
    ps5000a.AC('C')
    HAs[i], HBs[i], HCs[i] = CH.max_in_psd(5, ps5000a, 222000)
    #complexA[i], complexB[i], complexC[i]= CH.max_in_fft_debug(20, ps5000a, 222000, 'B')
    zs[i] = axis_z.Position()

def check_offset(TS, length, idx):
    result = fft_process(TS[:length], idx)
    plt.plot(10*np.log10(np.square(np.absolute(result))))
    plt.show()

def getResult_fft(avg, ps5000a, fcs):
    #only normalize to channel B's phase, ie, fftB is pure real
    nT = int(39063230/5)
    nT = 7812489
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
    all_ffts = [0]*avg
    for i in range(len(fcs)):
        idx[2*i] = np.argmin(np.absolute(ps5000a.f - fcs[i] +100))
        idx[2*i+1] = np.argmin(np.absolute(ps5000a.f - fcs[i] -100))
    for i in range(avg):
        fftB = fft_process(ps5000a.chs['B'].timeSignal[i*nT:(i+1)*nT], idx)
        fftC = fft_process(ps5000a.chs['C'].timeSignal[i*nT:(i+1)*nT], idx)
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
        all_ffts[i] = np.divide(fftC[1], norm)
    for j in range(len(fcs)):
        signal[j] = signal[j]/avg
        ref[j] = ref[j]/avg
    return signal, ref, idx, all_ffts

def fft_process(TS, idx):
    fftch = ms.Single_sided_fft(TS)/np.sqrt(len(TS)*ps5000a.getfs())
    ffts = []
    for i in range(int(len(idx)/2)):
        ffts = ffts +[fftch[idx[2*i]:idx[2*i+1]]]
    return ffts

def substract(TS0, TS1, offset):
    if offset == 0:
        re = TS0 - TS1
    if offset < 0:
        re = TS0[:offset] - TS1[-offset:]
    if offset > 0:
        re = TS0[offset:] - TS1[:-offset]
    return re

def plot_comp(signal, offset):
    if offset == 0:
        plt.plot(10*np.log10(np.square(np.absolute(signal[0]))))
        plt.plot(10*np.log10(np.square(np.absolute(signal[1]))))
    if offset < 0:
        plt.plot(10*np.log10(np.square(np.absolute(signal[0][:offset]))))
        plt.plot(10*np.log10(np.square(np.absolute(signal[1][-offset:]))))
    if offset > 0:
        plt.plot(10*np.log10(np.square(np.absolute(signal[0][offset:]))))
        plt.plot(10*np.log10(np.square(np.absolute(signal[1][:-offset]))))


signal, ref, idx = getResult_fft(1, ps5000a, fcs)
plt.plot(substract(10*np.log10(np.square(np.absolute(signal[0]))), 10*np.log10(np.square(np.absolute(signal[1]))), 2))
plt.plot(10*np.log10(np.square(np.absolute(signal[0]))))
result = FG.chirps('test2', [200000, 200500], [200100, 200600], int(1e6), int(10e6), [0.2]*2, [-0.2]*2)
fcs = [200050, 200550]
plt.plot(substract(10*np.log10(np.square(np.absolute(signal[0]))), 10*np.log10(np.square(np.absolute(signal[1]))), 0))
plt.show()


zs = np.linspace(12, 14.5, 11)
zs0 = [1, 5, 10, 15, 20]
steps = np.array(zs[1:]) - np.array(zs[0:-1])
axis_z.Quasi_MoveTo(zs[0])
zs = np.zeros(len(steps))
fftC0s = []
fftC1s = []
fftB0s = []
fftB1s = []
for step,i in zip(steps, range(len(steps))):
    axis_z.MoveBy(step)
    ps5000a.configurePSD(1, int(1953125/5))
    ps5000a.ChangeRangeTo('A', 1000)
    CH.balancer_prism(ps5000a, prism_x, 'A', -1, 0.005, offset = 0, debug = False, rel = 1.5)
    ps5000a.ChangeRangeTo('A', 2000)
    ps5000a.configurePSD(0.099, int(1953125/5))
    if step >= 0:
        signal, ref, idx, fftsall = getResult_fft(1, ps5000a, fcs)
    else:
        signal, ref, idx, fftsall  = getResult_fft(5, ps5000a, fcs)
    fftC0s = fftC0s + [signal[0]]
    fftC1s = fftC1s + [signal[1]]
    fftB0s = fftB0s + [ref[0]]
    fftB1s = fftB1s + [ref[1]]
    zs[i] = axis_z.Position()
    print(i)


index = [500, 1500]
complexC0s = np.zeros(len(zs), dtype = np.complex)
for i in range(len(zs)):
    complexC0s[i] = np.mean(fftC0s[i][index[0]:index[1]])

complexC1s = np.zeros(len(zs), dtype = np.complex)
for i in range(len(fftC1s)):
    maxI = np.argmax(np.absolute(fftC1s[i]))
    indexarra = np.array(list(range(index[0], maxI - 100)) + list(range(maxI + 100, index[1])))
    temp = fftC1s[i]
    complexC1s[i] = np.mean(temp[indexarra])

power = 10*np.log10(20)+10*np.log10(np.square(np.absolute(complexC0s)))
phase = np.angle(complexC0s)
lin_p = np.absolute(complexC0s)

power = 10*np.log10(20)+10*np.log10(np.square(np.absolute(complexC1s)))
phase = np.angle(complexC1s)
lin_p = np.absolute(complexC1s)

signal, ref, idx = getResult_fft(1, ps5000a, fcs)
power = 10*np.log10(20) + 10*np.log10(np.mean(np.square(np.absolute(signal[0][index[0]:index[1]]))))

RFf = np.linspace(79, 81, 11)
RFp = np.linspace(7, 10, 6)
mCs = []
for f in RFf.tolist():
#for p in RFp.tolist():
    E4436B.write(':FREQ:CW {f:.2f}MHz'.format(f = f))
    #E4436B.write(':POW:AMPL {p:.2f}dBm'.format(p = p))
    #CH.QPD_T_balance(ps5000a, UC8)
    signal, ref, idx, _ = getResult_fft(1, ps5000a, fcs)
    temp = 10*np.log10(20) + 10*np.log10(np.mean(np.square(np.absolute(signal[1][index[0]:index[1]]))))
    mCs = mCs + [temp]

axis_z.MoveBy(0.2)
ps5000a.DC('C')
ps5000a.configurePSD(1, int(1953125/5))
CH.balancer_prism(ps5000a, prism_x, 'C', -1, 0.01, offset = 0, debug = False)
ps5000a.AC('C')
ps5000a.configurePSD(0.099, int(1953125/5))
signal, ref, idx = getResult_fft(1, ps5000a, fcs)
power = 10*np.log10(20) + 10*np.log10(np.mean(np.square(np.absolute(signal[0][index[0]:index[1]]))))

ampsC = lin_p
anglesC = np.angle(complexC0s/complexC0s[0])
bnd = ((-0.01, 0.01), (-0.01, 0.01), (0, 25), (-6.283185307179586, 6.283185307179586), (-6.283185307179586, 6.283185307179586))
x0 = [5e-04, np.min(np.absolute(complexC0s)), 15, 0, 0.2]
res = CH.minimize(CH.cost_complex_test, x0, args = (zs, ampsC, anglesC), bounds = bnd, method = 'Powell', tol = 1e-7)

z = np.linspace(zs[0], zs[-1], 1000)
amp, phi = CH.fit_amp_angle(z, *res.x)
fig, ax1 = plt.subplots()
ax1.plot(zs, ampsC, '*')
ax1.plot(z, amp)
ax2 = ax1.twinx()
ax2.plot(zs, anglesC, '*')
ax2.plot(z, phi)
plt.show()

A, B, z0, phi0, phi1 = res.x
test = np.zeros(len(ampsC), dtype = np.complex)
test = np.multiply(ampsC, np.cos(anglesC)+1j* np.sin(anglesC))
residule = B*np.complex(np.cos(phi0 + phi1), np.sin(phi0+phi1))
plt.plot(np.angle(test - residule))
plt.plot(np.real((test - residule)/(test[0] - residule)))

fcs = [fc]
fstarts = (np.array(fcs) - 150).tolist()
fends = (np.array(fcs) + 150).tolist()
fs = int(1e6)
N = int(1e7)
Vmaxs = [1]*len(fcs)
Vmins = [-1]*len(fcs)

zs0 = [0.1, 5, 10, 12, 13, 14, 15, 18]
zs1 = [4.99999078,9.99997919,11.9999781,12.9999687,13.9999593,14.5998404,14.6998236,14.799789,14.8997721,14.999773,14.9998614,16.4398592,16.4998562,16.5598531,16.6198501,16.6798471,17.9998332]
dicts = []
dict_to_save = {}
t0 = time.time()
iter = 1
monitorTime = 8*3600
while time.time() - t0< monitorTime:
    if iter == 0:
        minz = 11
    else:
        minz = dict_to_save['fit'].x[2]
    sweep_range = ((np.random.rand(1)-0.5)*2*0.16 + 0.24)[0]
    target_zs = np.around(np.linspace(minz - sweep_range, minz + sweep_range, 8), decimals = 2)
    target_zs = target_zs.tolist()
    zs = zs0 + target_zs
    if iter > 0:
        minz = dict_to_save['zs'][np.argmin(np.absolute(dict_to_save['complexC']))]
        target_zs = (np.linspace(minz - 0.2, minz + 0.2, 5)).tolist()
        zs = zs + target_zs
    zs.sort()
    zs = np.array(zs)
    zs_index = np.argwhere(np.logical_and(zs>0, zs<25))
    zs = zs[zs_index.flatten()]
    zs = zs.tolist()
    zs = zs1
    axis_z.Quasi_MoveTo(zs[0])
    steps = np.array(zs[1:]) - np.array(zs[0:-1])
    zs = np.zeros(len(steps))
    signal0s = []
    signal1s = []
    idxs = []
    maxfs = []
    for step,i in zip(steps, range(len(steps))):
        axis_z.MoveBy(step)
        ps5000a.configurePSD(1, int(1953125/5))
        ps5000a.ChangeRangeTo('A', 500)
        CH.balancer_prism(ps5000a, prism_x, 'A', -1, 0.005, offset = 0, debug = False)
        ps5000a.ChangeRangeTo('A', 5000)
        if step > 0.32:
            signal, ref, idx, fftsall  = getResult_fft(1, ps5000a, fcs)
        else:
            signal, ref, idx, fftsall  = getResult_fft(5, ps5000a, fcs)
        tempf = ps5000a.f[idx[2]:idx[3]]
        maxf = tempf[np.argmax(np.absolute(signal[1]))]
        print(maxf)
        maxfs = maxfs + [maxf]
        if np.absolute(maxf - fcs[1]) > 20:
            fcs[1] = int(maxf)
            fstarts = (np.array(fcs) - 50).tolist()
            fends = (np.array(fcs) + 50).tolist()
            arb_written = FG.chirps('test', fstarts, fends, fs, N, Vmaxs, Vmins)
        signal0s = signal0s + [signal[0]]
        signal1s = signal1s + [signal[1]]
        idxs = idxs + [idx.tolist()]
        zs[i] = axis_z.Position()
        print(zs[i])
    dict_to_save = CH.process_signal_fft(zs, signal0s, signal1s, idxs)
    dict_to_save['maxfs'] = maxfs
    fileIndex = len(os.listdir(wd))
    dict_to_json(dict_to_save, wd, 'result_{:d}.json'.format(fileIndex))
    print(time.time() - t0)
    iter = iter + 1

files = os.listdir(wd)
for file in range(files):
    with open(wd + r'\\' + file) as f
    json_data = json.load(f)
for i in range(len())

maxfs = []
times = []
t0 = time.time()
while time.time() - t0 < 3600:
    signal, ref, idx = getResult_fft(1, ps5000a, fcs)
    tempf = ps5000a.f[idx[2]:idx[3]]
    maxf = tempf[np.argmax(np.absolute(signal[1]))]
    print(maxf)
    maxfs = maxfs + [maxf]
    times = times + [time.time() - t0]
    if np.absolute(maxf - fcs[1]) > 20:
        fcs[1] = int(maxf)
        fstarts = (np.array(fcs) - 50).tolist()
        fends = (np.array(fcs) + 50).tolist()
        arb_written = FG.chirps('test', fstarts, fends, fs, N, Vmaxs, Vmins)

signal, ref, idx = getResult_fft(1, ps5000a, fcs)
tempf = ps5000a.f[idx[2]:idx[3]]
maxf = tempf[np.argmax(np.absolute(signal[1]))]
print(maxf)

7812489 - 7891415
wd = r'Z:\data\optical lever project\NORCADA_NX53515C\169-chirp'

E4436B.write(':FREQ:CW {f:.2f}MHz'.format(f = 80.58))
signal, ref, idx = getResult_fft(1, ps5000a, fcs)
10*np.log10(20) + 10*np.log10(np.mean(np.square(np.absolute(signal[0][index[0]:index[1]]))))

avg = 10
binwidth = 0.25
ps5000a.configurePSD(binwidth, int(1953125/5))
psd025 = ps5000a.getPSD('C', avg = avg, offset = -1)
binwidth = 1
ps5000a.configurePSD(binwidth, int(1953125/5))
psd1 = ps5000a.getPSD('C', avg = avg, offset = -1)
binwidth = 4
ps5000a.configurePSD(binwidth, int(1953125/5))
psd4 = ps5000a.getPSD('C', avg = avg, offset = -1)

plt.plot(f025, psd025)
plt.plot(f1, psd1)
plt.plot(f4, psd4[:-1])
plt.show()

-37.8
-44.8

-83.1
-87.1

lowcut = 218400
highcut = 218700
fs = ps5000a.getfs()
T = 1
nsamples = int(T * fs)
t = np.linspace(0, T, nsamples, endpoint=False)
x = ps5000a.chs['B'].timeSignal
plt.figure(2)
plt.clf()
plt.plot(t, x, label='Noisy signal')

y = ms.butter_bandpass_filter(x, lowcut, highcut, fs, order=3)
plt.plot(t, y, label='Filtered signal (%g Hz)' % f0)
plt.xlabel('time (seconds)')
plt.hlines([-a, a], 0, T, linestyles='--')
plt.grid(True)
plt.axis('tight')
plt.legend(loc='upper left')

plt.show()



ps5000a = ps.picoscope()
ps5000a.defaultSetting()
ps5000a.DC('A')
ps5000a.AC('B')
ps5000a.AC('C')
ps5000a.DC('D')
ps5000a.ChangeRangeTo('A', 500)
ps5000a.ChangeRangeTo('B', 1000)
ps5000a.ChangeRangeTo('C', 1000)
ps5000a.ChangeRangeTo('D', 5000)
ps5000a.configurePSD(10, int(1953125/5))
ps5000a.LPF20M('A', 1)
ps5000a.LPF20M('B', 1)
ps5000a.LPF20M('C', 1)

rp.tx_txt('ACQ:DEC 32')
rp.tx_txt('ACQ:START')
time.sleep(0.1)
rp.tx_txt('ACQ:STOP')

rp.tx_txt('ACQ:SOUR1:DATA?')
buff_string = rp.rx_txt()
buff_string = buff_string.strip('{}\n\r').replace("  ", "").split(',')
ch1 = list(map(float, buff_string))

rp.tx_txt('ACQ:SOUR2:DATA?')
buff_string = rp.rx_txt()
buff_string = buff_string.strip('{}\n\r').replace("  ", "").split(',')
ch2 = list(map(float, buff_string))

CH.lock_in_amplifier(ch1, ch2, 1953125, 218726)

ps5000a.configurePSD(1, int(1953125/5))
ps5000a.ChangeRangeTo('A', 1000)
CH.balancer_prism(ps5000a, prism_x, 'A', -1, 0.005, offset = 0, debug = False)
ps5000a.ChangeRangeTo('A', 2000)

monitorTime = 11*3600
of = []
ps5000a.configurePSD(0.1, int(1953125/5))
ps5000a.getTimeSignal()
ps5000a.PSDfromTS(ps5000a.chs['C'].timeSignal, ps5000a.getfs())
idx1 = np.argmin(np.absolute(ps5000a.f - fc - 200))
idx0 = np.argmin(np.absolute(ps5000a.f - fc + 200))
fcs = []
t0 = time.time()
while time.time() - t0 < monitorTime:
    psdC = ps5000a.getPSD('C', avg = 10, offset = -1)
    if ps5000a.chs['A'].overflow:
        of = of + [time.time() - t0]
        ps5000a.configurePSD(1, int(1953125/5))
        ps5000a.ChangeRangeTo('A', 500)
        CH.balancer_prism(ps5000a, prism_x, 'A', -1, 0.005, offset = 0, debug = False)
        ps5000a.configurePSD(0.1, int(1953125/5))
        ps5000a.ChangeRangeTo('A', 1000)
    idx = np.argmax(psdC[idx0:idx1])
    fcs = fcs + [ps5000a.f[idx + idx0]]
    print(time.time() - t0)

amp, angle = CH.lock_in_amplifier(ps5000a.chs['A'].timeSignal, \
    ps5000a.chs['B'].timeSignal, ps5000a.getfs(), fc)

218607
269203
def tspsd(ts, n):
    m = int(len(ts)/n)
    for i in range(m):
        temp = ts[i*m : (i+1)*m]
        psdet = ms.Single_sided_PSD(temp, 10)
        plt.clf()
        plt.plot(np.linspace(0, 5, len(psdet)), 10*np.log10(psdet))
        plt.ylim([-80, 30])
        plt.pause(0.05)

sweep_back_forward(prism_x, axis_z, ps5000a, axis_z.Position(), axis_z.Position() + 0.5, 0.1, 2)


for i in range(len(data['zs'])):
     plt.plot(ps5000a.f[data['idxs'][i][0]:data['idxs'][i][1]], \
        10*np.log10(np.square(np.absolute(np.array(data['fftCs_Re'][i]) + 1j*np.array(data['fftCs_Im'][i]))))\
         - 10*np.log10(np.square(np.absolute(np.array(data['fftBs_Re'][i]) + 1j*np.array(data['fftBs_Im'][i])))))
     plt.show()

0.5 - kp*et[-1] - ki*np.trapz(et, t) - kd*(et[-1] - et[-2])/(t[-1] - t[-2])
def PID_check(t, et, kp, ki, kd, stop):
    Pt = np.zeros(stop)
    It = np.zeros(stop)
    Dt = np.zeros(stop)
    for i in range(stop-1):
        j = i + 1
        Pt[j] = - kp*et[j]
        It[j] = - ki*np.trapz(et[:j], t[:j])
        Dt[j] = - kd*(et[j] - et[j-1])/(t[j] - t[j-2])
    return Pt, It, Dt

os.chdir(wd)
filenames = sorted(filter(os.path.isfile, os.listdir('.')), key=os.path.getmtime)
datas = [[]]*len(filenames)
for i in range(len(filenames)):
    datas[i] = CH.json_to_dict_172(wd, filenames[i])

signals = [[]]
xs = []
for i in range(len(datas)):
    for j in range(len(datas[i]['processed_C'])):
        if datas[i]['fordata'][j]:
            signals = signals + [datas[i]['processed_C'][j]]
            xs = xs + [datas[i]['zs'][j] - datas[i]['fit']['z0']]

for i in range(len(datas)):
    for j in range():
        signals = signals + [data[i]['processed_C'][j]]
        xs = xs + [data[i]['zs'][j] - data[i]['fit'].x[2]]

del signals[0]
temp = [[]]*len(signals)
xs = np.array(xs)
idxs = xs.argsort()
for i in range(len(idxs)):
    temp[i] = signals[idxs[i]].copy()

del signals
signals = temp.copy()
del temp
xs = xs[idxs]

idx = [500, 1500]
for i in range(len(datas)):
    complexCs = []
    for j in range(len(datas[i]['zs'])):
        ref = datas[i]['fftBs'][j]
        complexCs = complexCs + [np.mean(ref[idx[0]:idx[1]])]
    datas[i]['complexCs_raw'] = complexCs

for i in range(len(datas)):
    complexC = datas[i]['complexCs_raw']
    zs = datas[i]['zs']
    inverted = False
    if zs[0] >zs[-1]:
        zs = np.flip(zs)
        complexC = np.flip(complexC)
        inverted = True
    p = np.polyfit(np.real(complexC), np.imag(complexC), 1)
    b = p[1]
    k = p[0]
    A = ((np.real(complexC[0])+k*np.imag(complexC[0]))/np.sqrt(1+k**2) + \
        (np.real(complexC[-1])+k*np.imag(complexC[-1]))/np.sqrt(1+k**2))/(zs[0] - zs[-1])
    phi0 = np.arctan(k)
    outphase_noise = b*np.cos(phi0)
    datas[i]['fftCs_inphase'] = [[]]*len(datas[i]['fftCs'])
    datas[i]['complexCs_inphase'] = complexC - 1j * outphase_noise * np.exp(1j*phi0)
    datas[i]['relative_x'] = np.zeros(len(datas[i]['fftCs']))
    datas[i]['refs_inphase'] = [[]]*len(datas[i]['fftBs'])
    #plt.plot(np.real(datas[i]['complexCs_raw']), np.imag(datas[i]['complexCs_raw']))
    #plt.plot(np.real(1j * outphase_noise * np.exp(1j*phi0)), np.imag(1j * outphase_noise * np.exp(1j*phi0)), '*')
    for j in range(len(datas[i]['fftCs'])):
        idx1 = int(len(datas[i]['fftBs'][j])/2)-500
        idx2 = int(len(datas[i]['fftBs'][j])/2)+500
        datas[i]['fftCs_inphase'][j] = datas[i]['fftCs'][j]
        datas[i]['refs_inphase'][j] = datas[i]['fftBs'][j]
        datas[i]['fftCs_inphase'][j][idx1:idx2] = datas[i]['fftCs_inphase'][j][idx1:idx2] - 1j * outphase_noise * np.exp(1j*phi0)
        datas[i]['refs_inphase'][j][idx1:idx2] = datas[i]['refs_inphase'][j][idx1:idx2] - 1j * outphase_noise * np.exp(1j*phi0)
        if inverted:
            datas[i]['relative_x'][j] = (np.real(datas[i]['complexCs_inphase'][-1-j])+k*np.imag(datas[i]['complexCs_inphase'][-1-j]))/np.sqrt(1+k**2)
        else:
            datas[i]['relative_x'][j] = (np.real(datas[i]['complexCs_inphase'][j])+k*np.imag(datas[i]['complexCs_inphase'][j]))/np.sqrt(1+k**2)

relative_z_temp = []
signals_temp = []
refs_temp = []
for i in range(len(datas)):
    for j in range(len(datas[i]['zs'])):
        relative_z_temp = relative_z_temp + [datas[i]['relative_x'][j]]
        signals_temp = signals_temp + [datas[i]['fftCs_inphase'][j]]
        refs_temp = refs_temp + [datas[i]['refs_inphase'][j]]
idx = np.argsort(relative_z_temp)
signals = [[]]*len(signals_temp)
refs = [[]]*len(signals_temp)
relative_z = np.zeros(len(relative_z_temp))
for i in range(len(relative_z_temp)):
    signals[i] = signals_temp[idx[i]].copy()
    relative_z[i] = relative_z_temp[idx[i]]
    refs[i] = refs_temp[idx[i]].copy()
del signals_temp
del refs_temp
del relative_z_temp

relative_z = xs
id1 = int(len(signals[0])/2) - 5000
id2 = int(len(signals[0])/2) + 5000
plt.ion()
figure, ax = plt.subplots(figsize=(10, 8))
line1, = ax.plot(10*np.log10(np.square(np.absolute(signals[0][id1:id2]))))
for i in range(len(relative_z)):
    temp = 10*np.log10(20) + 10*np.log10(np.square(np.absolute(signals[i][id1:id2]))) #- 10*np.log10(np.square(np.absolute(refs[i][id1:id2])))
    line1.set_ydata(temp)
    figure.canvas.draw()
    figure.canvas.flush_events()
    #plt.ylim([-20, 20])
    plt.ylim([-110, -30])
    time.sleep(0.05)

for i in range(len(datas)):
    plt.plot(np.real(datas[i]['complexCs_raw']), np.imag(datas[i]['complexCs_raw']))

for i in range(len(datas)):
    plt.plot(np.real(datas[i]['complexCs_inphase']), np.imag(datas[i]['complexCs_inphase']))


for i in range(10):
    complexC = datas[i]['complexCs_raw']
    zs = datas[i]['zs']
    p = np.polyfit(np.real(complexC), np.imag(complexC), 1)
    b = p[1]
    k = p[0]
    #estimates for A need improvement
    A = ((np.real(complexC[0])+k*np.imag(complexC[0]))/np.sqrt(1+k**2) + \
        (np.real(complexC[-1])+k*np.imag(complexC[-1]))/np.sqrt(1+k**2))/(zs[0] - zs[-1])
    #-----
    phi0 = np.arctan(k)
    b/np.cos(phi0)
    fig, ax1 = plt.subplots()
    #plt.plot(np.real(complexC), np.imag(complexC))
    ax1.plot(zs, np.real(complexC), '*')
    ax1.plot(z, np.real(complexC_fit))
    ax2 = ax1.twinx()
    ax2.plot(zs, np.imag(complexC), 'o')
    ax2.plot(z, np.imag(complexC_fit))
    plt.show()

zs = (np.linspace(8, 16, 9)).tolist()
zs0 = [1, 5, 10, 15, 20]
#zs = zs + zs0
zs.sort()
steps = np.array(zs[1:]) - np.array(zs[0:-1])
axis_z.Quasi_MoveTo(zs[0])
zs = np.zeros(len(steps))
signals = []
refs = []
for step,i in zip(steps, range(len(steps))):
    axis_z.MoveBy(step)
    ps5000a.configurePSD(1, int(1953125/5))
    ps5000a.ChangeRangeTo('A', 1000)
    CH.balancer_prism(ps5000a, prism_x, 'A', -1, 0.005, offset = 0, debug = False, rel = 1.5)
    ps5000a.ChangeRangeTo('A', 2000)
    ps5000a.configurePSD(0.099, int(1953125/5))
    if step >= 0:
        signal, ref, idx, fftsall = CH.getResult_fft(1, ps5000a, fcs)
    else:
        signal, ref, idx, fftsall  = CH.getResult_fft(5, ps5000a, fcs)
    signals = signals + [signal]
    refs = refs + [ref]
    zs[i] = axis_z.Position()
    print(i)

for i in range(len(signals)):
    for j in range(len(signals[i])):
        plt.plot(10*np.log10(np.square(np.absolute(signals[i][j]))))
    plt.show()

complexC = np.zeros([len(signals), len(signals[0])], dtype = np.complex)
for i in range(len(signals)):
    for j in range(len(signals[i])):
        complexC[i][j] = np.mean(signals[i][j])
complexC = np.delete(complexC, 3, 1)

complexCs = complexC.copy()
del complexC
ks = np.zeros(5)
bs = np.zeros(5)
As = np.zeros(5)
phi0s = np.zeros(5)
for i in range(5):
    complexC = complexCs[:, i]
    plt.plot(np.real(complexC), np.imag(complexC))
    p = np.polyfit(np.real(complexC), np.imag(complexC), 1)
    bs[i] = p[1]
    ks[i] = p[0]
    #estimates for A need improvement
    As[i] = ((np.real(complexC[0])+k*np.imag(complexC[0]))/np.sqrt(1+k**2) + \
        (np.real(complexC[-1])+k*np.imag(complexC[-1]))/np.sqrt(1+k**2))/(zs[0] - zs[-1])
    #-----
    phi0s[i] = np.arctan(ks[i])
    print(bs[i]*np.cos(phi0s[i]))

for i in range(5):
    complexC = complexCs[:, i] - 1j*bs[i]*np.cos(phi0s[i])*np.exp(1j*phi0s[i])
    plt.plot(np.real(complexC), np.imag(complexC))
    p = np.polyfit(zs, np.real(complexC), 1)
    z0 = -p[1]/p[0]
    print(z0)

def plot_fun(A, B, z0, phi0, phi1):
    z = np.linspace(0, 25, 1000)
    complexC = A*(z-z0)*np.exp(1j*phi0) + B*np.exp(1j*(phi0+phi1))
    plt.plot(np.real(complexC), np.imag(complexC))

lin0 = np.zeros(5)
lin1 = np.zeros(5)
for i in range(5):
    lin0[i] = (np.real(complexCs[0, i]) + ks[i] * np.imag(complexCs[0, i]))/np.sqrt(1+ks[i]**2)
    lin1[i] = (np.real(complexCs[-1, i]) + ks[i] * np.imag(complexCs[-1, i]))/np.sqrt(1+ks[i]**2)
plt.plot(As, lin0)
plt.plot(As, lin1)
plt.show()
p0 = np.polyfit(lin0, As, 1)
print(p0[1])
p1 = np.polyfit(lin1, As, 1)
print(p1[1])

idx = [id1, id2]
for i in range(len(datas)):
    complexCs = []
    for j in range(len(datas[i]['zs'])):
        ref = datas[i]['fftBs'][j]
        complexCs = complexCs + [np.mean(ref[idx[0]:idx[1]])]
    datas[i]['complexCs_raw'] = complexCs
    for i in range(len(datas)):
        complexC = datas[i]['complexCs_raw']
        zs = datas[i]['zs']
        inverted = False
        if zs[0] >zs[-1]:
            zs = np.flip(zs)
            complexC = np.flip(complexC)
            inverted = True
        p = np.polyfit(np.real(complexC), np.imag(complexC), 1)
        b = p[1]
        k = p[0]
        A = ((np.real(complexC[0])+k*np.imag(complexC[0]))/np.sqrt(1+k**2) + \
            (np.real(complexC[-1])+k*np.imag(complexC[-1]))/np.sqrt(1+k**2))/(zs[0] - zs[-1])
        phi0 = np.arctan(k)
        outphase_noise = b*np.cos(phi0)
        datas[i]['fftCs_inphase'] = [[]]*len(datas[i]['fftCs'])
        datas[i]['complexCs_inphase'] = complexC - 1j * outphase_noise * np.exp(1j*phi0)
        datas[i]['relative_x'] = np.zeros(len(datas[i]['fftCs']))
        datas[i]['refs_inphase'] = [[]]*len(datas[i]['fftBs'])
        #plt.plot(np.real(datas[i]['complexCs_raw']), np.imag(datas[i]['complexCs_raw']))
        #plt.plot(np.real(1j * outphase_noise * np.exp(1j*phi0)), np.imag(1j * outphase_noise * np.exp(1j*phi0)), '*')
        for j in range(len(datas[i]['fftCs'])):
            datas[i]['fftCs_inphase'][j] = datas[i]['fftCs'][j] - 1j * outphase_noise * np.exp(1j*phi0)
            datas[i]['refs_inphase'][j] = datas[i]['fftBs'][j] - 1j * outphase_noise * np.exp(1j*phi0)
            if inverted:
                datas[i]['relative_x'][j] = (np.real(datas[i]['complexCs_inphase'][-1-j])+k*np.imag(datas[i]['complexCs_inphase'][-1-j]))/np.sqrt(1+k**2)
            else:
                datas[i]['relative_x'][j] = (np.real(datas[i]['complexCs_inphase'][j])+k*np.imag(datas[i]['complexCs_inphase'][j]))/np.sqrt(1+k**2)

f = np.linspace(-50, 50, 1000))
angle  = relative_z+ np.pi/2
X,Y = np.meshgrid(f,angle)

fig = plt.figure()
ax = plt.axes(projection='3d')
cax = fig.add_axes([0.8, 0.27, 0.05, 0.5])
im = ax.plot_surface(x, y, z,cmap='gnuplot', edgecolor='none')
fig.colorbar(im, cax=cax, orientation = 'vertical')

fig,ax=plt.subplots(1,1)
levels = np.linspace(-27, 43, 36)
cp = ax.contourf(x, y + np.pi/2, z, levels = levels, cmap = 'jet')
cbar = plt.colorbar(cp)
cbar.ax.set_ylabel('power(dBm)')
plt.ylabel('Quadrature(not to scale)')
plt.xlabel('f-fm(Hz)')
plt.yticks([1.5708])

f = np.linspace(500,1500,10000) - 1000
temp = 10*np.log10(20) + 10*np.log10(np.square(np.absolute(signals[400][5000:15000])))
plt.plot(f, temp)
plt.ylabel('Power(dBm)')
plt.xlabel('f-fm(Hz)')

ax = gca()
fig, ax= plt.subplots()
ax.plot(np.real(complexC), np.imag(complexC),'*')
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.xlabel('Re{classical noise}')
plt.ylabel('Im{classical noise}')
ax.axhline(y=0, color='k')
ax.axvline(x=0, color='k')
plt.show()


fc = 120823
ps5000a.configurePSD(10, int(1953125/5))
freqs = np.linspace(fc-5, fc+5, 11)
lowcut = fc-1e3
highcut = fc+1e3
fs = ps5000a.getfs()
angles = []
amps = []
for f in freqs:
    FG.Sine(f, 10e-3, ch = 2)
    time.sleep(3)
    ps5000a.getTimeSignal()
    amp, angle = CH.lock_in_amplifier(ps5000a.chs['A'].timeSignal, \
        ps5000a.chs['B'].timeSignal, ps5000a.getfs(), f)
    angles = angles + [angle]
    amps = amps + [amp]

#angles = (np.array(angles) - angles[0]+1)%(2*np.pi)
fig, ax1 = plt.subplots()
ax1.plot(freqs, angles, 'r.')
ax2 = ax1.twinx()
ax2.plot(freqs, amps, 'b.')
plt.show()

ps5000a.configurePSD(10, int(1953125/5))
t, et, DCs, ps = CH.PID(-2.4, -0.01, -0.001, -0.01, fc, ps5000a, 30)

fc = 120823
FG.Sine(fc,15e-3, ch = 2)
time.sleep(3)
ps5000a.getTimeSignal()
CH.lock_in_amplifier(ps5000a.chs['A'].timeSignal, ps5000a.chs['B'].timeSignal, ps5000a.getfs(), fc)

ps5000a.configurePSD(10, int(1953125/5))
t, et, DCs, ps = CH.PID(0, -0.01, -0.005, -0.005, fc, ps5000a, 3600)
t, et, DCs, ps = CH.PID(0, 0.01, 0.003, 0.003, fc, ps5000a, 30)
ps5000a.configurePSD(1, int(1953125/5))
def report():
    HA, HB, HC = CH.max_in_psd(5, ps5000a, fc)
    to_report = 10*np.log10(20)+10*np.log10(HC)
    print(to_report)

def costfunc(popt, *arg):
    xdata, ydata = arg
    return np.sqrt(np.sum(np.square(ydata - ms.fitFunc(xdata, *popt)))/len(xdata))

bnd = ((0,5), (320890-5, 320890+5), (0, 300), (0, 0.01), (10,300))
specfit = {'1_10':{'x':xdata.tolist(), 'y':ydata.tolist(), 'popt':res.x.tolist(), 'cost':res.fun}}

fc = 320874.7
mode = '1_12'

FG.Sine(fc, 200e-3, ch = 1)
time.sleep(3)
ps5000a.configurePSD(1, int(1953125/5))
CH.balancer_prism(ps5000a, prism_x, 'A', -1, 0.005, offset = 0, debug = False, rel = 2)
ps5000a.configurePSD(0.1, int(1953125/5))
psd = ps5000a.getPSD('C', avg = 5, offset = -1)

f0 = fc - 750
f1 = fc + 750
i0 = int(f0*10)
i1 = int(f1*10)

ydata = psd[i0:i1]
xdata = range(i0, i1)
xdata = np.array(list(xdata))/10

specfit[mode] = {}
specfit[mode]['x'] = xdata.tolist()
specfit[mode]['y'] = ydata.tolist()

for key in specfit:
    t0 = specfit[key]['t0']
    t1 = specfit[key]['t1']
    p0 = specfit[key]['p0']
    p1 = specfit[key]['p1']
    specfit[key]['gamma'] = np.log(10)*(p1-p0)/(10*(t1-t0))

key = '1_12'
i0s = [ms.Approximate_index(specfit[key]['x'], 320520), ms.Approximate_index(specfit[key]['x'], 320691)]
i1s = [ms.Approximate_index(specfit[key]['x'], 320561), ms.Approximate_index(specfit[key]['x'], 320728)]
specfit[key]['x_filtered'], specfit[key]['y_filtered'] = ms.Exclude_data(specfit[key]['x'], specfit[key]['y'], i0s, i1s)

key = '1_2'
specfit[key]['x_filtered'] = specfit[key]['x']
specfit[key]['y_filtered'] = specfit[key]['y']

for key in specfit:
    if key=='1_8':
        continue
    x = specfit[key]['x_filtered']
    y = specfit[key]['y_filtered']
    p0[1] = x[np.argmax(y)]
    p0[0] = -specfit[key]['gamma']
    bnd = ((0, p0[0]*5), (p0[1]-5, p0[1]+5), (0, 300), (0, 1e-6), (10, 300))
    res = minimize(costfunc, p0, args = (x, y), bounds = bnd, method = 'Nelder-Mead', tol = 1e-7, options = {'maxiter':100000, 'maxfev':100000})
    specfit[key]['fitx'] = (res.x).tolist()
    specfit[key]['fitcost'] = res.fun

for key in specfit:
    if key=='1_8':
        continue
    print(key)
    print(specfit[key]['fitcost'])
    print(specfit[key]['fitx'])

key = '1_2'
x = specfit[key]['x_filtered']
y = specfit[key]['y_filtered']
p0[1] = x[np.argmax(y)]
p0[0] = -specfit[key]['gamma']
bnd = ((0, p0[0]*5), (p0[1]-5, p0[1]+5), (0, 300), (0, 1e-6), (10, 300))
res = minimize(costfunc, p0, args = (x, y), bounds = bnd, method = 'Nelder-Mead', tol = 1e-7, options = {'maxiter':100000, 'maxfev':100000})
plt.plot(x, y, '.')
plt.plot(x, ms.fitFunc(x, *res.x))
plt.show()

for key in specfit:
    if key=='1_8':
        continue
    x = specfit[key]['x_filtered']
    y = specfit[key]['y_filtered']
    fitx = specfit[key]['fitx']
    xdense = np.linspace(fitx[1]-1, fitx[1]+1, 1001)
    xfit = np.sort(x + xdense.tolist())
    print(key)
    print(fitx[0])
    plt.plot(x, y, '.')
    plt.plot(xfit, ms.fitFunc(xfit, *fitx))
    plt.show()

'{:f}, {:f}, {:f}, {:f}'.format(np.mean(ps5000a.chs['A'].timeSignal), np.mean(ps5000a.chs['B'].timeSignal), np.mean(ps5000a.chs['C'].timeSignal), np.mean(ps5000a.chs['D'].timeSignal))

ps5000a.configurePSD(1, int(1953125/5))
def try_PID(P, I, D, fc, phase, threshold, rfp, times):
    freqs = np.linspace(fc-20, fc+20, 11)
    lowcut = fc-1e3
    highcut = fc+1e3
    fs = ps5000a.getfs()
    t=0
    et=0
    DCs = 0
    ps= 0
    for f in freqs:
        FG.Sine(f, rfp, ch = 2)
        time.sleep(3)
        ps5000a.getTimeSignal()
        amp, angle = CH.lock_in_amplifier(ps5000a.chs['A'].timeSignal, \
            ps5000a.chs['B'].timeSignal, ps5000a.getfs(), f)
        print(f, amp, angle)
        if amp > threshold:
            t, et, DCs, ps = CH.PID(phase, P, I, D, f, ps5000a, times)
            break
    return t, et, DCs, ps

ps5000a.configurePSD(1, int(1953125/5*3))
fcs = [81803, 123872, 172525, 223504, 275648, 328107, 371770]
fcs = [326445, 376261, 468337, 609922]
span = 800
psds = [psd]*len(fcs)
monitorTime = 3600
t0 = time.time()
freqs = [[]]*len(fcs)
laserPB = []
laserPD = []
t = []
elapsed_time = 0
ps5000a.getTimeSignal()
ps5000a.PSDfromTS(ps5000a.chs['C'].timeSignal, ps5000a.getfs())
idxs = [[]]*len(fcs)
for i in range(len(fcs)):
    start = ms.Approximate_index(ps5000a.f, fcs[i] - span/2)
    end = ms.Approximate_index(ps5000a.f, fcs[i] + span/2)
    idxs[i] = [start, end]


while elapsed_time < monitorTime:
    ps5000a.getTimeSignal()
    if ps5000a.chs['A'].overflow is True:
        CH.balancer_prism(ps5000a, prism_x, 'A', -1, 0.005, offset = 0, debug = 0, rel = 2)
    psdC = ps5000a.getPSD('C', avg = 5, offset = -1)
    for fc, i in zip(fcs, range(len(fcs))):
        freqs[i] = freqs[i] + [fc + np.argmax(psdC[idxs[i][0]:idxs[i][1]])]
        peak = np.max(psdC[int(fc-span/2): int(fc + span/2)])
        if  peak > -40 and peak > np.max(psds[i][idxs[i][0]:idxs[i][1]]):
            psds[i] = psdC
    laserPB = laserPB + [np.mean(ps5000a.chs['B'].timeSignal)]
    laserPD = laserPD + [np.mean(ps5000a.chs['D'].timeSignal)]
    elapsed_time = time.time() - t0
    t = t + [elapsed_time]
    print(elapsed_time)


monitorTime = 900
laserPA = []
laserPB = []
laserPC = []
laserPD = []
t = []
t0 = time.time()
elapsed_time = time.time() - t0
while elapsed_time < monitorTime:
    ps5000a.getTimeSignal()
    laserPA = laserPA + [np.mean(ps5000a.chs['A'].timeSignal)]
    laserPB = laserPB + [np.mean(ps5000a.chs['B'].timeSignal)]
    laserPC = laserPC + [np.mean(ps5000a.chs['C'].timeSignal)]
    laserPD = laserPD + [np.mean(ps5000a.chs['D'].timeSignal)]
    elapsed_time = time.time() - t0
    t = t + [elapsed_time]
    aP = (laserPA[-1] - 1.1240649313864077)/1.1240649313864077*100
    bP = (laserPB[-1] - 3.6695497877200696)/3.6695497877200696*100
    cP = (laserPC[-1] - 0.2503204888527522)/0.2503204888527522*100
    dP = (laserPD[-1] - 1.3799933309198555)/1.3799933309198555*100
    print('t = {:.2f}; [{:.2f}%, {:.2f}%, {:.2f}%, {:.2f}%]'.format(elapsed_time, aP, bP, cP, dP))

monitorTime = 8*3600
laserPC = []
laserPD = []
t = []
t0 = time.time()
elapsed_time = time.time() - t0
while elapsed_time < monitorTime:
    ps5000a.getTimeSignal()
    laserPC = laserPC + [np.mean(ps5000a.chs['C'].timeSignal)]
    laserPD = laserPD + [np.mean(ps5000a.chs['D'].timeSignal)]
    elapsed_time = time.time() - t0
    t = t + [elapsed_time]
    print('t = {:.2f};'.format(elapsed_time))


times = [0]*len(laserPC)
for i in range(len(times)):
    times[i] = datetime.datetime.fromtimestamp(t0 + t[i])

plt.plot(times, (np.array(laserPA) - np.mean(laserPA))/np.mean(laserPA), label = 'A')
plt.plot(times, (np.array(laserPB) - np.mean(laserPB))/np.mean(laserPB), label = 'B')
plt.plot(times, (np.array(laserPC) - np.mean(laserPC))/np.mean(laserPC), label = 'C')
plt.plot(times, (np.array(laserPD) - np.mean(laserPD))/np.mean(laserPD), label = 'D')
plt.legend(loc = 'upper left')
plt.show()

plt.plot((np.array(laserPA) - np.mean(laserPA))/np.mean(laserPA), label = 'A')
plt.plot((np.array(laserPB) - np.mean(laserPB))/np.mean(laserPB), label = 'B')
plt.plot((np.array(laserPC) - np.mean(laserPC))/np.mean(laserPC), label = 'C')
plt.plot((np.array(laserPD) - np.mean(laserPD))/np.mean(laserPD), label = 'D')
plt.legend(loc = 'upper left')
plt.show()

end = 2200
plt.plot((np.array(laserPA[:end]) - np.mean(laserPA[:end]))/np.mean(laserPA[:end]), label = 'A')
plt.plot((np.array(laserPB[:end]) - np.mean(laserPB[:end]))/np.mean(laserPB[:end]), label = 'B')
plt.plot((np.array(laserPC[:end]) - np.mean(laserPC[:end]))/np.mean(laserPC[:end])+0.02, label = 'C')
plt.plot((np.array(laserPD[:end]) - np.mean(laserPD[:end]))/np.mean(laserPD[:end]), label = 'D')
plt.legend(loc = 'upper left')
plt.show()

plt.plot(-(np.array(laserPA) - np.mean(laserPA))/np.mean(laserPA)*25)
plt.plot((np.array(laserPB) - np.mean(laserPB))/np.mean(laserPB)*100)
plt.plot((np.array(laserPC) - np.mean(laserPC))/np.mean(laserPC)*100)
plt.plot((np.array(laserPD) - np.mean(laserPD))/np.mean(laserPD)*100)

times = [0]*len(freqs[0])
for i in range(len(times)):
    times[i] = datetime.datetime.fromtimestamp(t0 + t[i])

for i in range(len(fcs)):
    plt.plot(times, (freqs[i]-np.mean(freqs[i]))/np.mean(freqs[i])*1e5+3*i)
plt.plot(-(np.array(laserP) - np.mean(laserP))/np.mean(laserP)*100)

for i in range(len(fcs)):
    plt.plot((freqs[i]-np.mean(freqs[i]))/np.mean(freqs[i])*1e5+3*i)
plt.plot(-(np.array(laserP) - np.mean(laserP))/np.mean(laserP)*100)

for i in range(len(fcs)):
    plt.plot(psds[i])
    plt.show()

for i in range(len(fcs)):
    psds[i].tofile(wd + 'psd_C_to_1_{:d}_01.bin'.format(i), sep = '')

monitorTime = 4*3600
delay = 0
setPoint = 3.099942
refch = 'B'
sitPoint = 0.7
t0 = time.time()
laserP = []
errs = []
kp = 0
ki = 0.5
kd = 0
laserPA = []
laserPB = []
laserPC = []
laserPD = []
t = []
err = 0
DC = sitPoint
elapsedTime = time.time() - t0
err = sum(ps5000a.getTimeSignal(refch) - setPoint)/len(ps5000a.chs[refch].timeSignal)
while elapsedTime < monitorTime:
    ps5000a.getTimeSignal()
    if elaspedTime > delay:
        err = err + sum(ps5000a.chs[refch].timeSignal - setPoint)/len(ps5000a.chs[refch].timeSignal)
        errs = errs + [err]
        avgsig = np.mean(ps5000a.chs[refch].timeSignal)
        DC = sitPoint + err*ki + (avgsig - setPoint) * kp
        if np.abs(DC - sitPoint)>0.4:
            break
        else:
            FG.DC(DC)
    laserPA = laserPA + [np.mean(ps5000a.chs['A'].timeSignal)]
    laserPB = laserPB + [np.mean(ps5000a.chs['B'].timeSignal)]
    laserPC = laserPC + [np.mean(ps5000a.chs['C'].timeSignal)]
    laserPD = laserPD + [np.mean(ps5000a.chs['D'].timeSignal)]
    laserP = laserP + [np.mean(ps5000a.chs[refch].timeSignal)]
    elapsedTime = time.time() - t0
    t = t + [elapsedTime]
    print(DC, laserP[-1])


 with open('Book1.csv', 'r') as file:
     reader = csv.reader(file, delimiter = ' ')
     for row in reader:
         r = re.findall('([0-9]+:[0-9]+),([0-9]+\.[0-9]+).+,([0-9]+\.[0-9]+)%RH', row[1])
         T = T + [float(r[0][1])]
         H = H + [float(r[0][2])]
         temp = row[0] + ' ' + r[0][0]
         temp = time.mktime(datetime.datetime.strptime(temp, "%m/%d/%Y %H:%M").timetuple())
         t = t + [datetime.datetime.fromtimestamp(temp)]

hbar = 1.0545718e-34
kb = 1.380649e-23
T = 273.15+24.37
f = np.linspace(0, len(psd)/50, len(psd))
C_1 = 10*np.log10(kb*T/(hbar*2*np.pi*f)) - 105
plt.plot(f, C_1)
plt.plot(f, psd)
plt.show()

hbar = 1.0545718e-34
kb = 1.380649e-23
T = 273.15+24.37
f = np.linspace(0, len(sum)/50, len(sum))
C_1 = 10*np.log10(kb*T/(hbar*2*np.pi*f)) -100
plt.plot(f, C_1)
plt.plot(f, sum)
plt.show()

offset = 60
plt.plot(psdA)
plt.plot(psdB + offset)
plt.plot(psdC + 2*offset)
plt.plot(psdD + 3*offset)
plt.show()


take = 20
psds = [[]]*take
for i in range(take):
    ps5000a.configurePSD(1, int(1953125/5*3))
    if ps5000a.chs['A'].overflow is True:
        CH.balancer_prism(ps5000a, prism_x, 'A', -1, 0.005, offset = 0, debug = 0, rel = 2)
    ps5000a.configurePSD(0.1, int(1953125/5*3))
    temp = ps5000a.getPSD('C', avg = 5, offset = -1)
    psds[i] = temp[idx1:idx2]


avg = 50
sum = []
for i in range(avg):
    ps5000a.configurePSD(1, int(1953125/5*3))
    if ps5000a.chs['A'].overflow is True:
        CH.balancer_prism(ps5000a, prism_x, 'A', -1, 0.005, offset = 0, debug = 0, rel = 2)
    ps5000a.configurePSD(0.02, int(1953125/5*3))
    temp = ps5000a.getPSD('C', avg = 1, offset = -1)
    if ps5000a.chs['A'].overflow is True:
        i = i-1
        continue
    else:
        temp = np.power(10, temp/10)
        if i == 0:
            sum = temp
        else:
            sum = sum + temp
        print(i)

avg = 20
sum = []
i = 0
count = 0
while i < avg:
    ps5000a.configurePSD(1, int(1953125/5*3))
    if ps5000a.chs['A'].overflow is True:
        CH.balancer_prism(ps5000a, prism_x, 'A', -1, 0.005, offset = 0, debug = 0, rel = 2)
    ps5000a.configurePSD(0.02, int(1953125/5*3))
    temp = ps5000a.getPSD('C', avg = 1, offset = -1)
    count = count + 1
    if ps5000a.chs['A'].overflow is True:
        continue
    else:
        temp = np.power(10, temp/10)
        if i == 0:
            sum = temp
        else:
            sum = sum + temp
        print(i)
        i = i+1
sum = sum/avg
sum = 10*np.log10(sum)

CH.PID_ver2(-2.4, -0.03, -0.003, -0.004, fc, ps5000a, wd)

ps5000a.configurePSD(0.02, int(1953125/5))
plt.plot(-(ps5000a.chs['A'].timeSignal - np.mean(ps5000a.chs['A'].timeSignal))/np.mean(ps5000a.chs['A'].timeSignal)+0.05)
plt.plot((ps5000a.chs['B'].timeSignal - np.mean(ps5000a.chs['B'].timeSignal))/np.mean(ps5000a.chs['B'].timeSignal)+0.05)
plt.plot((ps5000a.chs['C'].timeSignal - np.mean(ps5000a.chs['C'].timeSignal))/np.mean(ps5000a.chs['C'].timeSignal)+0.1)
plt.plot((ps5000a.chs['D'].timeSignal - np.mean(ps5000a.chs['D'].timeSignal))/np.mean(ps5000a.chs['D'].timeSignal)+0.15)
plt.plot((temp - np.mean(temp))/np.mean(temp )+0.1)
plt.show()

plt.plot(np.linspace(0,1,len(ps5000a.chs['A'].timeSignal)), ps5000a.chs['A'].timeSignal)

fcs = [123763, 130000, 223324, 226000, 275306, 279000, 371437, 374000]
zs0 = [1, 5, 10, 15]
CH.sweep_back_forward(prism_x, axis_z, ps5000a, zs0, 4, wd, minzs, fcs)





N = len(datas_list)
relative_z_data_list = [[]]*N
check_list = [[]]*N
signals_list = [[]]*N
relative_z_list = [[]]*N
p_list = [0]*N
zs_data_list = [[]]*N
zs_list = [[]]*N
check_list = [[]]*N
check_data_list = [[]]*N
zs_data_list = [[]]*N
relative_z_data_list = [[]]*N
minzs = [0]*N

for i in range(N):
    relative_z_data_list[i], signals_list[i], relative_z_list[i], _, check_data_list[i], zs_data_list[i], zs_list[i],check_list[i]  = CH.BKAdata(datas_list[i], [1,2,3,4])
    #ind = np.argmin(relative_z_data_list[i])
    #zs_data_index = np.argwhere(np.absolute(zs_data_list[i] - zs_data_list[ind]<5))
    #zs_data = zs_data_list[i][zs_data_index.flatten()]
    #relative_z_data = relative_z_data_list[i][zs_data_index.flatten()]
    p = np.polyfit(relative_z_list[i], zs_list[i], 1)
    minzs[i] = p[1]
    plt.plot(relative_z_data_list[i], zs_data_list[i], '*')
    plt.plot(relative_z_list[i], relative_z_list[i]*p[0] + p[1])

iN = 0
signals = signals_list[iN]
relative_z_data = relative_z_data_list[iN]
check_data = check_data_list[iN]

plt.ion()
figure, ax = plt.subplots(figsize=(10, 8))
line1, = ax.plot(signals[0])
for i in range(len(relative_z_data)):
    line1.set_ydata(signals[i])
    figure.canvas.draw()
    figure.canvas.flush_events()
    #plt.ylim([-20, 20])
    plt.ylim([-110, -30])
    time.sleep(0.05)

f = np.linspace(-100, 100, 1000)
angle  = relative_z_data + np.pi/2
x,y = np.meshgrid(f,angle)
z = np.array(signals)
for i in range(z.shape[0]):
    z[i,:] = z[i,:] - 10*np.log10(np.square(np.absolute(relative_z_data[i])))- 10*np.log10(20)
    #z[i,:] = z[i,:] - 10*np.log10(np.square(np.absolute(check_data[i])))- 10*np.log10(20)
fig = plt.figure()
ax = plt.axes(projection='3d')
cax = fig.add_axes([0.8, 0.27, 0.05, 0.5])
im = ax.plot_surface(x, y, z,cmap='gnuplot', edgecolor='none')
fig.colorbar(im, cax=cax, orientation = 'vertical')

fig,ax=plt.subplots(1,1)
levels = np.linspace(-20, 30, 40)
cp = ax.contourf(x, y + np.pi/2, z,  levels = levels, cmap = 'jet')
cbar = plt.colorbar(cp)
cbar.ax.set_ylabel('power(dBm)')
plt.ylabel('Quadrature(not to scale)')
plt.xlabel('f-fm(Hz)')
plt.yticks([1.5708])

plt.ion()
figure, ax = plt.subplots(figsize=(10, 8))
line1, = ax.plot(z[0,:])
for i in range(len(relative_z_data)):
    line1.set_ydata(z[i,:])
    figure.canvas.draw()
    figure.canvas.flush_events()
    #plt.ylim([-20, 20])
    plt.ylim([-40, 50])
    time.sleep(0.05)

CH.sweep_back_forward(prism_x, axis_z, ps5000a, zs0, 10, wd, [0,0,0,0], fcs)


fig,ax=plt.subplots(1,1)
levels = np.linspace(-30, 90, 120)
cp = ax.contourf(x, y, z,  levels = levels, cmap = 'jet')
cbar = plt.colorbar(cp)
cbar.ax.set_ylabel('power(dBm)')
plt.ylabel('Quadrature(degree)')
plt.xlabel('f-fm(Hz)')
plt.yticks([89.6, 89.8, 90, 90.2, 90.4])

[123763, 130000, 223324, 226000, 275306, 279000, 371437, 374000]

fcs = [371571, 373300]
fstarts = (np.array(fcs) - 150).tolist()
fends = (np.array(fcs) + 150).tolist()
fs = int(3.2e6)
N = int(16e6)
Vmaxs = [2]*len(fcs)
Vmins = [-2]*len(fcs)
arb_written = FG.chirps('test', fstarts, fends, fs, N, Vmaxs, Vmins)
ps5000a.getTimeSignal()
ps5000a.Plot('B')

p = np.polyfit(relative_z_data, zs_data_list[0], 1)
real_z = relative_z_data*p[0]
real_z = real_z*1.8253/1e3
real_z = real_z/np.pi*180
angle = real_z + 90
x,y = np.meshgrid(f,angle)
