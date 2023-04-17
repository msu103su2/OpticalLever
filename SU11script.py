import lib.ps5000A as ps
import lib.FuncGen33522b as FuncGen
import lib.mathlib_shan as ms
import lib.splitPair as sp
import matplotlib.pyplot as plt
import numpy as np
import time
import sys
import pyvisa as visa
sys.path.append(r'C:\Program Files\Newport\Piezo Motion Control\Newport AG-UC2-UC8 Applet\Bin')
import clr
clr.AddReference('AgilisCmdLib')
from Newport.AgilisCmdLib import *
import lib.chirp_scan_lib as CH
nT = int(39063235/5)
from datetime import datetime
import json
import gc


SN1 = 'IY144/0152'
SN2 = 'GX150/0053'

fc = 259626
ps5000a.configurePSD(10, int(1953125/5))
freqs = np.linspace(fc-2, fc+2, 11)
lowcut = fc-1e3
highcut = fc+1e3
fs = ps5000a.getfs()
angles = []
amps = []
for f in freqs:
    FG.Sine(f, 50e-3, ch = 2)
    time.sleep(3)
    ps5000a.getTimeSignal()
    amp, angle = CH.lock_in_amplifier(ps5000a.chs['A'].timeSignal, \
        ps5000a.chs['B'].timeSignal, ps5000a.getfs(), f)
    angles = angles + [angle]
    amps = amps + [amp]

fig, ax1 = plt.subplots()
ax1.plot(freqs, angles, 'r.')
ax2 = ax1.twinx()
ax2.plot(freqs, amps, 'b.')
plt.show()

ps5000a.configurePSD(1, int(1953125/1))
plt.plot(ps5000a.getPSD('A', avg = 5, offset = -1))
plt.show()

FG.Sine(fc, 50e-3, ch = 2)
ps5000a.configurePSD(10, int(1953125/5))
FG.inst.close()
CH.PID_ver2(3.08, -0.02, -0.003, 0, fc, ps5000a, wd, initial = 1.06)

ps5000a = ps.picoscope()
ps5000a.defaultSetting()
ps5000a.DC('A')
ps5000a.AC('B')
ps5000a.AC('C')
ps5000a.DC('D')
ps5000a.ChangeRangeTo('A', 200)
ps5000a.ChangeRangeTo('B', 1000)
ps5000a.ChangeRangeTo('C', 50)
ps5000a.ChangeRangeTo('D', 5000)
ps5000a.configurePSD(1, int(1953125))
ps5000a.LPF20M('A', 1)
ps5000a.LPF20M('B', 1)
ps5000a.LPF20M('C', 1)
nT = 1953125
T = 0.5
ps5000a.configurePSD(T, nT)
FG = FuncGen.FuncGen()

FG.inst.write('SOUR1:BURS:STAT ON')
FG.inst.write('SOUR1:BURS:MODE TRIG')
FG.inst.write('SOUR1:BURS:NCYC 1')
FG.inst.write('SOUR1:BURS:PHAS 0')
FG.inst.write('TRIG1:SOUR BUS')

FG.inst.write('SOUR2:BURS:STAT ON')
FG.inst.write('SOUR2:BURS:MODE TRIG')
FG.inst.write('SOUR2:BURS:NCYC {NCYC:d}'.format(NCYC = int(400e3)))
FG.inst.write('SOUR2:BURS:PHAS 0')
FG.inst.write('TRIG2:SOUR BUS')

def drive_sin(t, f, amp, phi, t0, t1, DC = 0, window = None):
    t0_idx = np.argmin(np.abs(t - t0))
    t1_idx = np.argmin(np.abs(t - t1))
    drive = np.zeros(t.shape)
    if window is None:
        windowf = 1
    else:
        windowf = window(len(t[t0_idx:t1_idx+1]))
    drive[t0_idx : t1_idx+1] = DC + amp * np.multiply(windowf, np.sin(2*np.pi*f*t[t0_idx : t1_idx+1]+phi))
    return drive


N = int(6.4e6)
FGfs = int(3.2e6)
t0 = 0
t1 = 2 #in unit of second
FGt = (np.linspace(t0, t1, N+1))[:-1]

amp_s = 50e-3
mode_s_f = 201618.7
phi_s = 0/180*np.pi #at resoance there should be 90 phase shift between drive and response
ts0 = 0.1
ts1 = 0.11
drive_s = drive_sin(FGt, mode_s_f, amp_s, phi_s, ts0, ts1)

amp_i = 0
mode_i_f = 147618.1
phi_i = 0/180*np.pi
ti0 = 0.3
ti1 = 0.5
drive_i = drive_sin(FGt, mode_i_f, amp_i, phi_i, ti0, ti1)

amp_pa = 5
f_pa = mode_s_f + mode_i_f
phi_pa = 0/180*np.pi
tpa0 = ts1
tpa1 = ts1+1.8
drive_pa = drive_sin(FGt, f_pa, amp_pa, phi_pa, tpa0, tpa1)

amp_bs = 0
f_bs = mode_s_f - mode_i_f
phi_bs = 0/180*np.pi
tbs0 = ts1
tbs1 = ts1+1.8
drive_bs = drive_sin(FGt, f_bs, amp_bs, phi_bs, tbs0, tbs1)

drive = drive_s + drive_i + drive_pa + drive_bs
arb = FG.arb(drive, FGfs, 'test1', ch = 1)

lockinref = drive_sin(FGt, mode_s_f, 50e-3, 0, t0, t1) + drive_sin(FGt, mode_i_f, 50e-3, 0, t0, t1)
arb = FG.arb(lockinref, FGfs, 'test3', ch = 2)
FG.sync()


fs = ps5000a.getfs()
FG.inst.write('TRIG1:DEL 1')
FG.inst.write('TRIG2:DEL 1')
FG.inst.write('OUTP1:SYNC:MODE CARR')
FG.inst.write('SOUR:MARK:POIN {idx:d}'.format(idx = 4))
FG.inst.write('SOUR:MARK:POIN {idx:d}'.format(idx = np.argmin(np.abs(t - tbs1)) + 1))
FG.inst.write('SOUR:MARK:POIN {idx:d}'.format(idx = int(3.2e6)))

FG.inst.write('TRIG')
ps5000a.getTimeSignal('C', trigger = 'D', triggerThreshold = 8000)
t_ref = (np.linspace(0, 1, nT*2+1))[:-1]
ref_i = drive_sin(t_ref[:len(ps5000a.chs['C'].timeSignal)], mode_i_f, 1, phi_i, 0, 1)
splitN = 10
signal_is = np.array(np.split(ps5000a.chs['C'].timeSignal, splitN))
ref_is = np.array(np.split(ref_i, splitN))
amps = np.zeros(splitN)
angles = np.zeros(splitN)
for j in range(splitN):
    amps[j], angles[j] = CH.lock_in_amplifier(signal_is[j, :], ref_is[j,:], fs, mode_i_f)


FG.inst.write('TRIG')
FG.inst.write('*TRG')
ps5000a.getTimeSignal('C', trigger = 'D', triggerThreshold = 8000)
L = int(sum(ps5000a.chs['D'].timeSignal > 3)*2)
signal = ps5000a.chs['C'].timeSignal[:L]
ref = ps5000a.chs['A'].timeSignal[:L]
ampts, anglets = CH.lock_in_amplifier_t(signal, ref, fs, mode_s_f)
ampti, angleti = CH.lock_in_amplifier_t(signal, ref, fs, mode_i_f)
ps5000a_t = np.linspace(0,2,len(signal))

fig, ax1 = plt.subplots()
ax1.plot(ps5000a_t, ampts, 'b')
ax2 = ax1.twinx()
ax2.plot(ps5000a_t, anglets, 'r')
plt.show()

fig, ax1 = plt.subplots()
ax1.plot(ps5000a_t, ampti, 'b')
ax2 = ax1.twinx()
ax2.plot(ps5000a_t, angleti, 'r')
plt.show()

plt.plot(ps5000a_t, ampts)
plt.show()

plt.plot(ps5000a_t, ampis)
plt.show()

fig, ax1 = plt.subplots()
ax1.plot(ps5000a_t, ampts/(np.sqrt(sth)/mode_s_f), 'b')
ax2 = ax1.twinx()
ax2.plot(ps5000a_t, ampti/(np.sqrt(ith)/mode_i_f), 'r')
plt.show()

psd_signal = np.log10(ms.Single_sided_PSD(signal, fs))
psd_ref = np.log10(ms.Single_sided_PSD(ref, fs))
plt.plot(psd_signal)
plt.plot(psd_ref)
plt.show()

for t1 in range(tbs0, tbs0 + 0.1, 0.1/100):
    drive_bs = drive_sin(t, f_bs, amp_bs, phi_bs, tbs0, t1)
    drive = drive_s + drive_i + drive_pa + drive_bs
    arb = FG.arb(drive, fs, 'test', ch = 1)
    sleep(5)
    FG.inst.write('TRIG')
    ps5000a.getTimeSignal('C', trigger = 'D', triggerThreshold = 8000)

for phi in range(0, np.pi, np.pi/100):
    drive_bs = drive_sin(t, f_bs, amp_bs, phi, tbs0, tbs1)
    drive = drive_s + drive_i + drive_pa + drive_bs
    arb = FG.arb(drive, fs, 'test', ch = 1)
    sleep(10)
    FG.inst.write('TRIG')
    ps5000a.getTimeSignal('C', trigger = 'D', triggerThreshold = 8000)

[201618.7, 223303.2]
[201618.7, 0.9593, 1.3206e6] 2.6182968801034233e-11
[223303.2, 6.65, 2.1091e5]
[147618.1, 2.0845, 4.4495e+05] 2.4980384423344627e-10
[123746.3, 24.25, 3.2056e+04]

fs = 201618.7
fi = 147618.1
env = { 'N':int(6.4e6),
        'FGfs':int(3.2e6),
        't0':0,
        't1':2}
s = {   'amp':50e-3,
        'f':fs,
        'phi' : 0,
        't0':0.1,
        't1':0.11,
        'DC':0}
i = {   'amp':0,
        'f':fi,
        'phi' : 0,
        't0':0.1,
        't1':0.9,
        'DC':0}
pa = {  'amp':0,
        'f':fs + fi,
        'phi' : 0,
        't0':0.1,
        't1':0.9,
        'DC':0}
bs = {  'amp':1,
        'f':fs - fi,
        'phi' : 0,
        't0':0.11,
        't1':1.81,
        'DC':0}
params = {'env':env, 's':s, 'i':i, 'pa':pa, 'bs':bs}

def thermal_sum(freqs, psd, fc, bandwidth):
    idx0 = np.argmin(np.absolute(freqs - (fc - bandwidth/2)))
    idx1 = np.argmin(np.absolute(freqs - (fc + bandwidth/2)))
    psd_lin = np.power(10, psd[idx0:idx1]/10)/20
    E = sum(psd_lin)
    return E

def FGpulses(FG, params):
    N = params['env']['N']
    FGfs = params['env']['FGfs']
    t0 = params['env']['t0']
    t1 = params['env']['t1'] #in unit of second
    FGt = (np.linspace(t0, t1, N+1))[:-1]
    drive_s = drive_sin(FGt, params['s']['f'], params['s']['amp'], \
        params['s']['phi'], params['s']['t0'], params['s']['t1'], \
        DC = params['s']['DC'], window = lambda N : scipy.signal.tukey(N, alpha = 0.1))
    drive_i = drive_sin(FGt, params['i']['f'], params['i']['amp'], \
        params['i']['phi'], params['i']['t0'], params['i']['t1'], \
        DC = params['i']['DC'], window = lambda N : scipy.signal.tukey(N, alpha = 0.1))
    drive_pa = drive_sin(FGt, params['pa']['f'], params['pa']['amp'], \
        params['pa']['phi'], params['pa']['t0'], params['pa']['t1'], \
        DC = params['pa']['DC'], window = lambda N : scipy.signal.tukey(N, alpha = 0.1))
    drive_bs = drive_sin(FGt, params['bs']['f'], params['bs']['amp'], \
        params['bs']['phi'], params['bs']['t0'], params['bs']['t1'], \
        DC = params['bs']['DC'], window = lambda N : scipy.signal.tukey(N, alpha = 0.1))
    drive = drive_s + drive_i + drive_pa + drive_bs
    arb = FG.arb(drive, FGfs, 'test1', ch = 1)
    lockinref = drive_sin(FGt, params['s']['f'], 50e-3, 0, \
        params['env']['t0'], params['env']['t1']) + \
        drive_sin(FGt, params['i']['f'], 50e-3, 0, \
        params['env']['t0'], params['env']['t1'])
    arb = FG.arb(lockinref, FGfs, 'test3', ch = 2)

def FallEdge10(trig):
    idx = int(len(trig)/2)
    lo = 0
    hi = len(trig)
    while not np.array_equal(trig[idx-1:idx+2], [1,1,0]):
        idx = int((lo + hi)/2)
        if trig[idx] == 0:
            hi = idx
        else:
            lo = idx
    return idx

list = np.arange(-0.5, 0.6, 0.2)
ampti_list = [[]]*len(list)
ampts_list = [[]]*len(list)
angleti_list = [[]]*len(list)
anglets_list = [[]]*len(list)
i = 0

for deltaf in list:
    params['s']['f'] = fs + deltaf
    params['bs']['f'] = fs + deltaf - fi
    FGpulses(FG, params)
    FG.inst.write('*TRG')
    ps5000a.getTimeSignal('C', trigger = 'D', triggerThreshold = 8000)
    temp = ps5000a.chs['D'].timeSignal > 3
    temp = temp.astype(int)
    L = int((FallEdge10(temp)+1)*2)
    signal = ps5000a.chs['C'].timeSignal[:L]
    ref = ps5000a.chs['A'].timeSignal[:L]
    ampts, anglets = CH.lock_in_amplifier_t(signal, ref, ps5000a.getfs(), params['s']['f'])
    ampti, angleti = CH.lock_in_amplifier_t(signal, ref, ps5000a.getfs(), params['i']['f'])
    ampts_list[i] = ampts
    ampti_list[i] = ampti
    angleti_list[i] = angleti
    anglets_list[i] = anglets
    i = i+1
    #ps5000a_t = np.linspace(0,2,len(signal))
    time.sleep(5)


amp_list = np.arange(1, 1.2, 0.2)
deltaf_list = np.arange(0.5, 1, 0.1)
ampti_list = [[]]*len(list)
ampts_list = [[]]*len(list)
angleti_list = [[]]*len(list)
anglets_list = [[]]*len(list)
i = 0

for amp in list:
    params['bs']['amp'] = amp
    FGpulses(FG, params)
    FG.inst.write('*TRG')
    ps5000a.getTimeSignal('C', trigger = 'D', triggerThreshold = 8000)
    temp = ps5000a.chs['D'].timeSignal > 3
    temp = temp.astype(int)
    L = int((FallEdge10(temp)+1)*2)
    signal = ps5000a.chs['C'].timeSignal[:L]
    ref = ps5000a.chs['A'].timeSignal[:L]
    ampts, anglets = CH.lock_in_amplifier_t(signal, ref, ps5000a.getfs(), params['s']['f'])
    ampti, angleti = CH.lock_in_amplifier_t(signal, ref, ps5000a.getfs(), params['i']['f'])
    ampts_list[i] = ampts
    ampti_list[i] = ampti
    angleti_list[i] = angleti
    anglets_list[i] = anglets
    i = i+1
    #ps5000a_t = np.linspace(0,2,len(signal))
    time.sleep(5)

params['s']['amp'] = 0.05
params['bs']['DC'] = 0
amp_list = np.array([1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 0.1, 0.2, 0.5, 1, 1.5, 2, 3])
amp_list = [0.1, 0.5, 1, 2, 4]
deltaf_list = np.arange(-0.5, 0.6, 0.1)
phi_list = np.linspace(0, np.pi, 19)[:-1]
phi_list = [0]

L1 = len(amp_list)
L2 = len(phi_list)
ampti_list = np.zeros([L1, L2]).tolist()
ampts_list = np.zeros([L1, L2]).tolist()
angleti_list = np.zeros([L1, L2]).tolist()
anglets_list = np.zeros([L1, L2]).tolist()
i = 0
for amp in amp_list:
    params['bs']['amp'] = amp
    j = 0
    for phi in phi_list:
        params['i']['f'] = fi + deltaf
        params['bs']['phi'] = phi
        FGpulses(FG, params)
        FG.inst.write('*TRG')
        ps5000a.getTimeSignal(trigger = 'D', triggerThreshold = 8000)
        temp = ps5000a.chs['D'].timeSignal > 3
        temp = temp.astype(int)
        L = int((FallEdge10(temp)+1)*2)
        signal = ps5000a.chs['C'].timeSignal[:L]
        ref = ps5000a.chs['A'].timeSignal[:L]
        plt.plot(np.log10(ms.Single_sided_PSD(ps5000a.chs['B'].timeSignal[:L], ps5000a.getfs())))
        ampts, anglets = CH.lock_in_amplifier_t(signal, ref, ps5000a.getfs(), params['s']['f'])
        ampti, angleti = CH.lock_in_amplifier_t(signal, ref, ps5000a.getfs(), params['i']['f'])
        ampts_list[i][j] = ampts
        ampti_list[i][j] = ampti
        #angleti_list[i][j] = angleti
        #anglets_list[i][j] = anglets
        id = int(datetime.timestamp(datetime.now()))
        signal.tofile(wd + 'OLtimeSignal_id={id:d}.bin'.format(id = id), sep = '')
        ref.tofile(wd + 'OLtimeRef_id={id:d}.bin'.format(id = id), sep = '')
        ampts.tofile(wd + 'S_amp_id={id:d}.bin'.format(id = id), sep = '')
        ampti.tofile(wd + 'I_amp_id={id:d}.bin'.format(id = id), sep = '')
        anglets.tofile(wd + 'S_phi_id={id:d}.bin'.format(id = id), sep = '')
        angleti.tofile(wd + 'I_phi_id={id:d}.bin'.format(id = id), sep = '')
        with open(wd + 'params_id={id:d}.json'.format(id = id), "w") as outfile:
            json.dump(params, outfile)
        time.sleep(5)
        j = j+1
    i = i+1

fig, ax1 = plt.subplots()
ax1.plot(ampti/(np.sqrt(sth)/fs), 'b')
ax2 = ax1.twinx()
ax2.plot(ampts/(np.sqrt(ith)/fi), 'r')
plt.show()

fig, ax1 = plt.subplots()
i = 7
ax1.plot(ampti_list[0][i]/(np.sqrt(sth)/fs), 'b')
ax2 = ax1.twinx()
ax2.plot(ampts_list[0][i]/(np.sqrt(ith)/fi), 'r')
plt.show()

for i in range(len(phi_list)):
    plt.plot(ampti_list[0][i])

for i in range(len(phi_list)):
    plt.plot(ampts_list[0][i])

for i in range(len(deltaf_list)):
    plt.plot(ampti_list[0][i])

for i in range(len(amp_list)):
    plt.plot(ampti_list[i][0])

for i in range(len(deltaf_list)):
    plt.plot(ampts_list[0][i])

for i in range(len(amp_list)):
    plt.plot(ampts_list[i][0])

for i in range(len(amp_list)):
    plt.plot(ampti_list[i][5])

for i in range(15):
    plt.plot(ampti_list[0][5])

for i in range(len(amp_list)):
    plt.plot(ampts_list[i][2])
