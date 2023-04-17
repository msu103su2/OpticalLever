import pyvisa as py
from scipy.signal import chirp, spectrogram
import numpy as np
import os
import matplotlib.pyplot as plt

class FuncGen:

    def __init__(self):
        self.rm = py.ResourceManager()
        self.inst = self.rm.open_resource('USB0::0x0957::0x2C07::MY57802003::INSTR')

    def Sine(self, freq, amp, offset = 0, ch = 1):
        self.inst.write('SOUR'+str(ch)+':APPL:SIN '+str(freq)+','+str(amp)+','+str(offset))

    def Noise(self, freq, amp, offset = 0, ch = 1):
        self.inst.write('SOUR'+str(ch)+'APPL:NOIS '+str(freq)+','+str(amp)+','+str(offset))

    def chirp(self, arbname, fstart, fend, fs, N, Vmax, Vmin, phi = 0):
        t = np.array(list(range(int(N))))/fs
        data = chirp(t, fstart, t[-1], fend, method='linear', phi = phi)
        data = (data+1)/2*(Vmax- Vmin)+Vmin
        self.arb(data, fs, arbname)

    def chirps(self, arbname, fstarts, fends, fs, N, Vmaxs, Vmins, phis = None, ch = 1):
        t = np.array(list(range(int(N))))/fs
        data = np.zeros(int(N))
        if phis is None:
            phis = [0]*len(fstarts)
        for f1, f2, Vmax, Vmin, phi in zip(fstarts, fends, Vmaxs, Vmins, phis):
            temp = chirp(t, f1, t[-1], f2, method='linear', phi = phi)
            temp = (temp+1)/2*(Vmax- Vmin)-(Vmax-Vmin)/2
            data = data + temp
        data = data + (Vmaxs[0]+Vmins[0])/2
        return self.arb(data, fs, arbname, ch = ch)

    def sweep(self, fstart, fend, sweep_time, amp = 0.25, dc = 0.365, \
                hold_time = 0, return_time = 0, function = 'SIN', \
                trigger = False, triggerF = 0):
        self.inst.write('OUTPUT1 OFF')
        self.inst.write('SOUR1:FUNC {f:s}'.format(f = function))
        self.inst.write('SOUR1:FREQ {f:f}'.format(f = fstart))
        self.inst.write('SOUR1:FREQ:STAR {f:f}'.format(f = fstart))
        self.inst.write('SOUR1:FREQ:STOP {f:f}'.format(f = fend))
        self.inst.write('SOUR1:VOLT {amp:f}'.format(amp = amp))
        self.inst.write('SOUR1:VOLT:OFFSET {dc:f}'.format(dc = dc))
        self.inst.write('SOUR1:SWE:TIME {t:f}'.format(t = sweep_time))
        self.inst.write('SOUR1:FREQ:MODE SWE')
        self.inst.write('OUTPUT ON')

    def Sines(self, freqs, amplitudes):
        Tmax = 1/freqs.min()
        T = 20* Tmax
        fs = np.around(10*2*np.pi*freqs.max())
        x = np.arange(0, T, 1/fs)
        V = np.zeros(x.shape)
        for freq, amplitude in zip(freqs, amplitudes):
            V = V + amplitude*np.sin(2*np.pi*freq*x)
        return V, fs

    def arb(self, data, fs, arbname, ch = 1, clearMem = 1):
        #data in voltages
        self.inst.write('OUTPUT{:d} OFF'.format(ch))
        if clearMem:
            self.inst.write('SOUR{ch:d}:DATA:VOL:CLE'.format(ch = ch))

        self.inst.write('FORM:BORD SWAP')

        max = data.max()
        min = data.min()
        offset = (max + min)/2
        amp = max - offset

        data = (data - offset)/amp*32767
        data = data.astype(np.int16)
        data = data.tolist()

        del self.inst.timeout
        self.inst.write_binary_values('SOUR{ch:d}:DATA:ARB1:DAC {arbname},'.format(arbname = arbname, ch = ch), data, datatype = 'h')
        self.inst.write('*WAI')
        self.inst.write('SOUR{ch:d}:FUNC:ARB {arbname}'.format(arbname = arbname, ch = ch))
        self.inst.write('SOUR{ch:d}:FUNC:ARB:SRAT {fs:d}'.format(fs = int(fs), ch = ch))
        self.inst.write('SOUR{ch:d}:VOLT {amp:f}'.format(amp = amp, ch = ch))
        self.inst.write('SOUR{ch:d}:VOLT:OFFSET {offset:f}'.format(offset = offset, ch = ch))
        self.inst.write('SOUR{ch:d}:FUNC:ARB:FILT NORM'.format(ch = ch))
        self.inst.write('SOUR{ch:d}:FUNC ARB'.format(ch = ch))
        self.inst.write('OUTPUT{ch:d} ON'.format(ch = ch))
        self.inst.timeout = 2e3
        return data

    def sync(self):
        self.inst.write('OUTP:SYNC ON')

    def phi(self, degree, ch = 2):
        self.inst.write('SOUR{ch:d}:PHAS {phi:f}'.format(ch = ch, phi = degree))

    def DC(self, dc, ch = 1):
        self.inst.write('SOUR'+str(ch)+':APPL:DC DEF, DEF, {dc:.3f} V'.format(dc = dc))

    def sweep_trigger(self, fstart, fend, markerF, ch = 1):
        self.inst.write('SOUR1:FREQ:STAR {:f}'.format(fstart))
        self.inst.write('SOUR1:FREQ:STOP {:f}'.format(fend))
        self.inst.write('MARKer:FREQuency {:f}'.format(markerF))
        self.inst.write('TRIGger1:SOURce IMM')


def Ndigits(number):
    count = 0
    while number >0:
        number = number // 10
        count = count + 1
    return count

def utf8len(s):
    return len(s.encode('utf-8'))
