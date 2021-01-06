import sys, argparse
import os
import numpy as np
from picosdk.ps5000a import ps5000a as ps
from picosdk.functions import assert_pico_ok, adc2mV
from ctypes import cdll, c_long, c_ulong, c_uint32, byref, \
create_string_buffer, c_bool, c_char_p, c_int, c_int16, c_int32,\
c_double, sizeof, c_voidp, c_float
from scipy import signal
from math import floor, ceil
import matplotlib.pyplot as plt
import csv
import time

class ch:

    def __init__(self, channel, device_handle):
        self.channel = ps.PS5000A_CHANNEL["PS5000A_CHANNEL_"+channel]
        self.coupling_type = ps.PS5000A_COUPLING["PS5000A_AC"]
        self.range = ps.PS5000A_RANGE["PS5000A_50MV"]
        self.enabled = False
        self.offset = 0
        self.device_handle = device_handle
        self.buffer = None
        self.PSD = None
        self.timeSignal = None

    def Enable(self):
        self.enabled = True
        return ps.ps5000aSetChannel(self.device_handle, self.channel, self.enabled,\
            self.coupling_type, self.range, self.offset)

    def Disable(self):
        self.enabled = False
        return ps.ps5000aSetChannel(self.device_handle, self.channel, self.enabled,\
            self.coupling_type, self.range, self.offset)

    def set(self):
        return ps.ps5000aSetChannel(self.device_handle, self.channel, self.enabled,\
            self.coupling_type, self.range, self.offset)


class picoscope:

    def __init__(self, workingDirectory = r'D:\temp'):
        self.chandle = c_int16()
        self.status = {}
        self.resolution = ps.PS5000A_DEVICE_RESOLUTION["PS5000A_DR_14BIT"]
        self.status["openunit"] = ps.ps5000aOpenUnit(byref(self.chandle), \
            None, self.resolution)
        assert_pico_ok(self.status["openunit"])
        self.chs = {'A':ch('A', self.chandle), 'B':ch('B', self.chandle),\
            'C':ch('C', self.chandle), 'D':ch('D', self.chandle)}
        self.workingDirectory = workingDirectory + '\\'
        self.psVoltageRange = {0 : 10, 1 : 20, 2 : 50, 3 : 100, 4 : 200, \
            5 : 500, 6 : 1000, 7 : 2000, 8 : 5000, 9 : 10000, 10 : 20000}
        self.defaultSetting()

    def defaultSetting(self):
        self.timebase = int(6000)
        self.timeInternalns = c_float()
        self.bufferSize = c_int32()
        self.maxSamples = 10000
        self.status["getTimebase2"] = ps.ps5000aGetTimebase2(self.chandle, \
            self.timebase, self.maxSamples, byref(self.timeInternalns),\
            byref(self.bufferSize), 0)
        assert_pico_ok(self.status["getTimebase2"])

        self.chs['A'].coupling_type = ps.PS5000A_COUPLING["PS5000A_AC"]
        self.chs['A'].range = ps.PS5000A_RANGE["PS5000A_50MV"]
        self.status["setChA"] = self.chs['A'].Enable()
        assert_pico_ok(self.status["setChA"])
        self.AutoRange('A')

        self.chs['B'].coupling_type = ps.PS5000A_COUPLING["PS5000A_AC"]
        self.chs['B'].range = ps.PS5000A_RANGE["PS5000A_50MV"]
        self.status["setChB"] = self.chs['B'].Enable()
        assert_pico_ok(self.status["setChB"])
        self.AutoRange('B')

        self.chs['C'].coupling_type = ps.PS5000A_COUPLING["PS5000A_AC"]
        self.chs['C'].range = ps.PS5000A_RANGE["PS5000A_50MV"]
        self.status["setChC"] = self.chs['C'].Enable()
        assert_pico_ok(self.status["setChC"])
        self.AutoRange('C')

    def getPSD(self, binWidth, spectrumRange, nmbOfAvg):
        #PSD in unit of dBm/bin, as comparison, PICO gives dBm/bin
        #gives singles-sided PSD, but energy on such single side is full energy
        #which is the same as PICO
        spectrumRange = spectrumRange * 2
        requiredMeasureTime = 1 / binWidth
        requiredSamplingInterval = 1 / spectrumRange

        self.timebase = floor(requiredSamplingInterval * 62500000 + 3)
        self.timeInternalns = c_float()
        self.bufferSize = c_int32()
        self.maxSamples = ceil(spectrumRange / binWidth)
        self.status["getTimebase2"] = ps.ps5000aGetTimebase2(self.chandle, \
            self.timebase, self.maxSamples, byref(self.timeInternalns),\
            byref(self.bufferSize), 0)
        assert_pico_ok(self.status["getTimebase2"])
        assert self.timeInternalns.value < requiredSamplingInterval * 1e9

        for i in range(nmbOfAvg):
            self.getTimeSignal()
            for key in self.chs:
                if self.chs[key].enabled:
                    f, newPSD = signal.periodogram(self.chs[key].timeSignal, 1 / (self.timeInternalns.value / 1e9), \
                        window = signal.get_window('blackman', len(self.chs[key].timeSignal)))
                    if i == 0:
                        self.chs[key].PSD = np.array(newPSD)
                    else:
                        self.chs[key].PSD = self.chs[key].PSD + newPSD
        self.f = f
        for key in self.chs:
            if self.chs[key].enabled:
                self.chs[key].PSD = self.chs[key].PSD/nmbOfAvg
                self.chs[key].PSD = self.chs[key].PSD/(self.timeInternalns.value/1e9*len(self.chs[key].timeSignal))
                self.chs[key].PSD = 10 * np.log10(20 * self.chs[key].PSD)

    def savePSD(self, filename):
        if not os.path.exists(self.workingDirectory):
            os.makedirs(self.workingDirectory)
        rows = []
        fileds = ['f(Hz)', 'PSD_A(dBm)', 'PSD_B(dBm)', 'PSD_C(dBm)', 'PSD_D(dBm)']
        for key in self.chs:
            if not self.chs[key].enabled:
                self.chs[key].PSD = [None]* len(self.f)
        for i in range(0,len(self.f)):
            rows.append([self.f[i], self.chs['A'].PSD[i], self.chs['B'].PSD[i], self.chs['C'].PSD[i], self.chs['D'].PSD[i]])
        with open(self.workingDirectory+filename+'.csv', mode='w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(fileds)
            csvwriter.writerows(rows)

    def saveTS(self, filename):
        if not os.path.exists(self.workingDirectory):
            os.makedirs(self.workingDirectory)
        rows = []
        fileds = ['t(s)', 'A(V)', 'B(V)', 'C(V)', 'D(V)']
        for key in self.chs:
            if not self.chs[key].enabled:
                self.chs[key].timeSignal = [None]* len(self.t)
        for i in range(0,len(self.f)):
            rows.append([self.t[i], self.chs['A'].timeSignal[i], self.chs['B'].timeSignal[i], \
                self.chs['C'].timeSignal[i], self.chs['D'].timeSignal[i]])
        with open(self.workingDirectory+filename+'.csv', mode='w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(fileds)
            csvwriter.writerows(rows)

    def plot(self, channel):
        plt.plot(self.memory[0],self.memory[1])
        plt.xlabel('f(Hz)')
        plt.ylabel('Power(dBm)')
        plt.show()

    def close(self):
        self.status["stop"] = ps.ps5000aStop(self.chandle)
        assert_pico_ok(self.status["stop"])

        self.status["close"]=ps.ps5000aCloseUnit(self.chandle)
        assert_pico_ok(self.status["close"])

    def AutoRange(self, channel):
        orginalMaxSamples = self.maxSamples
        self.maxSamples = 10000
        while self.chs[channel].range < max(self.psVoltageRange.keys()):
            overrange = False
            timeSignal = self.getTimeSignal(channel)
            if max(map(abs, timeSignal)) > 0.95* self.psVoltageRange[self.chs[channel].range]/1000:
                overrange = True
            if overrange:
                self.chs[channel].range += 1
                self.status["setCh"+channel] = self.chs[channel].set()
                assert_pico_ok(self.status["setCh"+channel])
            else:
                break

        while self.chs[channel].range > min(self.psVoltageRange.keys()):
            toosmall = False
            self.chs[channel].range = self.chs[channel].range - 1
            self.status["setCh"+channel] = self.chs[channel].set()
            assert_pico_ok(self.status["setCh"+channel])
            timeSignal = self.getTimeSignal(channel)
            if max(map(abs, timeSignal)) < 0.95* self.psVoltageRange[self.chs[channel].range]/1000:
                toosmall = True
            if toosmall:
                self.chs[channel].range = self.chs[channel].range - 1
            else:
                self.chs[channel].range = self.chs[channel].range + 1
                self.status["setCh"+channel] = self.chs[channel].set()
                assert_pico_ok(self.status["setCh"+channel])
                break
        self.maxSamples = orginalMaxSamples

    def ChangeRangeTo(self, channel, Range):
        keys = list(self.psVoltageRange.keys())
        values = list(self.psVoltageRange.values())
        if Range in values:
            self.chs[channel].range = keys[values.index(Range)]
            self.status["setCh"+channel] = self.chs[channel].set()
            assert_pico_ok(self.status["setCh"+channel])
        else:
            print('Not a acceptable range. Select from (intger):')
            print(values)

    def DC(self, channel):
        self.chs[channel].coupling_type = ps.PS5000A_COUPLING["PS5000A_DC"]
        self.status["setCh"+channel] = self.chs[channel].set()
        assert_pico_ok(self.status["setCh"+channel])

    def AC(self, channel):
        self.chs[channel].coupling_type = ps.PS5000A_COUPLING["PS5000A_AC"]
        self.status["setCh"+channel] = self.chs[channel].set()
        assert_pico_ok(self.status["setCh"+channel])

    def getTimeSignal(self, channel = None):
        # get all channel data but only return required
        nMaxSamples = c_int32()
        ps.ps5000aMemorySegments(self.chandle, self.nEnabledChannels(), byref(nMaxSamples))

        if self.maxSamples > nMaxSamples.value:
            raise Exception("samples will be larger than memory")

        preTriggerSamples = 100
        self.status["runBlock"] = ps.ps5000aRunBlock(self.chandle, \
            preTriggerSamples, self.maxSamples, self.timebase, None, 0, None, None)
        assert_pico_ok(self.status["runBlock"])

        self.status["getTimebase2"] = ps.ps5000aGetTimebase2(self.chandle, \
            self.timebase, self.maxSamples, byref(self.timeInternalns),\
            None, 0)
        assert_pico_ok(self.status["getTimebase2"])
        self.t = np.linspace(0, (self.maxSamples-1)*self.timeInternalns.value/1e9, self.maxSamples)

        ready = c_int16(0)
        check = c_int16(0)
        while ready.value == check.value:
            self.status["isReady"] = ps.ps5000aIsReady(self.chandle, byref(ready))

        for key in self.chs:
            if self.chs[key].enabled:
                self.chs[key].buffer = (c_int16 * self.maxSamples)()
                self.status["setDataBuffer"+key] = ps.ps5000aSetDataBuffer(self.chandle, \
                    self.chs[key].channel, byref(self.chs[key].buffer), self.maxSamples, 0, 0)
                assert_pico_ok(self.status["setDataBufferA"])

        overflow = c_int16()
        cmaxSamples = c_uint32(self.maxSamples)
        self.status["getValues"] = ps.ps5000aGetValues(self.chandle, 0 , \
        byref(cmaxSamples), 0, 0, 0, byref(overflow))
        assert_pico_ok(self.status["getValues"])

        maxADC = c_int16()
        self.status["maximumValue"] = ps.ps5000aMaximumValue(self.chandle, byref(maxADC))
        assert_pico_ok(self.status["maximumValue"])

        for key in self.chs:
            if self.chs[key].enabled:
                self.chs[key].timeSignal = adc2mV(self.chs[key].buffer, \
                    self.chs[key].range, maxADC)
                self.chs[key].timeSignal = [x * 1e-3 for x in self.chs[key].timeSignal]
                self.chs[key].timeSignal = np.array(self.chs[key].timeSignal)
        if channel is not None:
            return self.chs[channel].timeSignal

    def timeSignalStable(self, channel):
        S = self.timeSignal(channel)
        return (np.std(S) < 0.01 * np.mean(S) )

    def nEnabledChannels(self):
        n = 0
        for key in self.chs:
            if self.chs[key].enabled:
                n = n + 1
        return n


def main(argv):
    flag = 0

if __name__ == "__main__":
    main(sys.argv[1:])
