import sys, argparse
import os
import numpy as np
from picosdk.ps5000a import ps5000a as ps
from picosdk.functions import assert_pico_ok, mV2adc
from ctypes import cdll, c_long, c_ulong, c_uint32, byref, \
create_string_buffer, c_bool, c_char_p, c_int, c_int16, c_int32,\
c_double, sizeof, c_voidp, c_float
from scipy import signal
from math import floor, ceil
import matplotlib.pyplot as plt
import csv
import time
import math
import re

class ch:

    def __init__(self, channel, device_handle):
        self.channel = ps.PS5000A_CHANNEL["PS5000A_CHANNEL_"+channel]
        self.coupling_type = ps.PS5000A_COUPLING["PS5000A_AC"]
        self.range = ps.PS5000A_RANGE["PS5000A_50MV"]
        self.enabled = False
        self.offset = 0
        self.device_handle = device_handle
        self.buffer = None
        self.timeSignal = None
        self.overflow = False

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

    def clearMem(self):
        self.timeSignal = None

    def getRange(self):
        s = list(ps.PS5000A_RANGE.keys())[list(ps.PS5000A_RANGE.values()).index(self.range)]
        result = re.findall("PS5000A_([0-9]+)([A-Z]+)", s)
        range = int(result[0][0])
        if result[0][1] == 'V':
            range = int(range*1e3)
        return range

    def getConfigureInfo(self):
        x = {
            'channel':list(ps.PS5000A_CHANNEL.keys())[list(ps.PS5000A_CHANNEL.values()).index(self.channel)],
            'range':list(ps.PS5000A_RANGE.keys())[list(ps.PS5000A_RANGE.values()).index(self.range)],
            'coupling_type':list(ps.PS5000A_COUPLING.keys())[list(ps.PS5000A_COUPLING.values()).index(self.coupling_type)],
            'enabled':self.enabled
        }
        return x

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
        self.timebase = int(6000)
        self.timeInternalns = c_float()
        self.maxSamples = 10000
        self.t = None
        self.f = None

    def defaultSetting(self):
        self.status["getTimebase2"] = ps.ps5000aGetTimebase2(self.chandle, \
            self.timebase, self.maxSamples, byref(self.timeInternalns),\
            None, 0)
        assert_pico_ok(self.status["getTimebase2"])

        self.chs['A'].coupling_type = ps.PS5000A_COUPLING["PS5000A_AC"]
        self.chs['A'].range = ps.PS5000A_RANGE["PS5000A_50MV"]
        self.status["setChA"] = self.chs['A'].Enable()
        assert_pico_ok(self.status["setChA"])

        self.chs['B'].coupling_type = ps.PS5000A_COUPLING["PS5000A_AC"]
        self.chs['B'].range = ps.PS5000A_RANGE["PS5000A_50MV"]
        self.status["setChB"] = self.chs['B'].Enable()
        assert_pico_ok(self.status["setChB"])

        self.chs['C'].coupling_type = ps.PS5000A_COUPLING["PS5000A_AC"]
        self.chs['C'].range = ps.PS5000A_RANGE["PS5000A_50MV"]
        self.status["setChC"] = self.chs['C'].Enable()
        assert_pico_ok(self.status["setChC"])

        self.chs['D'].coupling_type = ps.PS5000A_COUPLING["PS5000A_AC"]
        self.chs['D'].range = ps.PS5000A_RANGE["PS5000A_50MV"]
        self.status["setChD"] = self.chs['D'].Enable()
        assert_pico_ok(self.status["setChD"])

    def configurePSD(self, binWidth, spectrumRange, nTimes = None):
        spectrumRange = spectrumRange * 2
        requiredSamplingInterval = 1 / spectrumRange
        self.timebase = floor(requiredSamplingInterval * 125000000 + 2) #14bit
        #self.timebase = floor(requiredSamplingInterval * 62500000 + 3) #16bit

        self.status["getTimebase2"] = ps.ps5000aGetTimebase2(self.chandle, \
            self.timebase, self.maxSamples, byref(self.timeInternalns),\
            None, 0)
        assert_pico_ok(self.status["getTimebase2"])

        requiredMeasureTime = 1 / binWidth
        self.maxSamples = ceil(requiredMeasureTime *1e9 / self.timeInternalns.value)

        assert self.timeInternalns.value <= requiredSamplingInterval * 1e9

        if bool(nTimes):
            nMaxSamples = c_long()
            ps.ps5000aMemorySegments(self.chandle, 1, byref(nMaxSamples))
            print(nMaxSamples)
            nMaxSamples.value = math.floor(nMaxSamples.value/self.nEnabledChannels())
            if nTimes == 'Max':
                nTimes = math.floor(nMaxSamples.value/(self.maxSamples))-25
            print(nTimes, self.maxSamples)
            assert self.maxSamples * nTimes < nMaxSamples.value
            self.maxSamples = nTimes * self.maxSamples
            return nTimes, int(self.maxSamples/nTimes)

    def PSDfromTS(self, timeSignal, fs):
        #PSD in unit of dBm/bin, as comparison, PICO gives dBm/bin
        #gives singles-sided PSD, but energy on such single side is full energy
        #which is the same as PICO
        f, PSD = signal.periodogram(timeSignal, fs, \
            window = signal.get_window('blackman', len(timeSignal)))
        self.f = f
        PSD = PSD/(self.timeInternalns.value/1e9*len(timeSignal)) #converting from /Hz to /Bin
        return PSD

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
        fileds = ['A(V)', 'B(V)', 'C(V)']
        for i in range(0,len(self.t)):
            rows.append([self.chs['A'].timeSignal[i], self.chs['B'].timeSignal[i], \
                self.chs['C'].timeSignal[i]])
        nprows = np.array(rows)
        nprows.tofile(self.workingDirectory+filename+'.bin', sep = '')

    def close(self):
        self.status["stop"] = ps.ps5000aStop(self.chandle)
        assert_pico_ok(self.status["stop"])

        self.status["close"]=ps.ps5000aCloseUnit(self.chandle)
        assert_pico_ok(self.status["close"])

    def AutoRange(self, channel):
        orginalMaxSamples = self.maxSamples
        self.maxSamples = 10000
        while self.chs[channel].range < max(self.psVoltageRange.keys()):
            self.getTimeSignal(reportOverflow = False)
            if self.chs[channel].overflow:
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
            timeSignal = self.getTimeSignal(channel, reportOverflow = False)
            if max(map(abs, timeSignal)) < 0.95* self.psVoltageRange[self.chs[channel].range]/1000:
                toosmall = True
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

    def getTimeSignal(self, channel = None, check = True, trigger = None, triggerThreshold = 0, reportOverflow = True):
        # get all channel data but only return required
        if check:
            nMaxSamples = c_long()
            ps.ps5000aMemorySegments(self.chandle, 1, byref(nMaxSamples))
            nMaxSamples.value = math.floor(nMaxSamples.value/self.nEnabledChannels())
            if self.maxSamples > nMaxSamples.value:
                raise Exception("samples will be larger than memory")

            self.status["getTimebase2"] = ps.ps5000aGetTimebase2(self.chandle, \
                self.timebase, self.maxSamples, byref(self.timeInternalns),\
                None, 0)
            assert_pico_ok(self.status["getTimebase2"])

        if trigger is not None:
            source = ps.PS5000A_CHANNEL["PS5000A_CHANNEL_{trigCh}".format(trigCh = trigger)]
            self.status["trigger"] = ps.ps5000aSetSimpleTrigger(self.chandle, 1, source, triggerThreshold, 2, 0, 0)
            assert_pico_ok(self.status["trigger"])

        preTriggerSamples = 0
        self.status["runBlock"] = ps.ps5000aRunBlock(self.chandle, \
            preTriggerSamples, self.maxSamples, self.timebase, None, 0, None, None)
        assert_pico_ok(self.status["runBlock"])

        ready = c_int16(0)
        check = c_int16(0)
        while ready.value == check.value:
            self.status["isReady"] = ps.ps5000aIsReady(self.chandle, byref(ready))
            time.sleep(1e-3)

        for key in self.chs:
            if self.chs[key].enabled:
                self.chs[key].buffer = (c_int16 * self.maxSamples)()
                self.status["setDataBuffer"+key] = ps.ps5000aSetDataBuffer(self.chandle, \
                    self.chs[key].channel, byref(self.chs[key].buffer), self.maxSamples, 0, 0)
                assert_pico_ok(self.status["setDataBuffer"+key])


        overflow = c_int16()
        cmaxSamples = c_uint32(self.maxSamples)
        self.status["getValues"] = ps.ps5000aGetValues(self.chandle, 0 , \
        byref(cmaxSamples), 0, 0, 0, byref(overflow))
        assert_pico_ok(self.status["getValues"])

        overflow = '{0:04b}'.format(overflow.value)
        chOF = [bool(int(i)) for i in overflow]
        self.chs['A'].overflow = chOF[-1]
        self.chs['B'].overflow = chOF[-2]
        self.chs['C'].overflow = chOF[-3]
        self.chs['D'].overflow = chOF[-4]
        channels = ['A','B','C','D']
        if reportOverflow:
            for i in range(4):
                if chOF[-(i+1)]:
                    print('channel {0} overflow'.format(channels[i]))

        maxADC = c_int16()
        self.status["maximumValue"] = ps.ps5000aMaximumValue(self.chandle, byref(maxADC))
        assert_pico_ok(self.status["maximumValue"])

        for key in self.chs:
            if self.chs[key].enabled:
                self.chs[key].timeSignal = (np.array(self.chs[key].buffer) / maxADC.value) * \
                    self.psVoltageRange[self.chs[key].range] * 1e-3

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

    def clearMem(self):
        self.f = None
        self.t = None
        for channel in ['A', 'B', 'C', 'D']:
            self.chs[channel].clearMem()

    def getData(self, range = None):
        if not range:
            range = [0, len(self.chs['A'].timeSignal)]
        return np.stack((self.chs['A'].timeSignal[range[0]:range[1]], \
            self.chs['B'].timeSignal[range[0]:range[1]], \
            self.chs['C'].timeSignal[range[0]:range[1]], \
            self.chs['D'].timeSignal[range[0]:range[1]]), 1)

    def getfs(self):
        return 1e9/self.timeInternalns.value

    def getfreq(self):
        return self.f

    def getConfigureInfo(self):
        re = {
            'channels':[self.chs['A'].getConfigureInfo(),\
                        self.chs['B'].getConfigureInfo(),\
                        self.chs['C'].getConfigureInfo(),\
                        self.chs['D'].getConfigureInfo()],
            'resolution':list(ps.PS5000A_DEVICE_RESOLUTION.keys())[list(ps.PS5000A_DEVICE_RESOLUTION.values()).index(self.resolution)],
            'timebase':self.timebase,
            'samplingIntervals':{'value':self.timeInternalns.value, 'unit':'ns'},
            'nSamples':self.maxSamples
        }
        return re

    def AWG_clock(self):
        f = 10e6
        self.status["SetSigGenBuiltInV2"] = ps.ps5000aSetSigGenBuiltInV2(\
            self.chandle, c_int32(0), c_uint32(int(1e6)), 1, f, f,\
            0, 0, 0, 0, c_uint32(0), c_uint32(0), 0, 0, c_int16(0))
        assert_pico_ok(self.status["SetSigGenBuiltInV2"])

def main(argv):
    flag = 0

if __name__ == "__main__":
    main(sys.argv[1:])
