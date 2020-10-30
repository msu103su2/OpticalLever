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

class picoscope:

    def __init__(self, workingDirectory):
        self.chandle = c_int16()
        self.status = {}
        self.resolution = ps.PS5000A_DEVICE_RESOLUTION["PS5000A_DR_16BIT"]
        self.status["openunit"] = ps.ps5000aOpenUnit(byref(self.chandle), \
            None, self.resolution)
        assert_pico_ok(self.status["openunit"])
        self.defaultSetting()
        self.workingDirectory = workingDirectory + '\\'

    def defaultSetting(self):
        self.channel = ps.PS5000A_CHANNEL["PS5000A_CHANNEL_A"]
        self.chAcoupling_type = ps.PS5000A_COUPLING["PS5000A_AC"]
        self.chARange = ps.PS5000A_RANGE["PS5000A_50MV"]
        self.status["setChA"] = ps.ps5000aSetChannel(self.chandle, self.channel, 1,\
            self.chAcoupling_type, self.chARange, 0) #enabled = 1, analogue offset = 0 V
        self.timebase = int(6000)
        self.timeInternalns = c_float()
        self.bufferSize = c_int32()
        self.maxSamples = 10000
        self.status["getTimebase2"] = ps.ps5000aGetTimebase2(self.chandle, \
            self.timebase, self.maxSamples, byref(self.timeInternalns),\
            byref(self.bufferSize), 0)
        assert_pico_ok(self.status["getTimebase2"])

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
            preTriggerSamples = 100
            self.status["runBlock"] = ps.ps5000aRunBlock(self.chandle, \
            preTriggerSamples, self.maxSamples, self.timebase, None, 0, None, None)
            assert_pico_ok(self.status["runBlock"])

            ready = c_int16(0)
            check = c_int16(0)
            while ready.value == check.value:
                self.status["isReady"] = ps.ps5000aIsReady(self.chandle, byref(ready))

            bufferA = (c_int16 * self.maxSamples)()
            source = ps.PS5000A_CHANNEL["PS5000A_CHANNEL_A"]
            self.status["setDataBufferA"] = ps.ps5000aSetDataBuffer(self.chandle, \
                source, byref(bufferA), self.maxSamples, 0, 0)
            assert_pico_ok(self.status["setDataBufferA"])

            overflow = c_int16()
            cmaxSamples = c_uint32(self.maxSamples)
            self.status["getValues"] = ps.ps5000aGetValues(self.chandle, 0 , \
            byref(cmaxSamples), 0, 0, 0, byref(overflow))
            assert_pico_ok(self.status["getValues"])

            maxADC = c_int16()
            self.status["maximumValue"] = ps.ps5000aMaximumValue(self.chandle, byref(maxADC))
            assert_pico_ok(self.status["maximumValue"])

            timeSignal = adc2mV(bufferA, self.chARange, maxADC)
            timeSignal = [x * 1e-3 for x in timeSignal]
            f, newPSD = signal.periodogram(timeSignal, 1 / (self.timeInternalns.value / 1e9), \
                window = signal.get_window('blackman', len(timeSignal)))
            if i == 0:
                PSD = np.array(newPSD)
            else:
                PSD = PSD + newPSD
            energy = [x * x/50*self.timeInternalns.value/1e9 for x in timeSignal]
            print(sum(energy))
            energy = [x/50 for x in newPSD]
            print(sum(energy))

        PSD = PSD/nmbOfAvg
        PSD = PSD/(self.timeInternalns.value/1e9*len(timeSignal))
        PSD = 10 * np.log10(20 * PSD)
        self.memory = (f, PSD)

    def save(self, filename):
        if not os.path.exists(self.workingDirectory):
            os.makedirs(self.workingDirectory)
        rows = []
        fileds = ['f(Hz)','PSD(dBm)']
        for i in range(0,len(self.memory[1])):
            rows.append([self.memory[0][i],self.memory[1][i]])
        with open(self.workingDirectory+filename+'.csv', mode='w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(fileds)
            csvwriter.writerows(rows)

    def plot(self):
        plt.plot(self.memory[0],self.memory[1])
        plt.xlabel('f(Hz)')
        plt.ylabel('Power(dBm)')
        plt.show()

    def close(self):
        self.status["stop"] = ps.ps5000aStop(self.chandle)
        assert_pico_ok(self.status["stop"])

        self.status["close"]=ps.ps5000aCloseUnit(self.chandle)
        assert_pico_ok(self.status["close"])

def main(argv):
    flag = 0

if __name__ == "__main__":
    main(sys.argv[1:])
