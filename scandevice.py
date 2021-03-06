import lib.AgilentNA as AgilentNA
import lib.TLPM as TLPM
import sys
sys.path.append(r'C:\Program Files\Newport\Piezo Motion Control\Newport AG-UC2-UC8 Applet\Bin')
import clr
clr.AddReference('AgilisCmdLib')
from Newport.AgilisCmdLib import *
from ctypes import cdll, c_long, c_ulong, c_uint32, byref, create_string_buffer, c_bool, c_char_p, c_int, c_int16, c_int32, c_double, sizeof, c_voidp, c_float
import numpy as np
import time
from scipy.signal import find_peaks
from scipy import signal
import glob
import csv
from math import floor, ceil
import mysql.connector
from picosdk.ps5000a import ps5000a as ps
from picosdk.functions import assert_pico_ok, adc2mV
import matplotlib.pyplot as plt

workingDirectory = r'Z:\data\optical lever project\Die 000_21\03_scanPSDs'
psVoltageRange = {
    0 : 10, 1 : 20, 2 : 50, 3 : 100, 4 : 200, 5 : 500, 6 : 1000, 7 : 2000, 8 : 5000, 9 : 10000, 10 : 20000
}

#functions to simplify codes
def searchDevice(searchDirection):
    #start from void
    newLog = np.array([px[-1]])
    autoRange()
    minReflected = readVoltage()
    newLog = np.append(newLog, [[newLog[-1][0], minReflected]], axis = 0)
    print('step = %.0f; BPD voltage = %.2fmV;'%(newLog[-1][0],newLog[-1][1]))

    while (abs(newLog[-1][0] - newLog[0][0])  < 3000 and abs((newLog[-1][1] - minReflected)/minReflected) < 10):
        step = 5 * searchDirection
        UC2.PR(controllerAddress, step, err)
        time.sleep(1)
        autoRange()
        newLog = np.append(newLog, [[step+newLog[-1][0], readVoltage()]], axis = 0)
        print('step = %.0f; BPD voltage = %.2fmV;'%(newLog[-1][0],newLog[-1][1]))

    return newLog

def scanDevice(deviceCount):
    newLog = moveAround()
    BW = 30
    POIN = 801
    POWE = -30
    fstart = 1e5
    fend = 4.9e6
    fstep = 2.4e4
    filename = ('000_21_%02i'%(deviceCount))
    HP4395a.sweep([fstart, fend, fstep, BW, POIN, POWE], [filename])
    print('step = %.0f; BPD voltage = %.2fmV; One PSD obtained, device = %02i'%(newLog[-1][0],newLog[-1][1], deviceCount))
    return newLog
def travel(travelDirection):
    # start from last device
    firstJump = 550 * travelDirection
    UC2.PR(controllerAddress, firstJump, err)
    time.sleep(5)
    autoRange()
    newLog = np.array([[firstJump+px[-1][0], readVoltage()]])
    print('step = %.0f; BPD voltage = %.2fmV;'%(newLog[-1][0],newLog[-1][1]))
    return newLog

def hasnext(i, px):
    return i <= 17 and abs(px[-1][0]) < 14000;

def CollectModes(filenameRegx, directory = workingDirectory):
    pksLocation = []
    filelist = glob.glob(directory+'\\'+filenameRegx)
    for file in filelist:
        PSD = np.loadtxt(open(file),skiprows = 1, delimiter=",")
        mean = np.mean(PSD[:,1])
        std = np.std(PSD[:,1])
        pksIndex,_= find_peaks(PSD[:,1], height = mean + 3*std, distance = 33)
        pksLocation.extend(PSD[pksIndex,0])
    return pksLocation

def fineOnModes(pksLocation, deviceCount):
    BW = 2;
    POIN = 201;
    POWE = -50;
    fSTEP = 400;
    for pkF in pksLocation:
        HP4395a.sweep([int(pkF) - fSTEP/2, int(pkF) + fSTEP/2, fSTEP, BW, POIN, POWE], ['%02iCF=%iHz'%(deviceCount,pkF)])
        #HP4395a.sweep([int(pkF) - 100, int(pkF) + 100, fSTEP, BW, POIN, POWE+20], ['%02iCF=%iHz_Oc'%(pkF)]) #try overdrive
        HP4395a.inst.write('POWE '+ str(int(POWE)))

def readAmplitude(binWidth, spectrumRange, checkRange, nmbOfAvg):

    spectrumRange = spectrumRange * 2
    requiredMeasureTime = 2 / binWidth
    requiredSamplingInterval = 1 / spectrumRange
    timebase = floor(requiredSamplingInterval * 62500000 + 3)
    timeInternalns = c_float()
    returnedMaxSamples = c_int32()
    maxSamples = ceil(spectrumRange / binWidth)
    status["getTimebase2"] = ps.ps5000aGetTimebase2(chandle, timebase, maxSamples, byref(timeInternalns),byref(returnedMaxSamples), 0)
    assert_pico_ok(status["getTimebase2"])
    assert timeInternalns.value < requiredSamplingInterval * 1e9

    for i in range(nmbOfAvg):
        preTriggerSamples = 100
        status["runBlock"] = ps.ps5000aRunBlock(chandle, preTriggerSamples, maxSamples, timebase, None, 0, None, None)
        assert_pico_ok(status["runBlock"])

        ready = c_int16(0)
        check = c_int16(0)
        while ready.value == check.value:
            status["isReady"] = ps.ps5000aIsReady(chandle, byref(ready))

        bufferA = (c_int16 * maxSamples)()
        source = ps.PS5000A_CHANNEL["PS5000A_CHANNEL_A"]
        status["setDataBufferA"] = ps.ps5000aSetDataBuffer(chandle, source, byref(bufferA), maxSamples, 0, 0)
        assert_pico_ok(status["setDataBufferA"])

        overflow = c_int16()
        cmaxSamples = c_uint32(maxSamples)
        status["getValues"] = ps.ps5000aGetValues(chandle, 0 , byref(cmaxSamples), 0, 0, 0, byref(overflow))
        assert_pico_ok(status["getValues"])

        maxADC = c_int16()
        status["maximumValue"] = ps.ps5000aMaximumValue(chandle, byref(maxADC))
        assert_pico_ok(status["maximumValue"])

        timeSignal = adc2mV(bufferA, chARange, maxADC)
        f, newPSD = signal.periodogram(timeSignal, 1 / (timeInternalns.value / 1e9), window = signal.get_window('blackman', len(timeSignal)))
        startIndex = f.searchsorted(checkRange[0]) - 1
        endIndex = f.searchsorted(checkRange[1]) + 1
        f = f[startIndex:endIndex]
        newPSD = newPSD[startIndex:endIndex]
        if i == 0:
            PSD = np.array(newPSD)
        else:
            PSD = PSD + newPSD

    PSD = PSD/nmbOfAvg;

    PSD = 10 * np.log10(10 * PSD)
    PSD = PSD - PSD.mean()
    Index, _ = find_peaks(PSD, distance = PSD.size)
    return PSD[Index]

def readVoltage():

    maxSamples = 10000
    preTriggerSamples = 100
    status["runBlock"] = ps.ps5000aRunBlock(chandle, preTriggerSamples, maxSamples, timebase, None, 0, None, None)
    assert_pico_ok(status["runBlock"])

    ready = c_int16(0)
    check = c_int16(0)
    while ready.value == check.value:
        status["isReady"] = ps.ps5000aIsReady(chandle, byref(ready))

    bufferA = (c_int16 * 2)()
    source = ps.PS5000A_CHANNEL["PS5000A_CHANNEL_A"]
    downSampleTatioMode = ps.PS5000A_RATIO_MODE["PS5000A_RATIO_MODE_AVERAGE"]
    status["setDataBufferA"] = ps.ps5000aSetDataBuffer(chandle, source, byref(bufferA), 2, 0, downSampleTatioMode)
    assert_pico_ok(status["setDataBufferA"])

    maxDownSampleRatio = c_uint32()
    status["dwonSample"] = ps.ps5000aGetMaxDownSampleRatio(chandle, maxSamples, byref(maxDownSampleRatio), downSampleTatioMode, 0)
    assert_pico_ok(status["dwonSample"])

    overflow = c_int16()
    cmaxSamples = c_uint32(maxSamples)
    status["getValues"] = ps.ps5000aGetValues(chandle, 0 , byref(cmaxSamples), maxDownSampleRatio, downSampleTatioMode, 0, byref(overflow))
    assert_pico_ok(status["getValues"])

    maxADC = c_int16()
    status["maximumValue"] = ps.ps5000aMaximumValue(chandle, byref(maxADC))
    assert_pico_ok(status["maximumValue"])

    adc2mVChA = adc2mV(bufferA, chARange, maxADC)
    avg = adc2mVChA[0]
    return avg

def autoRange():
    global chARange
    maxSamples = 100
    while chARange < max(psVoltageRange.keys()):
        overrange = False
        preTriggerSamples = 100
        status["runBlock"] = ps.ps5000aRunBlock(chandle, preTriggerSamples, maxSamples, timebase, None, 0, None, None)
        assert_pico_ok(status["runBlock"])

        ready = c_int16(0)
        check = c_int16(0)
        while ready.value == check.value:
            status["isReady"] = ps.ps5000aIsReady(chandle, byref(ready))

        bufferA = (c_int16 * maxSamples)()
        source = ps.PS5000A_CHANNEL["PS5000A_CHANNEL_A"]
        status["setDataBufferA"] = ps.ps5000aSetDataBuffer(chandle, source, byref(bufferA), maxSamples, 0, 0)
        assert_pico_ok(status["setDataBufferA"])

        overflow = c_int16()
        cmaxSamples = c_uint32(maxSamples)
        status["getValues"] = ps.ps5000aGetValues(chandle, 0 , byref(cmaxSamples), 0, 0, 0, byref(overflow))
        assert_pico_ok(status["getValues"])

        maxADC = c_int16()
        status["maximumValue"] = ps.ps5000aMaximumValue(chandle, byref(maxADC))
        assert_pico_ok(status["maximumValue"])

        if max(map(abs, adc2mV(bufferA, chARange, maxADC))) == psVoltageRange[chARange]:
            overrange = True

        if overrange:
            chARange += 1
            status["setChA"] = ps.ps5000aSetChannel(chandle, channel, 1, coupling_type, chARange, 0) #enabled = 1, analogue offset = 0 V
            assert_pico_ok(status["setChA"])
        else:
            break

    while chARange > min(psVoltageRange.keys()):
        toosmall = False
        status["setChA"] = ps.ps5000aSetChannel(chandle, channel, 1, coupling_type, chARange - 1, 0) #enabled = 1, analogue offset = 0 V
        assert_pico_ok(status["setChA"])
        preTriggerSamples = 100
        status["runBlock"] = ps.ps5000aRunBlock(chandle, preTriggerSamples, maxSamples, timebase, None, 0, None, None)
        assert_pico_ok(status["runBlock"])

        ready = c_int16(0)
        check = c_int16(0)
        while ready.value == check.value:
            status["isReady"] = ps.ps5000aIsReady(chandle, byref(ready))

        bufferA = (c_int16 * maxSamples)()
        source = ps.PS5000A_CHANNEL["PS5000A_CHANNEL_A"]
        status["setDataBufferA"] = ps.ps5000aSetDataBuffer(chandle, source, byref(bufferA), maxSamples, 0, 0)
        assert_pico_ok(status["setDataBufferA"])

        overflow = c_int16()
        cmaxSamples = c_uint32(maxSamples)
        status["getValues"] = ps.ps5000aGetValues(chandle, 0 , byref(cmaxSamples), 0, 0, 0, byref(overflow))
        assert_pico_ok(status["getValues"])

        maxADC = c_int16()
        status["maximumValue"] = ps.ps5000aMaximumValue(chandle, byref(maxADC))
        assert_pico_ok(status["maximumValue"])

        if max(map(abs, adc2mV(bufferA, chARange, maxADC))) < psVoltageRange[chARange]:
            toosmall = True
        if toosmall:
            chARange = chARange - 1
        else:
            status["setChA"] = ps.ps5000aSetChannel(chandle, channel, 1, coupling_type, chARange, 0) #enabled = 1, analogue offset = 0 V
            assert_pico_ok(status["setChA"])
            break

def moveAround():
    autoRange()
    avg = readVoltage()
    totalTrial = 0
    newLog = np.array([[px[-1][0], avg]])
    while abs(avg + 0.65) > 1:
        if avg + 0.65 > 0 :
            UC2.PR(controllerAddress, 1, err)
            time.sleep(0.5)
            autoRange()
            avg = readVoltage()
            newLog = np.append(newLog, [[1 + newLog[-1][0], avg]], axis = 0)
            print('step = %.0f; BPD voltage = %.2fmV;'%(newLog[-1][0],newLog[-1][1]))
        else:
            UC2.PR(controllerAddress, -1, err)
            time.sleep(0.5)
            autoRange()
            avg = readVoltage()
            newLog = np.append(newLog, [[-1 + newLog[-1][0], avg]], axis = 0)
            print('step = %.0f; BPD voltage = %.2fmV;'%(newLog[-1][0],newLog[-1][1]))
        totalTrial += 1

    autoRange()
    avg = readVoltage()
    totalTrial = 0
    while readAmplitude(30, 15e6, [50e3, 200e3], 20) > 7:
        if avg + 0.65 > 0 :
            UC2.PR(controllerAddress, 1, err)
            time.sleep(1)
            autoRange()
            avg = readVoltage()
            newLog = np.append(newLog, [[1 + newLog[-1][0], avg]], axis = 0)
            print('step = %.0f; BPD voltage = %.2fmV;'%(newLog[-1][0],newLog[-1][1]))
        else:
            UC2.PR(controllerAddress, -1, err)
            time.sleep(1)
            autoRange()
            avg = readVoltage()
            newLog = np.append(newLog, [[-1 + newLog[-1][0], avg]], axis = 0)
            print('step = %.0f; BPD voltage = %.2fmV;'%(newLog[-1][0],newLog[-1][1]))
        totalTrial += 1
    return newLog

#initialize PM400
PM400 = TLPM.TLPM()
TLPMdeviceCount = c_uint32()
PM400.findRsrc(byref(TLPMdeviceCount))

TLPMresourceName = create_string_buffer(1024)
PM400.getRsrcName(c_int(0), TLPMresourceName)

PM400.open(TLPMresourceName, c_bool(True), c_bool(True))
PM400.setWavelength(c_double(1064))
PM400.setPowerAutoRange(c_bool(True))
PM400.setPowerUnit(c_int16(0))
time.sleep(1)

# initialize AG-UC2
UC2 = AgilisCmds()
DevicePorts = UC2.GetDevices()
UC2.OpenInstrument(DevicePorts[1])
controllerAddress = int(1);
err = ''


# initialize HP4395a
HP4395a = AgilentNA.AgilentNA(workingDirectory) # initialization of 4395a to be added later
#PM400.setAvgCnt(c_int16(10))
power = c_double()
PM400.measPower(byref(power))
time.sleep(1)
px = np.array([[0, power.value]])

#initialize picoscope
chandle = c_int16()
status = {}
resolution = ps.PS5000A_DEVICE_RESOLUTION["PS5000A_DR_16BIT"]
status["openunit"] = ps.ps5000aOpenUnit(byref(chandle), None, resolution)
channel = ps.PS5000A_CHANNEL["PS5000A_CHANNEL_A"]
coupling_type = ps.PS5000A_COUPLING["PS5000A_DC"]
chARange = ps.PS5000A_RANGE["PS5000A_50MV"]
status["setChA"] = ps.ps5000aSetChannel(chandle, channel, 1, coupling_type, chARange, 0) #enabled = 1, analogue offset = 0 V
assert_pico_ok(status["setChA"])

timebase = int(6000)
timeInternalns = c_float()
returnedMaxSamples = c_int32()
maxSamples = 10000
status["getTimebase2"] = ps.ps5000aGetTimebase2(chandle, timebase, maxSamples, byref(timeInternalns),byref(returnedMaxSamples), 0)
assert_pico_ok(status["getTimebase2"])


#connect to Database
DeviceDB = mysql.connector.connect(user='shanhao', password = 'SloanGW@138', host = '136.142.206.151', database='DeviceDB',port = '3307')


i = 16;
while (hasnext(i, px)):
    newLog = searchDevice(-1)
    px = np.append(px, newLog, axis = 0)
    newLog = scanDevice(i)
    px = np.append(px, newLog, axis = 0)
    #fineOnModes(CollectModes('%s_%02i_*'%('000_21', i)), i)
    newLog = travel(-1)
    px = np.append(px, newLog, axis = 0)
    i = i + 1;

with open(workingDirectory+'\\'+'StepLog.csv', mode='w', newline='') as csvfile:
    fileds = ['step','transmittedPower']
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(fileds)
    csvwriter.writerows(px)

HP4395a.close()
PM400.close()
UC2.CloseInstrument()

status["stop"] = ps.ps5000aStop(chandle)
assert_pico_ok(status["stop"])

status["close"]=ps.ps5000aCloseUnit(chandle)
assert_pico_ok(status["close"])
