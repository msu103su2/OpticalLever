from picosdk.ps5000a import ps5000a as ps
from picosdk.functions import assert_pico_ok, adc2mV
from ctypes import cdll, c_long, c_ulong, c_uint32, byref, create_string_buffer, c_bool, c_char_p, c_int, c_int16, c_int32, c_double, sizeof, c_voidp, c_float

#initialize picoscope
chandle = c_int16()
status = {}
resolution = ps.PS5000A_DEVICE_RESOLUTION["PS5000A_DR_16BIT"]
status["openunit"] = ps.ps5000aOpenUnit(byref(chandle), None, resolution)
channel = ps.PS5000A_CHANNEL["PS5000A_CHANNEL_A"]
coupling_type = ps.PS5000A_COUPLING["PS5000A_DC"]
chARange = ps.PS5000A_RANGE["PS5000A_100MV"]
status["setChA"] = ps.ps5000aSetChannel(chandle, channel, 1, coupling_type, chARange, 0) #enabled = 1, analogue offset = 0 V
assert_pico_ok(status["setChA"])

timebase = int(9)
timeInternalns = c_float()
returnedMaxSamples = c_int32()
maxSamples = 10000
status["getTimebase2"] = ps.ps5000aGetTimebase2(chandle, timebase, maxSamples, byref(timeInternalns),byref(returnedMaxSamples), 0)
assert_pico_ok(status["getTimebase2"])

preTriggerSamples = 100
status["runBlock"] = ps.ps5000aRunBlock(chandle, preTriggerSamples, maxSamples, timebase, None, 0, None, None)
assert_pico_ok(status["runBlock"])

ready = c_int16(0)
check = c_int16(0)
while ready.value == check.value:
    status["isReady"] = ps.ps5000aIsReady(chandle, byref(ready))

bufferA = (c_int16 * maxSamples)()
source = ps.PS5000A_CHANNEL["PS5000A_CHANNEL_A"]
downSampleTatioMode = ps.PS5000A_RATIO_MODE["PS5000A_RATIO_MODE_AVERAGE"]
status["setDataBufferA"] = ps.ps5000aSetDataBuffer(chandle, source, byref(bufferA), maxSamples, 0, downSampleTatioMode)
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

status["stop"] = ps.ps5000aStop(chandle)
assert_pico_ok(status["stop"])

status["close"]=ps.ps5000aCloseUnit(chandle)
assert_pico_ok(status["close"])
