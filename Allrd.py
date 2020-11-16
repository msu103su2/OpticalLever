import lib.redpitaya_scpi as scpi
import lib.Switch as Switch
import lib.AgilentNA as AG
import pyvisa as py
import numpy as np
import time
from scipy.signal import find_peaks

rp = scpi.scpi('192.168.137.209')
zx80 = Switch.zx80(rp)
NA = AG.AgilentNA()
rm = py.ResourceManager()
FG = rm.open_resource(rm.list_resources()[0])

PSD = np.loadtxt(open('D:\\20201109-0001.csv'),skiprows = 2, delimiter=",")
PSD = PSD[0:70000, :]
peaks, prop = find_peaks(PSD[:,1], prominence = 5)
locs = PSD[peaks, 0]
locs = locs[3:]

for loc in locs:
    zx80.RF2Off()
    zx80.RF1Off()
    NA.inst.write('BW 30')
    NA.inst.write('POIN 401')
    NA.inst.write('POWE -50')
    time.sleep(0.5)
    NA.center([loc*1e6, 2000])
    noiseFloor, noiseFloorStd = NA.getNoiseFloor([zx80])
    judge = True
    baseAMP = 0.001
    basePOWE = -50
    count = 0;
    zx80.RF2On()
    while judge and count <3:
        NA.inst.write('POWE '+str(basePOWE+10*count))
        NA.inst.write('CLES')
        NA.inst.write('SING')
        NA.auto()
        while not NA.IsFinish():
            time.sleep(2)
        PSD = NA.getTrace()
        peaks, _ = find_peaks(PSD, distance = 400)
        PowerCompensation = min(PSD) - (noiseFloor - 3*noiseFloorStd)
        judge = (PSD[peaks[0]] - noiseFloor - PowerCompensation)<5*noiseFloorStd
        count += 1

    if count < 2:
        NA.pkt()
        NA.auto()
        NA.ringdown([baseAMP*10**count, zx80, FG])
