import lib.redpitaya_scpi as scpi
import lib.Switch as Switch
import lib.AgilentNA as AG
import lib.ps5000A as ps
import lib.splitPair as sp
import lib.FuncGen33522b as FuncGen
import pyvisa as py
import numpy as np
import time
wd = r'Z:\data\optical lever project\NORCADA_NX53515C\16-SNR'
rpip = '192.168.137.222'
rm = py.ResourceManager()

FG = FuncGen.FuncGen()
NA = AG.AgilentNA(wd)
RP = scpi.scpi(rpip)
DS = sp.splitPair()
zx80 = Switch.zx80(RP)

Ax = DS.A.Position();
Bx = DS.B.Position();

avg = 100;
lowfRange = [200, 400]
ps5000a = ps.picoscope(wd)
ps5000a.defaultSetting()
ps5000a.AC()
ps5000a.ChangeRangeTo(2000)
ps5000a.getPSD(8.515, 4e6, 1)
ps5000a.DC()
ps5000a.AutoRange()
for i in range(50):
    ps5000a.DC()
    ps5000a.AutoRange()
    data = ps5000a.getTimeSignal();
    unbalance = -np.mean(data)
    std = np.std(data)

    flag = True
    if abs(unbalance)>0.5:
        step = 0.005
    else:
        step = 0.002
    while abs(unbalance)>0.3*std:
        if unbalance>0:
            if flag:
                DS.B.MoveBy(-step)
            else:
                DS.A.MoveBy(step)
        else:
            if not flag:
                DS.B.MoveBy(step)
            else:
                DS.A.MoveBy(-step)
        ps5000a.AutoRange()
        data = ps5000a.getTimeSignal();
        newUnbalance = -np.mean(data)
        std = np.std(data)
        print(DS.getQuasi_gapsize(),' ',newUnbalance)
        if newUnbalance*unbalance<0:
            flag = not flag
        unbalance = newUnbalance

    ps5000a.AC()
    ps5000a.ChangeRangeTo(2000)
    lowfPower = []
    count = 0;
    while count<20:
        count = count+1;
        ps5000a.getPSD(8.515, 4e6, 1)
        lowfPower.append(np.mean(ps5000a.memory[1][lowfRange[0]:lowfRange[1]]))
    upperBound = np.mean(lowfPower) - np.std(lowfPower)
    print(upperBound)

    j = 0;
    while j<avg:
        ps5000a.getPSD(8.515, 4e6, 1)
        lowfPower = np.mean(ps5000a.memory[1][lowfRange[0]:lowfRange[1]])
        if lowfPower < upperBound:
            name = 'j='+str(j)+'_SNR_D='+str(DS.getQuasi_gapsize())
            ps5000a.save(name)
            j = j+1;
            print('accepted and j = ',str(j),', power=',str(lowfPower))
        else:
            print('rejected',', power=',str(lowfPower))
    DS.OpenGapBy(0.03)

DS.A.MoveTo(Ax)
DS.B.MoveTo(Bx)
ps5000a.DC()
ps5000a.AutoRange()
data = ps5000a.getTimeSignal();
unbalance = -np.mean(data)
std = np.std(data)
flag = True
if abs(unbalance)>0.5:
    step = 0.005
else:
    step = 0.002
while abs(unbalance)>3*std:
    if unbalance>0:
        if flag:
            DS.B.MoveBy(-step)
        else:
            DS.A.MoveBy(step)
    else:
        if not flag:
            DS.B.MoveBy(step)
        else:
            DS.A.MoveBy(-step)
    data = ps5000a.getTimeSignal();
    newUnbalance = -np.mean(data)
    std = np.std(data)
    if newUnbalance*unbalance<0:
        flag = not flag
    unbalance = newUnbalance
