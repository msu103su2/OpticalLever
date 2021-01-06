import lib.redpitaya_scpi as scpi
import lib.Switch as Switch
import lib.AgilentNA as AG
import lib.ps5000A as ps
import lib.splitPair as sp
import lib.FuncGen33522b as FuncGen
import pyvisa as py
import numpy as np
import time
wd = r'Z:\data\optical lever project\NORCADA_NX53515C\18-SNR'
#wd = r'C:\data'
rpip = '192.168.137.14'
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
ps5000a.ChangeRangeTo('A', 5000)
ps5000a.ChangeRangeTo('B', 5000)
ps5000a.ChangeRangeTo('C', 5000)

for i in range(30):
    ps5000a.DC('A')
    ps5000a.DC('B')
    ps5000a.DC('C')
    ps5000a.AutoRange('A')
    ps5000a.AutoRange('B')
    ps5000a.AutoRange('C')

    ps5000a.getTimeSignal();
    data = ps5000a.chs['A'].timeSignal - ps5000a.chs['B'].timeSignal
    unbalance = np.mean(data)
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
        ps5000a.AutoRange('A')
        ps5000a.AutoRange('B')
        ps5000a.getTimeSignal();
        data = ps5000a.chs['A'].timeSignal - ps5000a.chs['B'].timeSignal
        newUnbalance = np.mean(data)
        std = np.std(data)
        print(DS.getQuasi_gapsize(),' ',newUnbalance)
        print(np.mean(ps5000a.chs['A'].timeSignal),' ',np.mean(ps5000a.chs['B'].timeSignal), ' ',np.mean(ps5000a.chs['C'].timeSignal))
        if newUnbalance*unbalance<0:
            flag = not flag
        unbalance = newUnbalance

    ps5000a.AC('A')
    ps5000a.AC('B')
    ps5000a.AC('C')
    ps5000a.AutoRange('A')
    ps5000a.AutoRange('B')
    ps5000a.AutoRange('C')
    ps5000a.chs['A'].range += 1
    ps5000a.chs['B'].range += 1
    ps5000a.chs['C'].range += 1
    ps5000a.chs['A'].set()
    ps5000a.chs['B'].set()
    ps5000a.chs['C'].set()

    for count in range(1):
        ps5000a.getPSD(8.515, 4e6, avg)
        #ps5000a.saveTS('TS_GS='+str(DS.getQuasi_gapsize())+'_count='+str(count))
        ps5000a.savePSD('PSD_GS='+str(DS.getQuasi_gapsize())+'_count='+str(count))

    '''j = 0;
    while j<avg:
        ps5000a.getPSD(8.515, 4e6, 1)
        lowfPower = np.mean(ps5000a.memory[1][lowfRange[0]:lowfRange[1]])
        if lowfPower < upperBound:
            name = 'j='+str(j)+'_SNR_D='+str(DS.getQuasi_gapsize())
            ps5000a.save(name)
            j = j+1;
            print('accepted and j = ',str(j),', power=',str(lowfPower))
        else:
            print('rejected',', power=',str(lowfPower))'''
    DS.OpenGapBy(0.03)

DS.A.MoveTo(Ax)
DS.B.MoveTo(Bx)
