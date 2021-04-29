import lib.ps5000A as ps
import lib.splitPair as sp
import pyvisa as py
import numpy as np
from scipy.optimize import minimize
from scipy import signal
import time
import os

if __name__ == '__main__':
    steps = 1
    safeRangeA = [-1, 9.5]
    safeRangeB = [-1, 9.5]
    wd = r'Z:\data\optical lever project\NORCADA_NX53515C\76-BeamWaists'
    DS = sp.splitPair(safeRangeA , safeRangeB)

    Ax = DS.A.Position();
    Bx = DS.B.Position();

    ps5000a = ps.picoscope(wd)
    ps5000a.defaultSetting()
    ps5000a.AC('A')
    ps5000a.AC('B')
    ps5000a.DC('C')

    ps5000a.AutoRange('A')
    ps5000a.AutoRange('B')
    ps5000a.AutoRange('C')
    ps5000a.configurePSD(8.515, 4e6)

    for j in range(steps):
        # split mirror balancing process
        ps5000a.DC('A')
        ps5000a.DC('B')
        ps5000a.AutoRange('A')
        ps5000a.AutoRange('B')
        sign = -1

        ps5000a.getTimeSignal();
        data = ps5000a.getData()
        data = data[:,0] #- data[:,1]
        unbalance = sign*np.mean(data)
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
            data = ps5000a.getData()
            diff = sign*data[:,0] #- data[:,1]
            newUnbalance = np.mean(diff)
            std = np.std(diff)
            print(DS.getQuasi_gapsize(),' ',newUnbalance)
            print(np.mean(data[:,0]),' ',np.mean(data[:,1]), ' ',np.mean(data[:,2]))
            if newUnbalance*unbalance<0:
                flag = not flag
            unbalance = newUnbalance
