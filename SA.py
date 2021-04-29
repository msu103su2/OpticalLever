import lib.ps5000A as ps
import lib.FuncGen33522b as FuncGen
import numpy as np
import time
import os

if __name__ == '__main__':
    wd = r'Z:\data\optical lever project\NORCADA_NX53515C\71-CM_test'
    ps5000a = ps.picoscope(wd)
    FG = FuncGen.FuncGen()

    ps5000a.defaultSetting()
    ps5000a.AC('A')
    ps5000a.DC('C')
    ps5000a.DC('D')

    ps5000a.AutoRange('A')
    ps5000a.AutoRange('D')
    ps5000a.AutoRange('C')

    ps5000a.ChangeRangeTo('D', 1000)
    ps5000a.ChangeRangeTo('A', 1000)

    ps5000a.configurePSD(0.1, 2e6)

    fstart = 176185
    fend = 176195
    df = 0.1
    amp = 0.25
    dc = 0.365

    wd = wd+'\\'
    if not os.path.exists(wd):
        os.makedirs(wd)

    ps5000a.getTimeSignal()
    ps5000a.PSDfromTS(ps5000a.chs['A'].timeSignal, ps5000a.getfs())
    nprows = np.array(ps5000a.getfreq())
    highCutI = np.argmin(np.absolute(nprows-500e3))
    nprows = nprows[0:highCutI]
    filename = 'PSDfreq'
    nprows.tofile(wd+filename+'.bin', sep = '')

    cPowers = np.zeros(len(np.arange(fstart, fend, df)))
    i = 0

    for f in np.arange(fstart, fend, df):
        FG.Sine(f, amp, offset = dc)
        time.sleep(20)
        ps5000a.getTimeSignal()
        PSD = ps5000a.PSDfromTS(ps5000a.chs['A'].timeSignal, ps5000a.getfs())
        PSD = np.array(PSD)
        PSD = PSD[0:highCutI]
        cPowers[i] = np.mean(ps5000a.chs['C'].timeSignal)
        filename = 'f={f:.2f}'.format(f = f)
        PSD.tofile(wd+filename+'.bin', sep = '')
        print(f)
        time.sleep(20)
        i = i+1

    cPowers.tofile(wd+'cPowers'+'.bin', sep = '')
