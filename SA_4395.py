import lib.FuncGen33522b as FuncGen
import lib.AgilentNA as AG
import numpy as np
import time
import os
import matplotlib.pyplot as plt

if __name__ == '__main__':
    wd = r'Z:\data\optical lever project\NORCADA_NX53515C\73-CM_test'
    FG = FuncGen.FuncGen()
    NA = AG.AgilentNA(wd)

    bw = int(10)
    span = int(100)
    poin = int(span/bw)
    center = 176094

    fstart = center - 10
    fend = center + 10
    df = 0.1
    amp = 0.25
    dc = 0.365

    NA.center([center, span])
    NA.inst.write('POIN {poin:d}'.format(poin = poin))
    NA.inst.write('bw {bw:d}'.format(bw = bw))

    wd = wd+'\\'
    if not os.path.exists(wd):
        os.makedirs(wd)
    avg = 10
    freqs = np.arange(fstart, fend, df)
    freqs.tofile(wd+'PSDfreq'+'.bin', sep = '')

    PSDs = np.zeros((poin, avg))
    powers = np.zeros(len(freqs))
    j= 0



    for f in freqs:
        FG.Sine(f, amp, offset = dc)
        time.sleep(20)

        for i in range(avg):
            NA.inst.write('CLES')
            NA.inst.write('SING')
            while not NA.IsFinish():
                time.sleep(0.1)
            NA.inst.write('AUTO')
            PSDs[: , i] = np.array(NA.getTrace())

        powers[j] = 10*np.log10(np.sum(np.power(10, (PSDs/10)))/avg)
        filename = 'f={f:.2f}'.format(f = f)
        PSDs.tofile(wd+filename+'.bin', sep = '')
        print(f)
        j = j+1
        time.sleep(20)

    plt.plot(freqs, powers)
    plt.show()
