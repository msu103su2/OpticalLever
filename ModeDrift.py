import lib.FuncGen33522b as FuncGen
import lib.AgilentNA as AG
import numpy as np
import time
import os
import matplotlib.pyplot as plt

if __name__ == '__main__':
    wd = r'Z:\data\optical lever project\NORCADA_NX53515C\75-drift_test'
    NA = AG.AgilentNA(wd)

    bw = int(2)
    span = int(200)
    poin = int(span/bw+1)
    center = 176094

    freqs = np.arange(center-span/2, center+span/2+1, bw)

    NA.center([center, span])
    NA.inst.write('POIN {poin:d}'.format(poin = poin))
    NA.inst.write('bw {bw:d}'.format(bw = bw))
    NA.rf([-50])

    wd = wd+'\\'
    if not os.path.exists(wd):
        os.makedirs(wd)

    PSDs = np.zeros((1, poin))
    powers = np.zeros(len(freqs))
    start = time.time()
    hour_passed = (time.time() - start)/3600
    f = []
    time_stamp = []

    while hour_passed < 12:
        NA.inst.write('CLES')
        NA.inst.write('SING')
        while not NA.IsFinish():
            time.sleep(0.1)
        NA.inst.write('AUTO')
        PSDs = np.append(PSDs, [np.array(NA.getTrace())], 0)
        hour_passed = (time.time() - start)/3600
        print(hour_passed)
        f.append(freqs[np.argmax(PSDs[-1,:])])
        time_stamp.append(hour_passed)

    PSDs = np.delete(PSDs, 0, 0)
    PSDs.tofile(wd+'PSDs'+'.bin', sep = '')
    f = np.array(f)
    time_stamp = np.array(time_stamp)
    f.tofile(wd+'fs'+'.bin', sep = '')
    time_stamp.tofile(wd+'time_stamp'+'.bin', sep = '')

    plt.plot(time_stamp, f)
    plt.show()
