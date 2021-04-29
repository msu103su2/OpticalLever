import lib.ps5000A as ps
import numpy as np
import matplotlib.pyplot as plt
wd = r'Z:\data\optical lever project\NORCADA_NX53515C\68-SNR'
ps5000a = ps.picoscope(wd)
ps5000a.defaultSetting()
ps5000a.DC('D')
ps5000a.AutoRange('D')
ps5000a.configurePSD(8.515, 4e6)
fs = ps5000a.getfs()
V = ps5000a.getTimeSignal(channel = 'D')
N = len(V)
T = N/fs
fftD = np.fft.fft(V)
if N % 2 == 0:
    f = np.arange(N/2+1)/T
    fftD = fftD[0:int(N/2)+1]
else:
    f = np.arange((N+1)/2)/T
    fftD = fftD[0:int((N+1)/2)]
phase = np.angle(fftD, deg=True)
amplitude = np.absolute(fftD)
plt.plot(f, phase, 'b-', f,amplitude,'r-')
plt.show()
