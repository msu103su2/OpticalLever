import sys, argparse
import pyvisa
import csv
import os
import numpy as np
import time
import matplotlib.pyplot as plt
import lib.redpitaya_scpi as scpi
from datetime import date

class AgilentNA:

    def __init__(self, workingDirectory = r'D:\temp'):
        self.rm = pyvisa.ResourceManager()
        self.inst = self.rm.open_resource('GPIB1::17::INSTR')
        self.inst.write_termination = '\r'
        self.inst.timeout = 20000
        self.workingDirectory = workingDirectory + '\\'
        self.inst.write('PRES')
        self.inst.write('NA')
        self.inst.write('MEAS R')
        self.inst.write('ESNB 1')
        self.inst.write('*SRE 4')
        self.inst.write('BW 30')
        self.inst.write('MKR ON')

    def savePSD(self, DeviceSN):
        self.inst.write('FORM3')
        if not os.path.exists(self.workingDirectory):
            os.makedirs(self.workingDirectory)

        y = self.inst.query_binary_values('OUTPDTRC?',datatype='d', is_big_endian=True)
        pw = y[::2]

        POIN = int(self.inst.query('POIN?'))
        STAR = float(self.inst.query('STAR?'))
        STOP = float(self.inst.query('STOP?'))
        if (STAR == STOP):
            SWET = float(self.inst.query('SWET?'))
            x = list(np.arange(0, SWET, SWET/POIN))
            xname = 'Time(s)'
            yname = 'Power(dB)'
            centF =str(float(self.inst.query('CENT?'))/1e6)
            filename =  DeviceSN + '_' + 'f='+centF+'MHz.csv'
        else:
            x = list(np.arange(STAR, STOP, (STOP-STAR)/POIN))
            xname = 'Frequencies(Hz)'
            yname = 'PSD(dB^2/Hz)'
            filename = DeviceSN + '_' + '%4.3f'%(STAR/1e6) + '~'+ '%4.3f'%(STOP/1e6) + 'MHz.csv'

        rows = []
        fileds = [xname,yname]
        for i in range(0,len(x)):
            rows.append([x[i],pw[i]])
        with open(self.workingDirectory+filename, mode='w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(fileds)
            csvwriter.writerows(rows)

    def getTrace(self):
        self.inst.write('FORM3')
        y = self.inst.query_binary_values('OUTPDTRC?',datatype='d', is_big_endian=True)
        pw = y[::2]
        return pw

    def rf(self, args):
        dB = args[0]
        self.inst.write('POWE '+str(int(dB)))

    def rfu(self):
        dB = float(self.inst.query('POWE?'))+10
        self.inst.write('POWE '+str(int(dB)))

    def rfd(self):
        dB = float(self.inst.query('POWE?'))-10
        self.inst.write('POWE '+str(int(dB)))

    def to_1st_peak(self):
        self.inst.write('SEAM PEAK')

    def next_peak(self):
        self.inst.write('SEANPK')

    def rough_span(self):
        self.inst.write('SPAN 20000')

    def fine_span(self):
        self.inst.write('SPAN 1000')

    def span(self, values):
        self.inst.write('SPAN '+str(int(values[0])))

    def check(self, values):
        start = values[0]
        end = values[1]
        if (len(values) < 3):
            power = -30
        else:
            power = values[2]
        self.inst.write('CLES')
        self.inst.write('POIN 401')
        self.inst.write('POWE '+str(power))
        self.inst.write('BW 30')
        self.inst.write('STAR '+str(int(start)))
        self.inst.write('STOP '+str(int(end)))
        self.inst.write('SING')
        self.inst.write('AUTO')

    def zoom(self):
        self.inst.write('MKRCENT')
        self.inst.write('SPAN 20000')

    def zoomzoom(self):
        #self.inst.write('SEANPK')
        self.inst.write('MKRCENT')
        self.inst.write('SPAN 1000')
        self.inst.write('CLES')
        self.inst.write('SING')

    def freqdomain(self):
        self.inst.write('CONT')
        self.inst.write('POWE -50')
        self.inst.write('BW 30')
        self.inst.write('SPAN 1000')

    def auto(self):
        self.inst.write('AUTO')

    def sweep(self, values, names):
        start = values[0]
        end = values[1]
        stepsize = values[2]
        BW = values[3]
        POIN = values[4]
        POWE = values[5]
        DeviceSN = names[0]
        print(start,end,stepsize,DeviceSN)
        self.inst.write('BW ' + str(int(BW)))
        self.inst.write('POIN '+ str(int(POIN)))
        for f in np.arange(start, end, stepsize):
            self.inst.write('POWE '+ str(int(POWE)))
            self.inst.write('STAR '+str(int(f)))
            self.inst.write('STOP '+str(int(min(f+stepsize,end))))
            self.inst.write('CLES')
            self.inst.write('SING')
            self.inst.write('AUTO')
            while not self.IsFinish():
                time.sleep(2)
            self.savePSD(DeviceSN)

    def IsFinish(self):
        return int(self.inst.query('*STB?')) == 68

    def showScreen(self):
        while not self.IsFinish():
            time.sleep(2)
        self.inst.write('FORM3')
        y = self.inst.query_binary_values('OUTPDTRC?',datatype='d', is_big_endian=True)
        POIN = int(self.inst.query('POIN?'))
        STAR = float(self.inst.query('STAR?'))
        STOP = float(self.inst.query('STOP?'))
        if (STAR == STOP):
            SWET = float(self.inst.query('SWET?'))
            x = list(np.arange(0, SWET, SWET/POIN))
        else:
            x = list(np.arange(STAR, STOP, (STOP-STAR)/POIN))
        plt.plot(x, y[::2])
        plt.show()

    def ringdown(self, values):
        amp = values[0]
        zx80 = values[1]
        FG = values[2]

        zx80.RF2Off()
        self.inst.write('SPAN 0')
        self.inst.write('POIN 801')
        self.inst.write('BW 100')
        noiseFloor,_ = self.getNoiseFloor([zx80])

        freq = str(float(self.inst.query('CENT?')))
        FG.Sine(freq, amp)
        zx80.RF1On()

        stable = False
        count = 0
        powerUpCount = 0
        while not stable and powerUpCount < 2:
            self.inst.write('CLES')
            self.inst.write('SING')
            self.auto()
            while not self.IsFinish():
                time.sleep(2)
            pw = self.getTrace()
            if np.std(pw) < 0.03*(np.mean(pw)-noiseFloor):
                stable = True
            count += 1
            if (count > 3):
                self.rfu()
                count = 0
                powerUpCount += 1

        if powerUpCount < 2:
            self.inst.write('CLES')
            self.inst.write('SING')
            time.sleep(2)
            zx80.RF1Off()
            FG33522b.write('OUTP OFF')
            self.auto()
            while not self.IsFinish():
                time.sleep(2)
            self.savePSD('rd')

    def close(self):
        self.rm.close()

    def center(self, values):
        self.inst.write('CENT '+str(int(values[0])))
        self.inst.write('SPAN '+str(int(values[1])))

    def pkt(self):
        self.inst.write('SEAM PEAK')
        self.inst.write('MKRCENT')
        self.inst.write('SPAN 0')

    def getNoiseFloor(self, values):
        zx80 = values[0]
        zx80.RF1Off()
        zx80.RF2Off()
        self.inst.write('CLES')
        self.inst.write('SING')
        self.auto()
        while not self.IsFinish():
            time.sleep(2)
        return np.mean(self.getTrace()), np.std(self.getTrace())

def main(argv):
    parser = argparse.ArgumentParser(description = 'Convenient command for ring down measurements')
    parser.add_argument('command')
    parser.add_argument('-v','--val',type = float, nargs = '+', default = [0,0,-30])
    parser.add_argument('-n','--name', nargs = '+')
    args = parser.parse_args()

    Agilent4935A = AgilentNA()

    if args.command == 's':
        Agilent4935A.save()
        print('saved')
    elif args.command == 'a':
        Agilent4935A.auto()
    elif args.command == 'rf':
        Agilent4935A.rf(args.val)
    elif args.command == 'rfu':
        Agilent4935A.rfu()
    elif args.command == 'rfd':
        Agilent4935A.rfd()
    elif args.command == 'c':
        Agilent4935A.check(args.val)
    elif args.command == 'z':
        Agilent4935A.zoom()
    elif args.command == 'zz':
        Agilent4935A.zoomzoom()
    elif args.command == 'rd':
        Agilent4935A.ringdown()
    elif args.command == 'rdr':
        Agilent4935A.ringdownRemote(args.val)
    elif args.command == 'p':
        Agilent4935A.to_1st_peak()
    elif args.command == 'np':
        Agilent4935A.next_peak()
    elif args.command == 'sp':
        Agilent4935A.span(args.val)
    elif args.command == 'fd':
        Agilent4935A.freqdomain()
    elif args.command == 'PSD':
        Agilent4935A.sweep(args.val, args.name)
    elif args.command == 'ss':
        Agilent4935A.showScreen()
    else:
        print("unknow command")

    Agilent4935A.close()

if __name__ == "__main__":
    main(sys.argv[1:])
