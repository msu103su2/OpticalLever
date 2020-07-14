import sys, argparse
import pyvisa
import csv
import os
import numpy as np
import time
import matplotlib.pyplot as plt
from datetime import date


rm = pyvisa.ResourceManager()
inst = rm.open_resource('GPIB0::17::INSTR')
inst.write_termination = '\r'

def save():
    directory = r'Z:\data\optical lever project\2020_06_15'+'\\'
    inst.write('FORM3')
    y = inst.query_binary_values('OUTPDTRC?',datatype='d', is_big_endian=True)
    SWET = float(inst.query('SWET?'))
    POIN = int(inst.query('POIN?'))
    times = list(np.arange(0,SWET,SWET/POIN))
    pw = y[::2]

    rows = []
    fileds = ['Time(s)','Power(dB)']
    for i in range(0,len(times)):
        rows.append([times[i],pw[i]])
    centF =str(float(inst.query('CENT?'))/1e6)
    filename =  'f='+centF+'MHz.csv'
    with open(directory+filename, mode='w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(fileds)
        csvwriter.writerows(rows)

def savePSD(f, stepsize, end, DeviceSN):
    today = date.today()
    directory = r'Z:\data\optical lever project'+'\\'+today.strftime("%Y_%m_%d")+'\\'
    if not os.path.exists(directory):
        os.makedirs(directory)
    inst.write('FORM3')
    y = inst.query_binary_values('OUTPDTRC?',datatype='d', is_big_endian=True)
    POIN = int(inst.query('POIN?'))
    Frequencies = list(np.arange(f,min(f+stepsize,end),min(stepsize,end-f)/POIN))
    pw = y[::2]

    rows = []
    fileds = ['Frequencies','PSD(dB^2/Hz)']
    for i in range(0,len(Frequencies)):
        rows.append([Frequencies[i],pw[i]])
    filename = DeviceSN + '_' + '%4.3f'%(f/1e6) + '~'+ '%4.3f'%(min(f+stepsize,end)/1e6) + 'MHz.csv'
    with open(directory+filename, mode='w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(fileds)
        csvwriter.writerows(rows)

def rf(args):
    dB = args[0]
    inst.write('POWE '+str(int(dB)))

def rfu():
    dB = float(inst.query('POWE?'))+10
    inst.write('POWE '+str(int(dB)))

def rfd():
    dB = float(inst.query('POWE?'))-10
    inst.write('POWE '+str(int(dB)))

def to_1st_peak():
    inst.write('SEAM PEAK')

def next_peak():
    inst.write('SEANPK')

def rough_span():
    inst.write('SPAN 20000')

def fine_span():
    inst.write('SPAN 1000')

def span(values):
    inst.write('SPAN '+str(int(values[0])))

def check(values):
    start = values[0]
    end = values[1]
    if (len(values) < 3):
        power = -30
    else:
        power = values[2]
    inst.write('CLES')
    inst.write('POIN 401')
    inst.write('POWE '+str(power))
    inst.write('BW 30')
    inst.write('STAR '+str(int(start)))
    inst.write('STOP '+str(int(end)))
    inst.write('SING')
    inst.write('AUTO')

def zoom():
    inst.write('MKRCENT')
    inst.write('SPAN 20000')

def zoomzoom():
    #inst.write('SEANPK')
    inst.write('MKRCENT')
    inst.write('SPAN 1000')
    inst.write('CLES')
    inst.write('SING')

def ringdown():
    c = 20+float(inst.query('POWE?'))
    inst.write('POWE '+str(c))
    inst.write('POIN 401')
    inst.write('SPAN 0')
    inst.write('BW 100')
    inst.write('CLES')
    inst.write('SING')

def freqdomain():
    inst.write('CONT')
    inst.write('POWE -50')
    inst.write('BW 30')
    inst.write('SPAN 1000')

def auto():
    inst.write('AUTO')

def sweep(values, names):
    start = values[0]
    end = values[1]
    stepsize = values[2]
    DeviceSN = names[0]
    print(start,end,stepsize,DeviceSN)
    inst.write('BW 30')
    inst.write('POIN 801')
    for f in np.arange(start, end, stepsize):
        inst.write('CLES')
        inst.write('POWE -30')
        inst.write('STAR '+str(int(f)))
        inst.write('STOP '+str(int(min(f+stepsize,end))))
        inst.write('SING')
        inst.write('AUTO')
        while not IsFinish():
            time.sleep(2)
        savePSD(f, stepsize, end, DeviceSN)

def IsFinish():
    return int(inst.query('*STB?')) == 68

def showScreen():
    while not IsFinish():
        time.sleep(2)
    inst.write('FORM3')
    y = inst.query_binary_values('OUTPDTRC?',datatype='d', is_big_endian=True)
    POIN = int(inst.query('POIN?'))
    STAR = float(inst.query('STAR?'))
    STOP = float(inst.query('STOP?'))
    if (STAR == STOP):
        SWET = float(inst.query('SWET?'))
        x = list(np.arange(0, SWET, SWET/POIN))
    else:
        x = list(np.arange(STAR, STOP, (STOP-STAR)/POIN))
    plt.plot(x, y[::2])
    plt.show()

def ringdownRemote(values):
    dB = values[0]
    inst.write('SPAN 0')
    inst.write('POWE '+str(int(dB)))
    inst.write('POIN 801')
    inst.write('BW 100')
    time.sleep(2)
    inst.write('POWE -70')
    inst.write('CLES')
    inst.write('SING')
    showScreen()

def main(argv):
    parser = argparse.ArgumentParser(description = 'Convenient command for ring down measurements')
    parser.add_argument('command')
    parser.add_argument('-v','--val',type = float, nargs = '+', default = [0,0,-30])
    parser.add_argument('-n','--name', nargs = '+')
    args = parser.parse_args()

    if args.command == 's':
        save()
        print('saved')
    elif args.command == 'a':
        auto()
    elif args.command == 'rf':
        rf(args.val)
    elif args.command == 'rfu':
        rfu()
    elif args.command == 'rfd':
        rfd()
    elif args.command == 'c':
        check(args.val)
    elif args.command == 'z':
        zoom()
    elif args.command == 'zz':
        zoomzoom()
    elif args.command == 'rd':
        ringdown()
    elif args.command == 'rdr':
        ringdownRemote(args.val)
    elif args.command == 'p':
        to_1st_peak()
    elif args.command == 'np':
        next_peak()
    elif args.command == 'sp':
        span(args.val)
    elif args.command == 'fd':
        freqdomain()
    elif args.command == 'PSD':
        sweep(args.val, args.name)
    elif args.command == 'ss':
        showScreen()
    else:
        print("unknow command")

    rm.close()

if __name__ == "__main__":
    main(sys.argv[1:])
