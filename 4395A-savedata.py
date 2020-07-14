import pyvisa
import csv
import numpy as np
directory = r'Z:\data\optical lever project\2020_06_15'+'\\'

rm = pyvisa.ResourceManager()
HP4395a = rm.open_resource('GPIB0::17::INSTR')
HP4395a.write_termination = '\r'
HP4395a.write('FORM3')
y = HP4395a.query_binary_values('OUTPDTRC?',datatype='d', is_big_endian=True)
SWET = float(HP4395a.query('SWET?'))
times = list(np.arange(0,SWET,SWET/200))
pw = y[::2]

rows = []
fileds = ['Time(s)','Power(dB)']
for i in range(0,len(times)):
    rows.append([times[i],pw[i]])
centF =str(float(HP4395a.query('CENT?'))/1e6)
filename =  'f='+centF+'MHz.csv'
with open(directory+filename, mode='w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(fileds)
    csvwriter.writerows(rows)
