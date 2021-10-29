import pyvisa
import time
import re
import numpy as np
ReadyT = ['1TS000036','1TS000037','1TS000038']
Ready = ['1TS000032','1TS000033','1TS000034']
Disable = ['1TS00003D']
DisableT = ['1TS00003E']
Tracking = ['1TS000046']
class CONEX_CC:
    """docstring for CONEX_CC."""

    def __init__(self, deviceKey, safeRange, reset = True):
        self.rm = pyvisa.ResourceManager()
        self.inst = self.rm.open_resource(deviceKey, baud_rate = 921600,\
            data_bits = 8, write_termination = '\r\n', read_termination = '\r\n',\
            timeout = 2000)
        self.min_step = 0.0001
        self.safeRange = safeRange #in unit of mm, [min, max]

        if self.inst.query('1TS') not in ReadyT:
            reset = True

        if reset:
            self.inst.write('1RS')
            while self.inst.query('1TS')!='1TS00000A':
                time.sleep(1)
            self.inst.write('1OR')
            while self.inst.query('1TS') not in Ready:
                time.sleep(1)
            self.inst.write('1TK1')
            while self.inst.query('1TS') not in ReadyT:
                time.sleep(1)
        self.x = self.Position()

    def IsReadyT(self):
        return self.inst.query('1TS') in ReadyT

    def IsReady(self):
        return self.inst.query('1TS') in Ready

    def Position(self):
        while not self.IsReadyT():
            time.sleep(1)
        while re.match('1TP-*[0-9]+([.]+[0-9]*)*',self.inst.query('1TP')) == None:
            time.sleep(1)
        return float(re.findall('1TP-*([0-9]+([.]+[0-9]*)*)',self.inst.query('1TP'))[0][0])

    def MoveTo(self, position):
        self.WaitForReady()
        if position > self.safeRange[0] and position < self.safeRange[1]:
            self.inst.write('1PA'+str(position))
            while not self.WaitForReady():
                self.inst.write('1PA'+str(position))
            self.x = self.Position()
        else:
            print('Warning! Mirrors might crash, out of safe range')

    def Quasi_MoveTo(self, position):
        self.WaitForReady()
        if position > self.safeRange[0] and position < self.safeRange[1]:
            while np.absolute(self.Position() - position) > 1:
                self.MoveBy(-np.sign(self.Position() - position))
                self.WaitForReady()
            self.MoveBy(-(self.Position() - position))
        else:
            print('Warning! Mirrors might crash, out of safe range')

    def ErrorHandler(self, move_direction = None):
        state = self.inst.query('1TS')
        print(state)
        if state in Disable+DisableT:
            self.inst.write('1MM1')
        elif state in Tracking:
            if move_direction == None:
                self.inst.write('1ST')
            else:
                self.inst.write('1PR'+str(move_direction*0.01))
        state = self.inst.query('1TS')
        print('Error handled, current state:{state:s}'.format(state = state))

    def WaitForReady(self, move_direction = None):
        waitTime = 0;
        queryInterval = 1;#in unit of sec
        JobDone = True
        while not (self.IsReadyT() or self.IsReady()):
            time.sleep(queryInterval)
            waitTime += queryInterval
            if waitTime > 10:
                self.ErrorHandler(move_direction = move_direction)
                JobDone = False
                break
        return JobDone

    def MoveBy(self, displacement):
        target = self.x + displacement
        if target > self.safeRange[0] and target < self.safeRange[1]:
            self.WaitForReady()
            self.inst.write('1PR'+str(displacement))
            if displacement > 0:
                self.WaitForReady(move_direction = 1)
            else:
                self.WaitForReady(move_direction = -1)
            self.x = self.Position()
        else:
            print('Warning! Mirrors might crash, out of safe range')

    def close(self):
        self.inst.close()
