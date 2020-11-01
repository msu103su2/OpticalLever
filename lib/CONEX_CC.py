import pyvisa
import time
import re
ReadyT = ['1TS000036','1TS000037','1TS000038']
Ready = ['1TS000032','1TS000033','1TS000034']
class CONEX_CC:
    """docstring for CONEX_CC."""

    def __init__(self, deviceKey, reset = True):
        self.rm = pyvisa.ResourceManager()
        self.inst = self.rm.open_resource(deviceKey, baud_rate = 921600,\
            data_bits = 8, write_termination = '\r\n', read_termination = '\r\n',\
            timeout = 2000)
        self.min_step = 0.0001

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
        while not (self.IsReadyT() or self.IsReady()):
            time.sleep(1)
        self.inst.write('1PA'+str(position))
        while not (self.IsReadyT() or self.IsReady()):
            time.sleep(1)
        self.x = self.Position()

    def MoveBy(self, displacement):
        while not (self.IsReadyT() or self.IsReady()):
            time.sleep(1)
        self.inst.write('1PR'+str(displacement))
        while not (self.IsReadyT() or self.IsReady()):
            time.sleep(1)
        self.x = self.Position()

    def close(self):
        self.inst.close()
