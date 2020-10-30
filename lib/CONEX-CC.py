import pyvisa
import time
import re
ReadyT = ['1TS000036','1TS000037','1TS000038']
Ready = ['1TS000032','1TS000033','1TS000034']
class CONEX_CC:
    """docstring for CONEX_CC."""

    def __init__(self, deviceKey, reset = True):
        self.rm = pyvis.ResourceManager()
        self.inst = self.rm.open_resource(deviceKey, baud_rate = 821600,\
         data_bits = 8, write_termination = '\r\n', timeout = 2000)

        if self.inst.query('1TS') not in ReadyT:
            reset = True

        if reset:
            self.inst.write('RS')
            while self.inst.query('1TS')!='1TS00000A':
                time.sleep(1)
            self.inst.write('OR')
            while self.inst.query('1TS')!='1TS000032':
                time.sleep(1)
            self.inst.write('1TK1')
            while self.inst.query('1TS')!='1TS000036':
                time.sleep(1)

    def IsReadyT(self):
        self.inst.query('1TS') in ReadyT

    def IsReady(self):
        self.inst.query('1TS') in Ready

    def Position(self):
        while self.inst.query('1TS')!='1TS000037':
            time.sleep(1)
        while re.match('1TP[0-9]+[.]+[0-9]*',self.inst.query('1TP')) == None:
            time.sleep(1)
        return float(re.finall('1TP([0-9]+[.]+[0-9]*)',self.inst.query('1TP'))[0])

    def MoveTo(self, position):
        while not (self.IsReadyT() or self.IsReady()):
            time.sleep(1)
        self.inst.write('1PA'+str(position))

    def MoveBy(self, displacement):
        while not (self.IsReadyT() or self.IsReady()):
            time.sleep(1)
        self.inst.write('1PR'+str(displacement))
