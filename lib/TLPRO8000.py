import pyvisa as py
import re

class PRO8000:

    def __init__(self):
        self.rm = py.ResourceManager()
        self.inst = self.rm.open_resource('ASRL13::INSTR', baud_rate = 19200, data_bits = 8, \
            write_termination = '\r\n', read_termination = '\r\n')

    def setTs(self, Ts):
        self.inst.write(':TEMP:SET {Ts:f}'.format(Ts = Ts))

    def setP(self, P):
        self.inst.write(':SHAREP:SET {P:f}'.format(P = P))

    def getP(self):
        text = self.inst.query(':SHAREP:SET?')
        match = re.findall(':SHAREP:SET (.+)', text)
        return float(match[0])

    def setI(self, I):
        self.inst.write(':SHAREI:SET {I:f}'.format(I = I))

    def getI(self):
        text = self.inst.query(':SHAREI:SET?')
        match = re.findall(':SHAREI:SET (.+)', text)
        return float(match[0])

    def setD(self, D):
        self.inst.write(':SHARED:SET {D:f}'.format(D = D))

    def getD(self):
        text = self.inst.query(':SHARED:SET?')
        match = re.findall(':SHARED:SET (.+)', text)
        return float(match[0])

    def setCH(self, ch):
        self.inst.write(':SLOT {ch:d}'.format(ch = int(ch)))

    def getTa(self):
        text = self.inst.query(':TEMP:ACT?')
        match = re.findall(':TEMP:ACT (.+)', text)
        return float(match[0])

    def getTs(self):
        text = self.inst.query(':TEMP:SET?')
        match = re.findall(':TEMP:SET (.+)', text)
        return float(match[0])

    def getRa(self):
        text = self.inst.query(':RESI:ACT?')
        match = re.findall(':RESI:ACT (.+)', text)
        return float(match[0])

    def setThermistor(B, R0, T0):
        self.inst.write(':CALTB:SET {B:f}'.format(B = B))
        self.inst.write(':CALTR:SET {R0:f}'.format(R0 = R0))
        self.inst.write(':CALTT:SET {T0:f}'.format(T0 = T0))

    def getThermistor(B, R0, T0):
        self.inst.write(':CALTB:SET {B:f}'.format(B = B))
        self.inst.write(':CALTR:SET {R0:f}'.format(R0 = R0))
        self.inst.write(':CALTT:SET {T0:f}'.format(T0 = T0))

    def Ishare(self, on):
        if on:
            self.inst.write(':INTEG ON')
        else:
            self.inst.write(':INTEG OFF')

    def getITE(self):
        text = self.inst.query(':ITE:ACT?')
        match = re.findall(':ITE:ACT (.+)', text)
        return float(match[0])

    def getVTE(self):
        text = self.inst.query(':VTE:ACT?')
        match = re.findall(':VTE:ACT (.+)', text)
        return float(match[0])

    def setIlim(self, Ilim):
        self.write(':LIMT:SET {Ilim:f}'.format(Ilim = Ilim))

    def SwitchTEC(self, on):
        if on:
            self.inst.write(':TEC ON')
        else:
            self.inst.write(':TEC OFF')

    def close(self):
        self.inst.close()
