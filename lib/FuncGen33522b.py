import pyvisa as py

class FuncGen:

    def __init__(self):
        self.rm = py.ResourceManager()
        self.inst = self.rm.open_resource('USB0::0x0957::0x2C07::MY57802003::INSTR')

    def Sine(self, freq, amp, offset = 0, ch = 1):
        self.inst.write('SOUR'+str(ch)+':APPL:SIN '+freq+','+str(amp)+','+str(offset))

    def Noise(self, freq, amp, offset = 0, ch = 1):
        self.inst.write('SOUR'+str(ch)+'APPL:NOIS '+freq+','+str(amp)+','+str(offset))
