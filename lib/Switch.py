import time

class zx80:

    def __init__(self, RP):
        self.rp = RP
        self.RF1_RD_Pin = 'DIO0_P'
        self.RF2_RD_Pin = 'DIO2_P'
        self.rp.tx_txt('DIG:PIN:DIR out,'+self.RF1_RD_Pin)
        self.rp.tx_txt('DIG:PIN:DIR out,'+self.RF2_RD_Pin)

    def RF1On(self):
        self.rp.tx_txt('DIG:PIN '+self.RF1_RD_Pin+', 1')
        time.sleep(1e-3)
        self.rp.tx_txt('DIG:PIN? '+self.RF1_RD_Pin)
        assert(self.rp.rx_txt() == '1')

    def RF2On(self):
        self.rp.tx_txt('DIG:PIN '+self.RF2_RD_Pin+', 1')
        time.sleep(1e-3)
        self.rp.tx_txt('DIG:PIN? '+self.RF2_RD_Pin)
        assert(self.rp.rx_txt() == '1')

    def RF1Off(self):
        self.rp.tx_txt('DIG:PIN '+self.RF1_RD_Pin+', 0')
        time.sleep(1e-3)
        self.rp.tx_txt('DIG:PIN? '+self.RF1_RD_Pin)
        assert(self.rp.rx_txt() == '0')

    def RF2Off(self):
        self.rp.tx_txt('DIG:PIN '+self.RF2_RD_Pin+', 0')
        time.sleep(1e-3)
        self.rp.tx_txt('DIG:PIN? '+self.RF2_RD_Pin)
        assert(self.rp.rx_txt() == '0')

    def ComGround(self):
        self.RF1Off()
        self.RF2Off()
