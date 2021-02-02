import lib.ps5000A as ps
import lib.splitPair as sp
ps5000a = ps.picoscope()
ps5000a.defaultSetting()
ps5000a.DC('A')
ps5000a.DC('B')
ps5000a.AutoRange('A')
ps5000a.AutoRange('B')
ps5000a.configurePSD(8.515, 4e6)
DS = sp.splitPair()
DS.Cali('A', ps5000a, 300)
