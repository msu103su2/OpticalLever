import lib.redpitaya_scpi as scpi
import lib.Switch as Switch
import lib.AgilentNA as AG
import lib.ps5000A as ps
import lib.splitPair as sp
import lib.FuncGen33522b as FuncGen
import pyvisa as py
import numpy as np

wd = r'Z:\data\optical lever project\NORCADA_NX53515C\02-balancing'
rpip = '192.168.137.222'
rm = py.ResourceManager()

FG = FuncGen.FuncGen()
NA = AG.AgilentNA(wd)
RP = scpi.scpi(rpip)
DS = sp.splitPair()
zx80 = Switch.zx80(RP)
