import lib.CONEX_CC as cc
import numpy as np
from scipy.optimize import curve_fit
from scipy import special
import lib.ps5000A as ps
import matplotlib.pyplot as plt

class splitPair:

    def __init__(self):
        self.A = cc.CONEX_CC('ASRL3::INSTR', False)
        self.B = cc.CONEX_CC('ASRL6::INSTR', False)
        self.Aangle = np.pi/2
        self.Bangle = np.pi/2
        self.centerXinA = 0
        self.centerXinB = 0
        # A B positive position is 'opposite'
        self.Quasi_gapsize = abs(self.A.x*np.sin(self.Aangle) +\
            self.B.x*np.sin(self.Bangle))

    def MoveGapBy(self, distance):
        #positive direction is the same as A
        self.A.MoveBy(distance/np.sin(self.Aangle))
        self.B.MoveBy(-distance/np.sin(self.Bangle))
        self.Quasi_gapsize = abs(self.A.x*np.sin(self.Aangle) +\
            self.B.x*np.sin(self.Bangle))

    def OpenGapBy(self, size):
        self.A.MoveBy(-size/(2*np.sin(self.Aangle)))
        self.B.MoveBy(-size/(2*np.sin(self.Bangle)))
        self.Quasi_gapsize = abs(self.A.x*np.sin(self.Aangle) +\
            self.B.x*np.sin(self.Bangle))

    def CaliAngleA(self, stepscale):
        input("Caution! Calibration need split D mirror pair to be far aprat!\
         \n Press Enter to continue...")
        x0 = self.A.x
        rows = np.ndarray((0,2), dtype = float)
        tempPicoHandle = ps.picoscope()
        tempPicoHandle.defaultSetting()
        tempPicoHandle.DC()
        tempPicoHandle.AutoRange()

        allCaptured = False
        while not allCaptured:
            self.A.MoveBy(stepscale*self.A.min_step)
            tempPicoHandle.AutoRange()
            p = np.mean(tempPicoHandle.getTimeSignal())
            rows = np.append(rows,[[self.A.x, p]], axis = 0)
            print(rows[-1])
            middleP = (rows[:,1].min()+rows[:,1].max())/2
            if (len(rows) > 10):
                if(abs(rows[-1,1] - rows[-10,1])/rows[-1,1] < 0.01 \
                    and abs(rows[-1,1]-2*middleP)/(2*middleP) < 0.05 ):
                    allCaptured = True
            else:
                allCaptured = False

        middleP = (rows[:,1].min()+rows[:,1].max())/2
        middleIndex = np.argmin(abs(rows[:,1]-middleP))
        middleX = rows[middleIndex,0]
        stepsToCenter = len(rows[:,1]) - middleIndex
        cutIndex = min(middleIndex - stepsToCenter,0)
        x = rows[cutIndex:,0]
        P = rows[cutIndex:,1]
        popt, pcov = curve_fit(func, x, P, p0 =[middleX, np.pi/2.7, 2*middleP], \
            bounds = (0,[25,np.pi,np.inf]))
        plt.plot(x, P, 'b-', label='data')
        plt.plot(x, func(x, *popt),'g--', label='fit')
        plt.show()
        self.centerXinA = popt[0][0]
        self.Aangle = popt[0][1]
        self.A.MoveTo(x0)
        tempPicoHandle.close()

    def CaliAngleB(self, stepscale):
        input("Caution! Calibration need split D mirror pair to be far aprat!\
         \n Press Enter to continue...")
        x0 = self.B.x
        rows = np.ndarray((0,2), dtype = float)
        tempPicoHandle = ps.picoscope()
        tempPicoHandle.defaultSetting()
        tempPicoHandle.DC()
        tempPicoHandle.AutoRange()

        while not allCaptured:
            self.B.MoveBy(stepscale*self.B.min_step)
            tempPicoHandle.AutoRange()
            p = -np.mean(tempPicoHandle.getTimeSignal())
            rows = np.append(rows,[[self.B.x, p]], axis = 0)
            middleP = (rows[:,1].min()+rows[:,1].max())/2
            if (len(rows) > 10):
                if(abs(rows[-1,1] - rows[-10,1])/rows[-1,1] < 0.01 \
                    and abs(rows[-1,1]-2*middleP)/(2*middleP) < 0.05 ):
                    allCaptured = True
            else:
                allCaptured = False

        middleP = (rows[:,1].min()+rows[:,1].max())/2
        middleIndex = np.argmin(abs(rows[:,1]-middleP))
        middleX = rows[middleIndex,0]
        stepsToCenter = len(rows[:,1]) - middleIndex
        cutIndex = min(middleIndex - stepsToCenter,0)
        x = rows[cutIndex:,0]
        P = rows[cutIndex:,1]
        popt, pcov = curve_fit(func, x, P , p0 =[middleX, np.pi/2.7, 2*middleP], \
            bounds = (0,[25,np.pi,np.inf]))
        plt.plot(x, P, 'b-', label='data')
        plt.plot(x, func(x, *popt),'g--', label='fit')
        plt.show()
        self.centerXinB = popt[0][0]
        self.Bangle = popt[0][1]
        self.B.MoveTo(x0)
        tempPicoHandle.close()

    def close(self):
        self.A.close()
        self.B.close()

def func(x, x0, angle, power):
    #numbers from 55mm focal lens with 30um beam waist on focal point
        return 0.5*power*(1+special.erf(2274.97*(x-x0)*np.sin(angle)))
