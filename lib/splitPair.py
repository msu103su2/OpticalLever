import lib.CONEX_CC as cc
import numpy as np
from scipy.optimize import curve_fit
from scipy import special
import matplotlib.pyplot as plt
import datetime
import json

class splitPair:

    def __init__(self, safeRangeA, safeRangeB):
        self.A = cc.CONEX_CC('ASRL3::INSTR', safeRangeA, False)
        self.B = cc.CONEX_CC('ASRL6::INSTR', safeRangeB, False)
        self.Aangle = 72.5042/180*np.pi
        self.Bangle = 63.8/180*np.pi
        self.tiltAngle = 0.041
        self.centerXinA = 0
        self.centerXinB = 0
        # A B positive position is 'opposite'
        self.Quasi_gapsize = np.cos(self.tiltAngle)*abs(self.A.x*np.sin(self.Aangle) +\
            self.B.x*np.sin(self.Bangle))

    def MoveGapBy(self, displacement):
        #positive direction is the same as A
        if displacement > 0:
            self.B.MoveBy(-displacement/np.sin(self.Bangle))
            self.A.MoveBy(displacement/np.sin(self.Aangle))
        else:
            self.A.MoveBy(displacement/np.sin(self.Aangle))
            self.B.MoveBy(-displacement/np.sin(self.Bangle))
        self.Quasi_gapsize =np.cos(self.tiltAngle)*abs(self.A.x*np.sin(self.Aangle) +\
            self.B.x*np.sin(self.Bangle))

    def OpenGapBy(self, size):
        self.A.MoveBy(-size/(2*np.sin(self.Aangle)))
        self.B.MoveBy(-size/(2*np.sin(self.Bangle)))
        self.Quasi_gapsize = np.cos(self.tiltAngle)*abs(self.A.x*np.sin(self.Aangle) +\
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
        rows[:,0] = rows[:,0]/1000
        middleP = (rows[:,1].min()+rows[:,1].max())/2
        middleIndex = np.argmin(abs(rows[:,1]-middleP))
        middleX = rows[middleIndex,0]
        stepsToCenter = len(rows[:,1]) - middleIndex
        cutIndex = max(middleIndex - stepsToCenter,0)
        x = rows[cutIndex:,0]
        P = rows[cutIndex:,1]
        popt= curve_fit(func, x, P , p0 =[middleX, 75/180*np.pi, 2*middleP, self.tiltAngle], \
            bounds = ([0,0,0,self.tiltAngle*(1-1e-3)],[25,np.pi,3*middleP,self.tiltAngle*(1+1e-3)]))
        plt.plot(x, P, 'b-', label='data')
        plt.plot(x, func(x, *popt[0]),'g--', label='fit')
        plt.show()
        self.centerXinA = popt[0][0]
        self.Aangle = popt[0][1]
        print(self.Aangle*180/np.pi)
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

        allCaptured = False
        while not allCaptured:
            self.B.MoveBy(stepscale*self.B.min_step)
            tempPicoHandle.AutoRange()
            p = -np.mean(tempPicoHandle.getTimeSignal())
            rows = np.append(rows,[[self.B.x, p]], axis = 0)
            print(rows[-1])
            middleP = (rows[:,1].min()+rows[:,1].max())/2
            if (len(rows) > 10):
                if(abs(rows[-1,1] - rows[-10,1])/rows[-1,1] < 0.01 \
                    and abs(rows[-1,1]-2*middleP)/(2*middleP) < 0.05 ):
                    allCaptured = True
            else:
                allCaptured = False
        rows[:,0] = rows[:,0]/1000
        middleP = (rows[:,1].min()+rows[:,1].max())/2
        middleIndex = np.argmin(abs(rows[:,1]-middleP))
        middleX = rows[middleIndex,0]
        stepsToCenter = len(rows[:,1]) - middleIndex
        cutIndex = max(middleIndex - stepsToCenter,0)
        x = rows[cutIndex:,0]
        P = rows[cutIndex:,1]
        popt= curve_fit(func, x, P , p0 =[middleX, 66/180*np.pi, 2*middleP, self.tiltAngle], \
            bounds = ([0,0,0,self.tiltAngle*(1-1e-3)],[25,np.pi,3*middleP,self.tiltAngle*(1+1e-3)]))
        plt.plot(x, P, 'b-', label='data')
        plt.plot(x, func(x, *popt[0]),'g--', label='fit')
        plt.show()
        self.centerXinB = popt[0][0]
        self.Bangle = popt[0][1]
        print(self.Bangle*180/np.pi)
        self.B.MoveTo(x0)
        tempPicoHandle.close()

    def close(self):
        self.A.close()
        self.B.close()

    def getQuasi_gapsize(self):
        self.Quasi_gapsize = np.cos(self.tiltAngle)*abs(self.A.x*np.sin(self.Aangle) +\
            self.B.x*np.sin(self.Bangle))
        return self.Quasi_gapsize

    def Cali(self, mirror, ps5000a, stepscale, dir_filename, endpoint = None, direction = 1, sign = 1):
        if mirror == 'A':
            tg = self.A
            channel = 'A'
            angle = self.Aangle
        elif mirror =='B':
            tg = self.B
            channel = 'A'
            angle = self.Bangle
        laserP_minitor_channel = 'B'
        x0 = tg.x
        rows = np.ndarray((0,3), dtype = float)
        ps5000a.DC(channel)
        ps5000a.DC(laserP_minitor_channel)
        ps5000a.AutoRange(laserP_minitor_channel)
        angle_msnt_err = 0.05 # relative angle measuremnt error
        laserP_msnt_err = 0.1 #relative laser power measuerment error
        startTime = datetime.datetime.now()

        if not endpoint:
            allCaptured = False
            while not allCaptured:
                tg.MoveBy(stepscale*tg.min_step)
                ps5000a.AutoRange(channel)
                p = sign*np.mean(ps5000a.getTimeSignal(channel))
                cP = np.mean(ps5000a.chs[laserP_minitor_channel].timeSignal)
                rows = np.append(rows,[[tg.x, p, cP]], axis = 0)
                print(rows[-1])
                middleP = (rows[:,1].min()+rows[:,1].max())/2
                if (len(rows) > 10):
                    if(abs(rows[-1,1] - rows[-10,1])/rows[-1,1] < 0.01 \
                        and abs(rows[-1,1]-2*middleP)/(2*middleP) < 0.05 )\
                        and abs(rows[-1,1]) > 1:
                        allCaptured = True
                else:
                    allCaptured = False
        else:
            allCaptured = False
            while not allCaptured:
                tg.MoveBy(stepscale*tg.min_step*direction)
                ps5000a.AutoRange(channel)
                p = sign*np.mean(ps5000a.getTimeSignal(channel))
                cP = np.mean(ps5000a.chs[laserP_minitor_channel].timeSignal)
                rows = np.append(rows,[[tg.x, p, cP]], axis = 0)
                print(rows[-1])
                middleP = (rows[:,1].min()+rows[:,1].max())/2
                if (tg.x < endpoint and direction == 1 or tg.x > endpoint and direction == -1):
                    allCaptured = False
                else:
                    allCaptured = True
        endTime = datetime.datetime.now()
        rows[:,0] = rows[:,0]/1000
        middleP = (rows[:,1].min()+rows[:,1].max())/2
        middleIndex = np.argmin(abs(rows[:,1]-middleP))
        middleX = rows[middleIndex,0]
        stepsToCenter = len(rows[:,1]) - middleIndex
        cutIndex = max(middleIndex - 3*stepsToCenter,0)
        x = rows[cutIndex:,0]
        P = rows[cutIndex:,1]
        cP = rows[cutIndex:,2]
        popt, pcov = curve_fit(fitFunc, np.array(x), np.array(P)/np.array(cP)*np.mean(cP) , p0 =[middleX, angle, 2*middleP, self.tiltAngle, 1.25e-3, 0], \
            bounds = ([middleX*0.9, angle*(1-angle_msnt_err), 2*middleP*0.9, self.tiltAngle*(1-1e-3), 3e-4, -0.012],\
            [middleX*1.1, angle*(1+angle_msnt_err), 2*middleP*1.1, self.tiltAngle*(1+1e-3), 5e-3, 0.012]))
        print(popt)
        fitInfo = {
            'CaliMirror':mirror,
            'stagePositions':x.tolist(),
            'voltages':P.tolist(),
            'laserPowerMonitor':cP.tolist(),
            'fitResult':popt.tolist(),
            'fitCov':pcov.tolist(),
            'fitVariableNames':['middleX', 'angle', 'LaserPower','tiltAngle','beamWaist', 'offset'],
            'picoChannel':channel,
            'Mirror_tilt_angle':{'value':self.tiltAngle, 'unit':'radians'},
            'picoscopeInfo':ps5000a.getConfigureInfo(),
            'startTime':startTime.isoformat(),
            'endTime':endTime.isoformat()
            }
        with open(dir_filename,'w') as file:
            json.dump(fitInfo, file)
        plt.plot(x, P, 'b-', label='data')
        plt.plot(x, fitFunc(x, *popt),'g--', label='fit')
        plt.show()
        tg.MoveTo(x0)

    def RawCali(self, mirror, ps5000a, stepscale, dir_filename, endpoint = None, direction = 1):
        if mirror == 'A':
            tg = self.A
            channel = 'C'
            angle = self.Aangle
        elif mirror =='B':
            tg = self.B
            channel = 'C'
            angle = self.Bangle
        x0 = tg.x
        rows = np.ndarray((0,2), dtype = float)
        ps5000a.DC(channel)
        ps5000a.AutoRange(channel)
        angle_msnt_err = 0.05 # relative angle measuremnt error
        laserP_msnt_err = 0.1 #relative laser power measuerment error
        startTime = datetime.datetime.now()

        if not endpoint:
            allCaptured = False
            while not allCaptured:
                tg.MoveBy(stepscale*tg.min_step)
                ps5000a.AutoRange(channel)
                p = np.mean(ps5000a.getTimeSignal(channel))
                rows = np.append(rows,[[tg.x, p]], axis = 0)
                #print(rows[-1])
                middleP = (rows[:,1].min()+rows[:,1].max())/2
                if (len(rows) > 10):
                    if(abs(rows[-1,1] - rows[-10,1])/rows[-1,1] < 0.01 \
                        and abs(rows[-1,1]-2*middleP)/(2*middleP) < 0.05 )\
                        and abs(rows[-1,1]) > 1:
                        allCaptured = True
                else:
                    allCaptured = False
        else:
            allCaptured = False
            while not allCaptured:
                tg.MoveBy(stepscale*tg.min_step*direction)
                ps5000a.AutoRange(channel)
                p = np.mean(ps5000a.getTimeSignal(channel))
                rows = np.append(rows,[[tg.x, p]], axis = 0)
                #print(rows[-1])
                middleP = (rows[:,1].min()+rows[:,1].max())/2
                if (tg.x < endpoint and direction == 1 or tg.x > endpoint and direction == -1):
                    allCaptured = False
                else:
                    allCaptured = True
        endTime = datetime.datetime.now()
        rows[:,0] = rows[:,0]/1000
        middleP = (rows[:,1].min()+rows[:,1].max())/2
        middleIndex = np.argmin(abs(rows[:,1]-middleP))
        middleX = rows[middleIndex,0]
        stepsToCenter = len(rows[:,1]) - middleIndex
        cutIndex = max(middleIndex - 3*stepsToCenter,0)
        x = rows[cutIndex:,0]
        P = rows[cutIndex:,1]
        popt, pcov = curve_fit(fitFunc, np.array(x), np.array(P) , p0 =[middleX, np.pi/2, 2*middleP, 0, 1e-3, 0], \
            bounds = ([middleX*0.9, np.pi/2*0.999, 2*middleP*0.9, 0, 0, -0.012],\
            [middleX*1.1, np.pi/2*1.001, 2*middleP*1.1, 0.001, 0.1, 0.012]))
        print(popt)
        fitInfo = {
            'CaliMirror':mirror,
            'stagePositions':x.tolist(),
            'voltages':P.tolist(),
            'fitResult':popt.tolist(),
            'fitCov':pcov.tolist(),
            'fitVariableNames':['middleX', 'angle', 'LaserPower','tiltAngle','beamWaist', 'offset'],
            'picoChannel':channel,
            'Mirror_tilt_angle':{'value':self.tiltAngle, 'unit':'radians'},
            'picoscopeInfo':ps5000a.getConfigureInfo(),
            'startTime':startTime.isoformat(),
            'endTime':endTime.isoformat()
            }
        with open(dir_filename,'w') as file:
            json.dump(fitInfo, file)
        #plt.plot(x, P, 'b-', label='data')
        #plt.plot(x, fitFunc(x, *popt),'g--', label='fit')
        #plt.show()
        tg.MoveTo(x0)

    def AutoBalance(self, ps5000a, debug = False):
        ps5000a.DC('A')
        ps5000a.AutoRange('A')
        sign = -1

        ps5000a.getTimeSignal()
        data = ps5000a.getData()
        data = data[:,0]
        unbalance = sign*np.mean(data)
        std = np.std(data)

        while abs(unbalance)>3*std:
            self.MoveGapBy(sign*np.sign(unbalance)*0.005)
            ps5000a.getTimeSignal()
            data = ps5000a.getData()
            data = data[:,0]
            unbalance = np.mean(data)
            std = np.std(data)

        flag = True
        if abs(unbalance)>0.5:
            step = 0.005
        else:
            step = 0.002

        while abs(unbalance)>std:
            if unbalance>0:
                if flag:
                    self.B.MoveBy(-step)
                else:
                    self.A.MoveBy(step)
            else:
                if not flag:
                    self.B.MoveBy(step)
                else:
                    self.A.MoveBy(-step)
            ps5000a.AutoRange('A')
            ps5000a.AutoRange('B')
            ps5000a.getTimeSignal();
            data = ps5000a.getData()
            diff = sign*data[:,0]
            newUnbalance = np.mean(diff)
            std = np.std(diff)
            if debug:
                print(self.getQuasi_gapsize(),' ',newUnbalance)
                print(np.mean(data[:,0]),' ',np.mean(data[:,1]), ' ',np.mean(data[:,2]))
            if newUnbalance*unbalance<0:
                flag = not flag
            unbalance = newUnbalance

def func(x, x0, angle, power, tiltAngle):
    #numbers from 55mm focal lens with 30um beam waist on focal point
    return 0.5*power*(1+special.erf(1221.4*(x-x0)*np.sin(angle)*np.cos(tiltAngle)))

def fitFunc(x, x0, angle, power, tiltAngle, wx, offset):
    return offset + 0.5*power*(1+special.erf((x-x0)*np.sin(angle)*np.cos(tiltAngle)/(wx/np.sqrt(2))))
