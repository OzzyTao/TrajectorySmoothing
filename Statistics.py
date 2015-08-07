'''
Created on 14/06/2013

@author: ljguan
'''

import random
import numpy
import math
import fileinput
import os
import re

class Statistics():
    '''
    This function samples trajectories along each path in segls
    and determines the distance along each path at the measurement sampling
    times.
    '''


    def __init__(self):
        self.sig = 5    
    
    
    
    def densityest(self, x,y):
        '''A normal kernel with optimal bandwidth is used.
        estimate density of x at points y'''
        
        m = len(y)
        n = len(x)
        h = 1.06*numpy.std(x)/(n**(1.0/5.0)) 
        
        # Initialize fhat as list of zeros
        fhat = numpy.zeros((m));
          
        for k in range(m):
            sumvalue = 0
            for i in range(n):
                epxvalue = x[i]-y[k]
                if (abs(epxvalue)/h)<15:
                    epxvalue = x[i]-y[k];
                    nu = epxvalue * epxvalue / (h * h)
                    npdf = math.exp(-0.5*nu) / (math.sqrt(2*math.pi)*h)
                    sumvalue = sumvalue + npdf;
            fhat[k] = sumvalue/n;
              
        return fhat;
    
    def getMeanSegLength(self, segls):
        seglength = 0
        for seg in segls:
            seglength += seg.sum(axis=0)
        return seglength/segls.shape[0]
    
    def getMaxSegNum(self, segls):
        num = 0
        for seg in segls:
            if seg.shape[0] > num:
                num = seg.shape[0]
        return num
    
    
    def getReshapeMatrix(self, segls):
        '''fill a new array with a predefined row and column by an old array'''
        segls1 = numpy.zeros((segls.shape[0],self.getMaxSegNum(segls)))
        for i in range(len(segls)):
            for j in range(len(segls[i])):
                #print str(i) + "," + str(j)
                segls1[i,j] = segls[i][j]
        return segls1
    
    
    def arclgen3(self, obstimes,segls):
        '''This function samples trajectories along each path in segls and determines
           the distance along each path at the measurement sampling times.'''
        # [datamat,isvalid,pathprior] = arclgen3(obstimes,segls)
        # Samples paths and determines the distance along each path at the measurement sampling times
        # Inputs
        # obstimes: times at which way point observations are available (includes start and end points)
        # segls: lengths for each segment of each path
        # Outputs
        # datamat: distance along each sampled path at the measurement sampling times
        # isvalid: indicator for validity of the sampled path
        # pathprior: prior probability for each path
                
        #reshape the matrix
        segls = self.getReshapeMatrix(segls)
    
        nobs = len(obstimes)
        ssize = 2
        szsegls = segls.shape
        ndsegls = segls.ndim
        npaths = szsegls[0]
        if ndsegls==2:
            npaths = szsegls[0]
            maxsegs = szsegls[1]
            sumsegls = segls.sum(axis=1)
            cumsumsegls = segls.cumsum(axis=1)
            mnsegls = sumsegls.sum(axis=0)/npaths
        else:
            npaths = 1
            maxsegs = szsegls[0]
            sumsegls = segls.sum(axis=0)
            cumsumsegls = segls.cumsum(axis=0)
            mnsegls = sumsegls
        
#        print cumsumsegls
        # Velocity profile
        t = obstimes[nobs-1]-obstimes[0]
        vmax = 60/3.6
        r = (mnsegls/t)/vmax
        b = 1.5
        a = b*r/(1-r)
        
        # Ratio for sampling start and end points
        rho = 0.05
        
        # Fudge factor for path acceptance
        eps = max(1,t/100);
        
        maxiter = 10;
        
        # This is needed for input to densityest (seems a bit clumsy)
        tin = numpy.zeros((1))
        tin[0] = t
        
        # Variable initialisation
        vs = numpy.zeros((maxsegs))
        ts = numpy.zeros((maxsegs))
        tdata = numpy.zeros((ssize))
        tkerni = numpy.zeros((npaths))
        datamat = numpy.zeros((ssize*npaths,nobs))
        isvalid = numpy.zeros((ssize*npaths))
        for i in range(npaths):
            for j in range(ssize):
                count = 0
                isok = 0
                d0 = random.random()*rho*sumsegls[i]
                dmp1 = (1-rho+rho*random.random())*sumsegls[i]
                while ((count<maxiter) and (not isok)):
                    count = count+1
                    # Determine travel time
                    # Start node of starting segment
                    nd0 = 0
                    while (d0>cumsumsegls[i,nd0]):
                        nd0 = nd0+1
                    # Start node of end segment
                    ndmp1 = maxsegs-1
                    while (dmp1<cumsumsegls[i,ndmp1]):
                        ndmp1 = ndmp1-1
                    ndmp1 = ndmp1+1
                    for k in range(nd0):
                        ts[k] = 0
                    # Starting segment
                    vs[nd0] = vmax*random.betavariate(a,b)
                    ts[nd0] = (cumsumsegls[i,nd0]-d0)/vs[nd0]
                    # Middle segments
                    for k in numpy.arange(nd0+1,ndmp1):
                        vs[k] = vmax*random.betavariate(a,b)
                        ts[k] = segls[i,k]/vs[k]
                    # Final segment
                    vs[ndmp1] = vmax*random.betavariate(a,b)
                    ts[ndmp1] = (dmp1-cumsumsegls[i,ndmp1-1])/vs[ndmp1]
                    for k in numpy.arange(ndmp1+1,maxsegs):
                        ts[k] = 0
                
                    tcum = ts.cumsum(axis=0)   
                    tdata[j] = tcum[maxsegs-1]
                    if ((tdata[j]>(t-eps)) and (tdata[j]<(t+eps))):
                        isvalid[i*ssize+j] = 1
                        isok = 1;
                        datamat[i*ssize+j,0] = d0
                        for k in range(nobs-2):        
                            seg = 0
                            d = d0
                            while obstimes[k]>tcum[seg]:
                                d = d+segls[i,seg]
                                seg = seg + 1;
                            if seg>0:
                                d = d+(obstimes[k]-tcum[seg-1])*vs[seg]
                            else:
                                d = d0+obstimes[k]*vs[0]
                            datamat[i*ssize+j,k+1] = d
                        datamat[i*ssize+j,nobs-1] = dmp1 
            tkerni[i] = self.densityest(tdata,tin)
        priornorm = tkerni.sum(axis=0)
        pathprior = numpy.zeros((npaths))
        for i in range(npaths):
            pathprior[i] = tkerni[i]/priornorm
        
        return datamat, isvalid, pathprior
    
    
    # [pP1,pP2] = calcpostprobs2(z,x,isvalid,pP0,sig)
    # Computes the posterior probability of each path
    # Start and end points are included
    # Inputs
    # z: measurements (includes the start and end points)
    # x: positions along each sample path (should include start and end points)
    # isvalid: indicator for the validity of each sample path
    # pP0: path prior
    # sig: measurement noise standard deviation
    # Outputs
    # pP1: path posterior using only start and end point
    # pP2: path posterior using all measurements

    def calcpostprobs2(self, z,x,isvalid,pP0,sig):
        szz = z.shape
        npts = szz[1]
        npaths = len(pP0)
        szx = x.shape
        ssize = szx[0]/npaths  
        
#        sumx = abs(x).sum(axis=1)    
        logfs1 = numpy.zeros((npaths,ssize))
        logfs2 = numpy.zeros((npaths,ssize))
#        PPinum = numpy.zeros((npaths))
        nsuccess = numpy.zeros((npaths))
        maxlog1 = -1e10
        maxlog2 = -1e10
        for i in range(npaths):
            for j in range(ssize):
                idx = i*ssize+j
                if isvalid[idx]:
                    epx = z[0,0]-x[idx,0]
                    nux = epx*epx/(sig*sig)
                    epy = z[1,0]-x[idx,0]
                    nuy = epy*epy/(sig*sig)
                    logfs1[i,nsuccess[i]] = -0.5*(nux+nuy)
                    logfs2[i,nsuccess[i]] = logfs1[i,nsuccess[i]]
                    for k in numpy.arange(1,npts-1):
                        epx = z[0,k]-x[idx,2*k]
                        nux = epx*epx/(sig*sig)
                        epy = z[1,k]-x[idx,2*k+1]
                        nuy = epy*epy/(sig*sig)
                        logfs2[i,nsuccess[i]] = logfs2[i,nsuccess[i]]-0.5*(nux+nuy)
                    k = npts-1
                    epx = z[0,k]-x[idx,2*k]
                    nux = epx*epx/(sig*sig)
                    epy = z[1,k]-x[idx,2*k+1]
                    nuy = epy*epy/(sig*sig)
                    logfs2[i,nsuccess[i]] = logfs2[i,nsuccess[i]]-0.5*(nux+nuy)
                    logfs1[i,nsuccess[i]] = logfs1[i,nsuccess[i]]-0.5*(nux+nuy)
                    if logfs1[i,nsuccess[i]]>maxlog1:
                        maxlog1 = logfs1[i,nsuccess[i]]
                    if logfs2[i,nsuccess[i]]>maxlog2:
                        maxlog2 = logfs2[i,nsuccess[i]]
                    nsuccess[i] = nsuccess[i]+1
              
        pPnum1 = numpy.zeros(npaths)
        pPnum2 = numpy.zeros(npaths)
        for i in range(npaths):
            for j in range(nsuccess[i]):
                pPnum1[i] = pPnum1[i]+math.exp(logfs1[i,j]-maxlog1)
                pPnum2[i] = pPnum2[i]+math.exp(logfs2[i,j]-maxlog2)
            pPnum1[i] = pPnum1[i]*pP0[i]/ssize
            pPnum2[i] = pPnum2[i]*pP0[i]/ssize
    
        pPnorm1 = pPnum1.sum(axis=0)
        pPnorm2 = pPnum2.sum(axis=0)
        pP1 = []
        pP2 = []
        for i in range(npaths):
            pP1[i] = pPnum1[i]/pPnorm1
            pP2[i] = pPnum2[i]/pPnorm2
          
        return pP1, pP2
    
    def arclgen2(self, obstimes,t,segls):
        # [datamat,pathprior] = arclgen2(obstimes,t,segls)
        # Inputs:
        # obstimes: times at which way point observations are available (numpy.array([20,45]))
        # t: duration of surveillance, (t=80)
        # segls: lengths for each segment of each path (segls = numpy.array([[125,225,500,150],[250,360,200,0]]))
        # Outputs:
        # datamat: distance along each path for each measurement sampling time
        # pathprior: prior probability for each path
        nobs = len(obstimes)
        ssize = 100
#        ndsegls = segls.ndim
        ndsegls = 2 #ndsegls should always bigger than 1 if more than one path exist. 

        #TODO:check if only one path was generated.
        npaths = segls.shape[0]
        maxsegs = self.getMaxSegNum(segls)
#        sumsegls = segls.sum(axis=1)
#        mnsegls = sumsegls.sum(axis=0)/npaths
        mnsegls = self.getMeanSegLength(segls)
        
        # Velocity profile
        vmax = 60/3.6
        r = (mnsegls/t)/vmax  #ratio of average driving speed to max one
        b = 1.5
        a = b*r/(1-r)
        
        
        # Fudge factor for path acceptance
        eps = max(1,t/100);
        
        maxiter = 10;
        
        # This is needed for input to densityest (seems a bit clumsy)
        tin = numpy.zeros((1))
        tin[0] = t
        
        
        #Variable initialisation
        vs = numpy.zeros((maxsegs))
        ts = numpy.zeros((maxsegs))
        tdata = numpy.zeros((ssize))
        tkerni = numpy.zeros((npaths))
        datamat = numpy.zeros((ssize*npaths,nobs))
        
        #reshape the matrix
        segls = self.getReshapeMatrix(segls)
        
        
        for i in range(npaths):
#            i
            for j in range(ssize):
                count = 0
                isok = 0
                while ((count<maxiter) and (not isok)):
                    count = count+1
                    for k in range(maxsegs):
                        vs[k] = vmax*random.betavariate(a,b)
                        if ndsegls>1:
                            ts[k] = segls[i,k]/vs[k]
                        else:
                            ts[k] = segls[k]/vs[k]
                    tcum = ts.cumsum(axis=0)   
                    tdata[j] = tcum[maxsegs-1]
                    if ((tdata[j]>(t-eps)) and (tdata[j]<(t+eps))):
                        isok = 1;
                        for k in range(nobs):
                            seg = 0
                            d = 0
                            while obstimes[k]>tcum[seg]:
                                if ndsegls>1:
                                    d = d+segls[i,seg]
                                else:
                                    d = d+segls[seg]
                                seg = seg + 1;
                            if seg>0:
                                d = d + (obstimes[k]-tcum[seg-1])*vs[seg]
                            else:
                                d = obstimes[k]*vs[0]
                            datamat[i*ssize+j,k] = d 
            tkerni[i] = self.densityest(tdata,tin)
        priornorm = tkerni.sum(axis=0)
        pathprior = numpy.zeros((npaths))
        for i in range(npaths):
            pathprior[i] = tkerni[i]/priornorm

        return datamat, pathprior
        
      
    # [PPi] = calcpostprobs(z,x,pP0)
    # Computes the posterior probability of each path
    # Inputs
    # z: measurements (includes the start and end points)
    # x: positions along each sample path (should include start and end points)
    # pP0: path prior
    # PPi: posterior probability for each path   
        
    def calcpostprobs(self,z,x,pP0):
        szz = z.shape
        npts = szz[1]
        npaths = len(pP0)
        szx = x.shape
        ssize = szx[0]/npaths  

        sumx = abs(x).sum(axis=1)    
        logfs = -2e10*numpy.ones((npaths,ssize))
        PPinum = numpy.zeros((npaths))
        nsuccess = numpy.zeros((npaths))
        for i in range(npaths):
            for j in range(ssize):
                idx = i*ssize+j
                if sumx[idx]>1e-12:
                    nsuccess[i] = nsuccess[i]+1
                    logfs[i,nsuccess[i]] = 0
                    for k in range(npts):
                        epx = z[0,k]-x[idx,2*k]
                        nux = epx*epx/(self.sig*self.sig)
                        epy = z[1,k]-x[idx,2*k+1]
                        nuy = epy*epy/(self.sig*self.sig)
                        logfs[i,nsuccess[i]] = logfs[i,nsuccess[i]]-0.5*(2*math.log(2*math.pi*self.sig*self.sig)+nux+nuy)
        
        maxlog = numpy.max(logfs)
        pPnum = numpy.zeros(npaths)
        for i in range(npaths):
            for j in range(nsuccess[i]):
                pPnum[i] = pPnum[i]+math.exp(logfs[i,j]-maxlog)
            pPnum[i] = pPnum[i]*pP0[i]/ssize
        
        PPinorm = pPnum.sum(axis=0)
        pP = []
        for i in range(npaths):
            pP[i] = PPinum[i]/PPinorm;
        return pP
    


    
        #function [datamat,pathprior] = arclgen2(meas,est,segls)
        #
        #% [datamat,pathprior] = arclgen2(meas,est,segls)
        #% This function samples trajectories along each path in segls
        #% and determines the distance along each path at the measurement sampling
        #% times.
        #% Inputs:
        #% meas: measurement data
        #% est: estimator data
        #% segls: lengths for each segment of each path
        #% Outputs:
        #% datamat: distance along each path for each measurement sampling time
        #% pathprior: prior probability for each path
        #    def arclgen2(self, meas, est, segls):

    
    def readoutput1(self):
        fp = os.path.join("/Users/mrmore/Dropbox/DSTO/case_one/", "output1.txt")
        dList = []
        d = []
        for line in fileinput.input(fp, inplace=0):
            d = re.split('[(\)[\]\n]', line)
            d = filter(None, d)
            d.pop(0)
            dList.append(self.extractArcLength(d))
        
        self.dPath = numpy.array(dList)
#        self.dPath = dList    
#            dtime = re.split('[- :]', str(feat.attribute('timestamp').toPyObject()))
    
    def extractArcLength(self, strList):
        dist = []
        for dstr in strList:
            dist.append(float(dstr.split(",").pop()))
        return numpy.array(dist, dtype=float)
    
    def readRawPoints(self):
#        d = csvread(fname); #read matrix out from file
    
        fp = os.path.join("/Users/mrmore/Dropbox/DSTO/case_one/", "output0.txt")
        dList = []
        for line in fileinput.input(fp, inplace=0):
            dList.append(line.split(","))
        
        self.dArray = numpy.array(dList, dtype=float)
        self.duration = self.dArray[-1,2] - self.dArray[0,2]
        
#        npts = size(d,1);   #the number of points
#        xs = [0 0]';         #init xs
#        x0 = d(1,1:2).';     #start point(x,y)
#        ts = 0;              #
#        xe = (d(npts,1:2)-d(1,1:2)).'; #
#        te = d(npts,3)-d(1,3);         #
#        waypts = (d(2:npts-1,1:2)-ones(npts-2,1)*d(1,1:2)).'; get relative coordinates for all waypoints
#        tpts = d(2:npts-1,3).'-d(1,3);
        
#        self.npts = fileinput.filelineno()
        self.npts = self.dArray.shape[0]
        self.xs = [0,0]
        
    
def test_readRawPoints(sts):
    sts.readRawPoints()
    print sts.npts
#    print sts.dArray
    print sts.duration
    print sts.dArray[-1] - sts.dArray[0]
    
sts = Statistics()
#test_readRawPoints(sts)

sts.readRawPoints()
sts.readoutput1()
#print sts.dPath
#print sts.dPath.ndim
obstimes = numpy.array([sts.dArray[0][2],sts.dArray[-1][2]])
#print sts.getMeanSegLength(sts.dPath)
#print sts.getMaxSegNum(sts.dPath)
(datamat, isvalid, pathprior) = sts.arclgen3(obstimes, sts.dPath)
#print sts.getReshapeMatrix(sts.dPath)
print datamat.shape
#print datamat
print isvalid.sum()
print pathprior







 
 
 
 
 
 
 
 
 
 
    
    
    
        
        
        
        
        