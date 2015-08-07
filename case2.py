import random
import roadnet
import numpy
# import densityest
import math
import mystats
from scipy.stats import norm
import case1
from math import pi
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d.axes3d import Axes3D
import time
import testdb


def calcinter(G,z1,t1,z2,t2,n,n2,K,tax,ngrid,vmean,vvar,ups):
    
# [H,P1,P2,s1,eidx1] = calcinter(G,z1,t1,z2,t2,n,n2,K,tax,ngrid,vmean,vvar,ups)
# Calculates the Helliner affinity between the posterior position densities of two 
# objects at specified times
# Inputs:
# G- road network
# z1- m1-length list of position measurements for object 1 (list of list of floats)
# t1- list of measurement times for object 1 (list of floats)
# z2- m2-length list of position measurements for object 2 (list of list of floats)
# t2- m2-length list of measurement times for object 2 (list of floats)
# n- sample size for posterior density approximation (integer)
# n2- sample size for Hellinger affinity approximation (integer)
# K- number of candidate paths (integer)
# tax- times at which to compute Hellinger affinity (list of floats)
# ngrid- number of intervals used for numerical approximation to Hellinger affinity
# (vmean,vvar)- statistics of vehicle movement (floats)
# ups- measurement noise variance (float)
# Outputs:
# H- Hellinger affinity at times (t1,t2) (ntax x ntax array of floats)
# P1- position probabilities for object 1 in each graph interval (ngrid x ntax 
# array of floats)
# P2- position probabilities for object 2 in each graph interval (ngrid x ntax 
# array of floats)
# s1- sampled positions for object 1 at each time (ntax-length list of n2-length 
# arrays of floats)
# eidx1- sampled edges for object 1 at each time (ntax-length list of n2-length 
# lists of integers)
    
    # Sample paths and positions for both objects
    [pp1,ppos1,wp1,segp1,sp1] = case1.calcpostprobs_case2(G,z1,t1,n,K,vmean,vvar,ups)
    [pp2,ppos2,wp2,segp2,sp2] = case1.calcpostprobs_case2(G,z2,t2,n,K,vmean,vvar,ups)
    #w1cdf = numpy.cumsum(wp1)
    #w2cdf = numpy.cumsum(wp2)    
    
    # Pre-compute cumulative length along the network edges
    # This is used for probability calculations for graph intervals
    e = G.edges()
    nedges = len(e)
    ell = numpy.zeros(nedges)
    ellcum = numpy.zeros(nedges)
    i = 0
    pos0 = G.node[e[i][0]]["pos"]
    pos1 = G.node[e[i][1]]["pos"]
    ell[0] = numpy.linalg.norm(pos1-pos0)
    ellcum[0] = 0.
    for i in range(1,nedges):
        pos0 = G.node[e[i][0]]["pos"]
        pos1 = G.node[e[i][1]]["pos"]
        ell[i] = numpy.linalg.norm(pos1-pos0)
        ellcum[i] = ellcum[i-1]+ell[i-1]
    totell = ellcum[nedges-1]+ell[nedges-1]
    
    # Grid size
    gsize = totell/ngrid

    # Pre-allocate required quantities
    nt = len(tax)
    s1 = [numpy.zeros(n2) for j in range(nt)]
    eidx1 = [[0 for i in range(n2)] for j in range(nt)]
    w1 = [[] for j in range(nt)]
    lw1 = numpy.zeros(n2)
    gpos1 = numpy.zeros(n2)
    w1tilde = numpy.zeros(n2)
    e1 = [[0 for i in range(2)] for j in range(n2)]
    s2 = numpy.zeros(n2)
    lw2 = numpy.zeros(n2)
    gpos2 = numpy.zeros(n2)
    w2tilde = numpy.zeros(n2)
    e2 = [[0 for i in range(2)] for j in range(n2)]
    H = numpy.zeros((nt,nt))
    P1 = numpy.zeros((ngrid,nt))
    P2 = numpy.zeros((ngrid,nt))
    # Loop over times at which to compute the Hellinegr affinity
    # Here the position probabilities are calculated for each time
    for j in range(nt):
        t = tax[j]
        # Find the bracketing measurement times for this time
        k1 = 0
        while t1[k1]<t:
            k1 += 1
        t1k = t1[k1]
        t1km1 = t1[k1-1]
        k2 = 0
        while t2[k2]<t:
            k2 += 1
        t2k = t2[k2]
        t2km1 = t2[k2-1]
        maxw1 = -1e100
        maxw2 = -1e100
        # Draw sample positions for each object at this time
        for i in range(n2):
            # Sample a position for object 1
            idx1 = mystats.drawmultinom(wp1)
            seg1k = segp1[k1][idx1]
            s1k = sp1[k1,idx1]
            seg1km1 = segp1[k1-1][idx1]
            s1km1 = sp1[k1-1,idx1]
            [s1[j][i],seg1,lw1[i]] = draw_location(t,ppos1[idx1],s1k,seg1k,s1km1,seg1km1,t1k,t1km1,vmean,vvar)
            e1[i][0] = pp1[idx1][seg1]
            e1[i][1] = pp1[idx1][seg1+1]
            # Find the edge index and convert to a position along the line 
            # representation of the network
            if tuple(e1[i]) in e:
                eidx1[j][i] = e.index(tuple(e1[i]))
                gpos1[i] = ellcum[eidx1[j][i]]+s1[j][i]
            else:
                eidx1[j][i] = e.index(tuple([e1[i][1],e1[i][0]]))
                gpos1[i] = ellcum[eidx1[j][i]]+ell[eidx1[j][i]]-s1[j][i]
            if lw1[i]>maxw1:
                maxw1 = lw1[i]
            # Do the same for object 2            
            idx2 = mystats.drawmultinom(wp2)
            seg2k = segp2[k2][idx2]
            s2k = sp2[k2,idx2]
            seg2km1 = segp2[k2-1][idx2]
            s2km1 = sp2[k2-1,idx2]
            [s2[i],seg2,lw2[i]] = draw_location(t,ppos2[idx2],s2k,seg2k,s2km1,seg2km1,t2k,t2km1,vmean,vvar)
            e2[i][0] = pp2[idx2][seg2]
            e2[i][1] = pp2[idx2][seg2+1]
            if tuple(e2[i]) in e:
                eidx2 = e.index(tuple(e2[i]))
                gpos2[i] = ellcum[eidx2]+s2[i]
            else:
                eidx2 = e.index(tuple([e2[i][1],e2[i][0]]))
                gpos2[i] = ellcum[eidx2]+ell[eidx2]-s2[i]
            if lw2[i]>maxw2:
                maxw2 = lw2[i]
        # Normalise the weights
        for i in range(n2):
            w1tilde[i] = math.exp(lw1[i]-maxw1)
            w2tilde[i] = math.exp(lw2[i]-maxw2)    
        w1[j] = w1tilde/sum(w1tilde)
        w2 = w2tilde/sum(w2tilde)
        # Compute position probabilities for each interval
        for i in range(n2):
            idx = int(math.floor(gpos1[i]/gsize))
            P1[idx,j] += w1[j][i]
            idx = int(math.floor(gpos2[i]/gsize))
            P2[idx,j] += w2[i]
    
    # Compute Hellinger affinity for each pair of times
    for j1 in range(nt):
        for j2 in range(nt):
            for k in range(ngrid):
                H[j1,j2] += math.sqrt(P1[k,j1]*P2[k,j2])
                
    
    return H, P1, P2, s1, eidx1
                
                

def draw_location(t,path_pos,sk,segk,skm1,segkm1,tk,tkm1,vmean,vvar):
    
# [ssamp,segidx,logw] =  
# draw_location(t,path_pos,sk,segk,skm1,segkm1,tk,tkm1,vmean,vvar):
# Sample a location between two points on the network
# Inputs:
# t- time at which sample is being drawn (float)
# path_pos- array of positions of path nodes (2 x q array of floats)
# (sk,segk)- position at time tk>t (float, integer)
# (skm1,segkm1)- position at time tkm1<t (float, integer)
# tk- next time (float)
# tkm1- previous time (float)
# (vmean,vvar)- - statistics of vehicle movement (floats)
# Outputs:
# ssamp- sample position on edge (float)
# segidx- index of sample along path (integer)
# logw- sample weight (float)
    

    # Pre-allocate
    nsegs = segk-segkm1+1 # This is the number of segments between the bracketing positions
    Tbar = numpy.zeros((nsegs))
    kap = numpy.zeros((nsegs))
    ell = numpy.zeros((nsegs)) # total length of each segment
    lb = numpy.zeros((nsegs))
    ub = numpy.zeros((nsegs))
    dist_cs = numpy.zeros((nsegs))
    
    totell = 0
    # For each segment, compute its length, mean duration and duration variance
    if nsegs==1:
        pos1 = path_pos[:,segkm1]
        pos2 = path_pos[:,segkm1+1]
        ell[0] = numpy.linalg.norm(pos1-pos2)
        Tbar[0] = ell[0]/vmean
        kap[0] = ell[0]*ell[0]*vvar/(vmean**4)
        lb[0] = skm1
        ub[0] = sk
        totell = sk-skm1
        dist_cs[0] = totell
    else:
        pos1 = path_pos[:,segkm1]
        pos2 = path_pos[:,segkm1+1]
        ell[0] = numpy.linalg.norm(pos1-pos2)
        Tbar[0] = ell[0]/vmean
        kap[0] = ell[0]*ell[0]*vvar/(vmean**4)
        lb[0] = skm1
        ub[0] = ell[0]
        totell = ell[0]-skm1
        dist_cs[0] = totell
        for j in range(1,nsegs-1):
            pos1 = path_pos[:,segkm1+j]
            pos2 = path_pos[:,segkm1+j+1]
            ell[j] = numpy.linalg.norm(pos1-pos2)
            Tbar[j] = ell[j]/vmean
            kap[j] = ell[j]*ell[j]*vvar/(vmean**4)
            lb[j] = 0
            ub[j] = ell[j]
            totell = totell+ell[j]
            dist_cs[j] = dist_cs[j-1]+ell[j]
        j = nsegs-1
        pos1 = path_pos[:,segk-1]
        pos2 = path_pos[:,segk]
        ell[j] = numpy.linalg.norm(pos1-pos2)
        Tbar[j] = ell[j]/vmean
        kap[j] = ell[j]*ell[j]*vvar/(vmean**4)
        lb[j] = 0
        ub[j] = sk
        totell = totell+sk
        dist_cs[j] = dist_cs[j-1]+sk
    
    
    shat = numpy.zeros((nsegs))
    lam = numpy.zeros((nsegs))
    etilde = numpy.zeros((nsegs))
    p2 = numpy.zeros((nsegs))
    xi1 = numpy.zeros((nsegs))
    xi2 = numpy.zeros((nsegs))
    a1 = numpy.zeros((nsegs))
    a2 = numpy.zeros((nsegs))
    Ttilde1 = numpy.zeros((nsegs))
    Ttilde2 = numpy.zeros((nsegs))
    scfact = numpy.zeros(nsegs)
    scfact[0] = (ell[0]-skm1)*(ell[0]-skm1)/(ell[0]*ell[0])
    for b in range(1,nsegs-1):
        scfact[b] = 1.
    scfact[nsegs-1] = sk*sk/(ell[nsegs-1]*ell[nsegs-1])
    
    # Compute the weights and sampling density for each segment
    # First segment
    a1[0] = -Tbar[0]/ell[0]
    Ttilde1[0] = sk*Tbar[nsegs-1]/ell[nsegs-1]
    for b in range(0,nsegs-1):
        Ttilde1[0] += Tbar[b]
    xi1[0] = kap[0]
    for b in range(1,nsegs):
        xi1[0] += scfact[b]*kap[b]
    a2[0] = Tbar[0]/ell[0]
    Ttilde2[0] = -skm1*Tbar[0]/ell[0]
    xi2[0] = kap[0]
    db = (a1[0]*a1[0]*xi2[0]+a2[0]*a2[0]*xi1[0])
    shat[0] = (a2[0]*xi1[0]*(t-tkm1-Ttilde2[0])+a1[0]*xi2[0]*(tk-t-Ttilde1[0]))/db
    lam[0] = xi1[0]*xi2[0]/db
    p1 = mystats.normpdf(a2[0]*(tk-t-Ttilde1[0]),a1[0]*(t-tkm1-Ttilde2[0]),math.sqrt(db))
    p2[0] = mystats.normprob(shat[0],math.sqrt(lam[0]),lb[0],ub[0])
    etilde[0] = p1*p2[0]
    for b in range(1,nsegs):
        a1[b] = -Tbar[b]/ell[b]
        a2[b] = Tbar[b]/ell[b]
        Ttilde2[b] = Ttilde2[b-1]+Tbar[b-1]
        Ttilde1[b] = Ttilde1[b-1]-Tbar[b-1]
        xi2[b] = xi2[b-1]+(scfact[b-1]-1)*kap[b-1]+kap[b]
        xi1[b] = xi1[b-1]-kap[b-1]+(1-scfact[b])*kap[b]
        db = (a1[b]*a1[b]*xi2[b]+a2[b]*a2[b]*xi1[b])
        shat[b] = (a2[b]*xi1[b]*(t-tkm1-Ttilde2[b])+a1[b]*xi2[b]*(tk-t-Ttilde1[b]))/db
        lam[b] = xi1[b]*xi2[b]/db
        p1 = mystats.normpdf(a2[b]*(tk-t-Ttilde1[b]),a1[b]*(t-tkm1-Ttilde2[b]),math.sqrt(db))
        p2[b] = mystats.normprob(shat[b],math.sqrt(lam[b]),lb[b],ub[b])
        etilde[b] = p1*p2[b]

    # Normalise weights
    e = etilde/numpy.sum(etilde)
    # Draw a segment
    segidx = mystats.drawmultinom(e)
    # Draw a position in the segment
    if segidx==0:
        z = mystats.tnorm(shat[segidx],lam[segidx],skm1,ell[segidx])
    elif segidx<(nsegs-1):
        z = mystats.tnorm(shat[segidx],lam[segidx],0,ell[segidx])
    else:
        z = mystats.tnorm(shat[segidx],lam[segidx],0,sk)
    
    # Sample weight (allow for approximate prior)
    T1 = Ttilde1[segidx]+a1[segidx]*z
    T2 = Ttilde2[segidx]+a2[segidx]*z
    kap2 = 0
    for b in range(0,segidx):
        kap2 = kap2+(ub[b]-lb[b])*(ub[b]-lb[b])*kap[b]/(ell[b]*ell[b])
    kap2 = kap2+(z-lb[segidx])*(z-lb[segidx])*kap[segidx]/(ell[segidx]*ell[segidx])
    kap1 = (ub[segidx]-z)*(ub[segidx]-z)*kap[segidx]/(ell[segidx]*ell[segidx])
    for b in range(segidx+1,nsegs):
        kap1 = kap1+(ub[b]-lb[b])*(ub[b]-lb[b])*kap[b]/(ell[b]*ell[b])
    logw = lognormpdf(tk-t,T1,kap1)+lognormpdf(t-tkm1,T2,kap2)
    logw = logw-(lognormpdf(tk-t,T1,xi1[segidx])+lognormpdf(t-tkm1,T2,xi2[segidx]))

    return z, segidx+segkm1, logw
        
def eg1():
    
# Example to demonstrate interaction detection
    
    # Environment parameters
    sidelen = 100.
    
    # Measurement parameters
    sig = 5.
    g = 6
    t = 60.
    obstimes1 = [0.,20.,40.,t]
    obstimes2 = [0.,20.,40.,t]
          
    # Generate graph
    G = roadnet.makegridgraph(g,sidelen)
    
    # Vehicle parameters
    # Object 1
    p1 = [0,1,2,8,9,15,16,22] # path
    nseg1 = len(p1)-1
    # Start and end point
    s01 = 24.   
    se1 = 36.
    # Object 2
    p2 = [27,21,15,9,8,7,6,12] # path
    # Start and end point
    s02 = 55.
    se2 = 30.
    nseg2 = len(p2)-1
    
    # Statistics for path generation
    vmax = 50/3.6
    alph = 5.
    bet = 1.25
    
    
    # Meeting parameters
    segn = [4,4]    # Meeting segment for objects
    dres = 50.       # Distance along meeting segment for object 1
    xM = roadnet.getpos(G,p1,(segn[0]-1)*sidelen+dres,sidelen)
    tM = 28.
    T = 10.
    for j in range(len(obstimes1)):
        if obstimes1[j]>(T+tM):
            obstimes1[j] += T
    for j in range(len(obstimes2)):
        if obstimes2[j]>(T+tM):
            obstimes2[j] += T
    
    p1bm = p1[0:segn[0]+1]
    p1am = p1[(segn[0]-1):(nseg1+1)]
    p2bm = p2[0:segn[1]+1] 
    p2am = p2[(segn[1]-1):(nseg2+1)]
        
    # Algorithm parameters
    n = 200 # Sample size for path posterior
    n2 = 200    # sample size for position probability approximation
    K = 100 # Number of candidate paths
    # Movement statistics
    vmean = 36/3.6  
    vvar = 2.*2.
            
    # Generate data
    z1 = roadnet.makedata2(G,p1bm,p1am,xM,tM,T,obstimes1,s01,se1,sig,vmax,alph,bet) 
    z2 = roadnet.makedata2(G,p2bm,p2am,xM,tM,T,obstimes2,s02,se2,sig,vmax,alph,bet)
    
    plt.figure(1)
    plt.clf()
    roadnet.plotpath(G,[p1,p2],'-',0.)
    plt.plot(xM[0],xM[1],'s',markerfacecolor='orange',markeredgecolor='orange',markersize=9.,markeredgewidth=3.)
    plt.xlabel('x-position (m)',fontsize=16)
    plt.ylabel('y-position (m)',fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=14)
    
    # Times at which interaction is considered    
    ntpts = 51
    tax = numpy.linspace(10,60,ntpts)
    
    t0 = time.time()
    [H,P1,P2,s1,eidx1] = calcinter(G,z1,obstimes1,z2,obstimes2,n,n2,K,tax,1000,vmean,vvar,sig*sig)
    et = time.time()-t0
    
    # Plot results
    plt.figure(2)
    plt.clf()
    roadnet.plotpath(G,[p1,p2],'-',0.)
    for j in range(len(z1)):
        plt.plot(z1[j][0],z1[j][1],'bo',mew=0.,ms=8.)
    for j in range(len(z2)):
        plt.plot(z2[j][0],z2[j][1],'ro',mew=0.,ms=8.)
        
    plt.figure(3)
    plt.clf()
    Hdiag = numpy.diag(H)
    plt.plot(tax,Hdiag,linewidth=2.5)
    plt.xlabel('t (s)',fontsize=16)
    plt.ylabel('H(t,t)',fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.axis([tax[0],tax[-1],0,1])
    plt.grid(True)
    
    # Meeting place
    idxmax = numpy.argmax(Hdiag)
    plt.figure(2)
    plot_samp_pos(G,s1[idxmax],eidx1[idxmax])
    plt.plot(xM[0],xM[1],'s',markerfacecolor='orange',markeredgecolor='orange',markersize=9.,markeredgewidth=3.)
    plt.xlabel('x-position (m)',fontsize=16)
    plt.ylabel('y-position (m)',fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=14)
    
    fig = plt.figure(4)
    plt.clf()
    [tx,ty] = numpy.meshgrid(tax,tax)
    ax = fig.gca(projection='3d')
    ax.plot_surface(tx,ty,H,rstride=1,cstride=1,cmap=cm.jet)
    ax.set_xlabel('t1 (s)',fontsize=16)
    ax.set_ylabel('t2 (s)',fontsize=16)
    ax.set_zlabel('H(t1,t2)',fontsize=16)
    ax.set_zlim3d(0,1)
    for label in ax.get_xticklabels() + ax.get_yticklabels() + ax.get_zticklabels():
        label.set_fontsize(14)
    # plt.show()
    return H,tax

def real_eg1():
    
# Example to demonstrate interaction detection
    boundingboxratio = 0.1
    z1,obstimes1,z2,obstimes2,extent=testdb.ecourier_pair_data((72,'motorbike'),(44,'motorbike'),boundingboxratio)
    print "Found measurements."
    print obstimes1,obstimes2
    print z1
    print z2
    G=testdb.london_roadmap(extent)
    print "Graph generated. "
    # Environment parameters
    obstimebottom = max([obstimes1[0],obstimes2[0]])
    obstimetop = min([obstimes1[1],obstimes2[1]])
    print "path 1 length:", math.sqrt((z1[0][0]-z1[1][0])**2+(z1[1][1]-z1[0][1])**2)
    print "path 2 length:", math.sqrt((z2[0][0]-z2[1][0])**2+(z2[1][1]-z2[0][1])**2)
    print "obstimes1:",obstimes1
    print "obstimes2:",obstimes2
    # Measurement parameters
    sig = 5.
        
    # Algorithm parameters
    n = 100 # Sample size for path posterior
    n2 = 100    # sample size for position probability approximation
    K = 50 # Number of candidate paths
    # Movement statistics
    vmean = 36/3.6  
    vvar = 2.*2.
    
    # Times at which interaction is considered    
    ntpts = int(obstimetop-obstimebottom)-1
    tax = numpy.linspace(obstimebottom+1,obstimetop-1,ntpts)
    print "tax:",tax
    
    t0 = time.time()
    [H,P1,P2,s1,eidx1] = calcinter(G,z1,obstimes1,z2,obstimes2,n,n2,K,tax,1000,vmean,vvar,sig*sig)
    et = time.time()-t0
    print "Algorithm took "+str(et)+" seconds."
    return H,tax

def plot_samp_pos(G,s,eidx):
    
# plot_samp_pos(G,s,eidx)
# Plot a position on the network
# G- road network
# (s,eidx)- position (float and integer)

    elist = G.edges()
    n = len(s)
    for i in range(n):
        e = elist[eidx[i]]
        pos0 = G.node[e[0]]["pos"]
        pos1 = G.node[e[1]]["pos"]
        pos = pos0+s[i]*(pos1-pos0)/roadnet.mynorm(pos1-pos0)
        plt.plot(pos[0],pos[1],'.',color=(0,0.7,0),markersize=7.)
        
def plotres(tx,ty,H,fig):
    
# plotres(tx,ty,H,fig)
# Surface plot of the Hellinger affinity
# (tx,ty)- time axes (arrays of floats)
# H- Hellinger affinity (2D array of floats)
# fig- figure (integer)
    
    ax = fig.gca(projection='3d')
    ax.plot_surface(tx,ty,H,rstride=1,cstride=1,cmap=cm.jet)
    
    return

def lognormpdf(x,mu,kap):

# p = lognormpdf(x,mu,kap)
# Computes the log of the normal pdf
    
    p = -0.5*((x-mu)*(x-mu)/kap+math.log(2*pi*kap))
    
    return p
