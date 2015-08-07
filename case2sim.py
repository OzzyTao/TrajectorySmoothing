import roadnet
import networkx as nx
import case2
import numpy
import random
import scipy.stats
import math
import matplotlib.pyplot as plt

def sim1(nsim):

    # Environment parameters
    sidelen = 100.
    
    # Measurement parameters
    sig = 5.
    g = 5
    t = 60
    t = 40.
    obstimes1_w = [10]#[10,15,20] # way points only at this point
    obstimes2_w = [5,20]
    
    # Vehicle parameters
    # p = [0,1,7,13,14,15,21,22]
    p1 = [0,1,6,7,12,11,16]
    p2 = [13,12,7,6,1,2]
    #p = [0,1,2,7]
    vmax = 60/3.6
    alph = 5
    bet = 1.25
    
    # Algorithm parameters
    n = 100
    K = 5
    m1 = 10
    m2 = 10
    
    # Generate graph
    G = roadnet.makegridgraph(g,sidelen)
    
    # Generate candidate paths
    numcandp1 = 0
    d1 = 0
    while numcandp1<K:
        candp1 = list(nx.all_simple_paths(G,p1[0],p1[-1],d1))
        numcandp1 = len(candp1)
        d1 = d1+1
        
    numcandp2 = 0
    d2 = 0
    while numcandp2<K:
        candp2 = list(nx.all_simple_paths(G,p2[0],p2[-1],d2))
        numcandp2 = len(candp2)
        d2 = d2+1
        
        #print numcandp 
    # Put in a suitable format for processing
    segls1 = numpy.zeros((numcandp1,d1-1))
    for i in range(numcandp1):
        print candp1[i]
        for j in range(len(candp1[i])-1):
            segls1[i,j] = sidelen
    segls2 = numpy.zeros((numcandp2,d2-1))
    for i in range(numcandp2):
        print candp2[i]
        for j in range(len(candp2[i])-1):
            segls2[i,j] = sidelen
            
    # Generate data
    tact = 0
    while abs(tact-t)>0.01:
        [z1,tact] = roadnet.makedata(G,p1,obstimes1_w,sig,vmax,alph,bet)
    tact = 0
    while abs(tact-t)>0.01:
        [z2,tact] = roadnet.makedata(G,p2,obstimes2_w,sig,vmax,alph,bet)

    for cnt in range((nsim)):
        print cnt
        
        # Generate data
        tact = 0
        while abs(tact-t)>0.01:
            [z1,tact] = roadnet.makedata(G,p1,obstimes1_w,sig,vmax,alph,bet)
        tact = 0
        while abs(tact-t)>0.01:
            [z2,tact] = roadnet.makedata(G,p2,obstimes2_w,sig,vmax,alph,bet)
        
        # First stage of processing
        obstimes1 = numpy.zeros((len(obstimes1_w)+2))
        obstimes1[0] = 0.
        for i in range(len(obstimes1_w)):
            obstimes1[i+1] = obstimes1_w[i]
        obstimes1[-1] = t
        obstimes2 = numpy.zeros((len(obstimes2_w)+2))
        obstimes2[0] = 0.
        for i in range(len(obstimes2_w)):
            obstimes2[i+1] = obstimes2_w[i]
        obstimes2[-1] = t
        [datamat1,isvalid1,pP01,v1,tcum1] = case2.arclgen(obstimes1,segls1,n,vmax,alph,bet)
        [datamat2,isvalid2,pP02,v2,tcum2] = case2.arclgen(obstimes2,segls2,n,vmax,alph,bet)
        
        # Convert to positions in 2D space
        xpts1 = numpy.zeros((n*numcandp1,2*len(obstimes1)))
        for i in range((numcandp1)):
            for k in range((n)):
                cnt = k+i*n
                if isvalid1[cnt]:
                    #print candp[i]
                    for j in range(len(obstimes1)):
                        pos = roadnet.getpos(G,candp1[i],datamat1[cnt,j],sidelen)
                        xpts1[cnt,2*j] = pos[0]
                        xpts1[cnt,2*j+1] = pos[1]
        xpts2 = numpy.zeros((n*numcandp2,2*len(obstimes2)))
        for i in range((numcandp2)):
            for k in range((n)):
                cnt = k+i*n
                if isvalid2[cnt]:
                    #print candp[i]
                    for j in range(len(obstimes2)):
                        pos = roadnet.getpos(G,candp2[i],datamat2[cnt,j],sidelen)
                        xpts2[cnt,2*j] = pos[0]
                        xpts2[cnt,2*j+1] = pos[1]
                    #print datamat[cnt,:]
                    #print xpts[cnt,:]

        [pP1,w1,pP1w,ww1] = case2.calcpostprobs(z1,xpts1,isvalid1,pP01,sig)
        
        [pP2,w2,pP2w,ww2] = case2.calcpostprobs(z2,xpts2,isvalid2,pP02,sig)
        
        print pP1w
        print pP2w
                
        npts = 20
        tax = numpy.linspace(1,35,npts)
        
        bcoeff = numpy.zeros(npts)
        w1cdf = numpy.cumsum(ww1)
        w2cdf = numpy.cumsum(ww2)
        for a in range((npts)):
            #print a
            uj = random.uniform(0,1./m1)
            i = 0
            for j in range((m1)):
                while w2cdf[i]<uj:
                    i = i+1
                # Find s
                pathidx = i/n
                [snum,seg,pos] = roadnet.getseg(G,candp2[pathidx],v2[i,:],tcum2[i,:],tax[a],sidelen)
                res2 = abs((pos-G.node[candp2[pathidx][snum]]["pos"]).sum())
                #print seg
                # MC approximation of p1 and p2
                uk = random.uniform(0,1./m2)
                l1 = 0
                l2 = 0
                post1 = 0
                post2 = 0
                #print i, ww2[i]
                for k in range((m2)):
                    while w1cdf[l1]<uk:
                        l1 = l1+1
                    while w2cdf[l2]<uk:
                        l2 = l2+1
                    pathidx1 = l1/n
                    [snum1,seg1,pos1] = roadnet.getseg(G,candp1[pathidx1],v1[l1,:],tcum1[l1,:],tax[a],sidelen)
                    if (seg1[0]==seg[0]) and (seg1[1]==seg[1]):
                        if snum1>0:
                            tleft = tax[a]-tcum1[l1,snum1-1]
                        else:
                            tleft = tax[a]
                        post1 = post1+scipy.stats.beta.pdf(res2/(vmax*tleft),alph,bet)/(m2*vmax*tleft)
                    elif (seg1[0]==seg[1]) and (seg1[1]==seg[0]):
                        if snum1>0:
                            tleft = tax[a]-tcum1[l1,snum1-1]
                        else:
                            tleft = tax[a]
                        post1 = post1+scipy.stats.beta.pdf((sidelen-res2)/(vmax*tleft),alph,bet)/(m2*vmax*tleft)
                    pathidx2 = l2/n
                    [snum2,seg2,pos2] = roadnet.getseg(G,candp2[pathidx2],v2[l2,:],tcum2[l2,:],tax[a],sidelen)
                    #print seg, seg2
                    if (seg2[0]==seg[0]) and (seg2[1]==seg[1]):
                        if snum2>0:
                            tleft = tax[a]-tcum2[l2,snum2-1]
                        else:
                            tleft = tax[a]
                        post2 = post2+scipy.stats.beta.pdf(res2/(vmax*tleft),alph,bet)/(m2*vmax*tleft)
                    uk = uk+1./m2
                #print p1, p2
                if post2>1e-20:
                    bcoeff[a] = bcoeff[a]+math.sqrt(post1/post2)/m1
                uj = uj+1./m1
    
        print bcoeff
    
    return ww1, ww2
    
def sim2(nsim):
    
    # Environment parameters
    sidelen = 100.
    
    # Measurement parameters
    sig = 5.
    g = 5
    t = 45.
    obstimes1_w = [10.,16.,30.,40.]#[10,15,20] # way points only at this point
    obstimes2_w = [6.,13.,27.,36.]
    
    # Generate graph
    G = roadnet.makegridgraph(g,sidelen)
    
    # Vehicle parameters
    p1 = [0,1,6,7,12,11,16]
    nseg1 = len(p1)-1
    p2 = [13,12,7,6,1,2,3]
    nseg2 = len(p2)-1
    vmax = 60/3.6
    alph = 5.
    bet = 1.25
    
    # Meeting parameters
    segn = [3,3]    # Meeting segment for objects
    dres = 50.       # Distance along meeting segment for object 1
    xM = roadnet.getpos(G,p1,(segn[0]-1)*sidelen+dres,sidelen)
    tM = 21.
    T = 0.
    t = t+T
    
    p1bm = p1[0:segn[0]+1]
    p1am = p1[(segn[0]-1):(nseg1+1)]
    p2bm = p2[0:segn[1]+1] 
    p2am = p2[(segn[1]-1):(nseg2+1)]
    
    print p2bm
    print p1bm
    print p2am
    print p1am
    print xM
    
    # Algorithm parameters
    n = 400
    K = 10
    m1 = 2
    m2 = 2
    
    # Generate candidate paths
    numcandp1 = 0
    d1 = 0
    while numcandp1<K:
        candp1 = list(nx.all_simple_paths(G,p1[0],p1[-1],d1))
        numcandp1 = len(candp1)
        d1 = d1+1
        
    numcandp2 = 0
    d2 = 0
    while numcandp2<K:
        candp2 = list(nx.all_simple_paths(G,p2[0],p2[-1],d2))
        numcandp2 = len(candp2)
        d2 = d2+1
        
        #print numcandp 
    # Put in a suitable format for processing
    segls1 = numpy.zeros((numcandp1,d1-1))
    for i in range(numcandp1):
        print candp1[i]
        for j in range(len(candp1[i])-1):
            segls1[i,j] = sidelen
    segls2 = numpy.zeros((numcandp2,d2-1))
    for i in range(numcandp2):
        print candp2[i]
        for j in range(len(candp2[i])-1):
            segls2[i,j] = sidelen
            
    # Generate data
    random.seed(2.)
    z1 = roadnet.makedata2(G,p1bm,p1am,xM,tM,T,t,obstimes1_w,sig,vmax,alph,bet) 
    z2 = roadnet.makedata2(G,p2bm,p2am,xM,tM,T,t,obstimes2_w,sig,vmax,alph,bet)
    random.seed()
    
    roadnet.plotpath(G,p2)
    xpos = z2[0,:]
    ypos = z2[1,:]
    plt.plot(xpos,ypos,'rx')
    
    print z1
    print z2
    
    ntpts = 7
    tax = numpy.linspace(16,22,ntpts)
    bcoeff = numpy.zeros((nsim,ntpts))
    tmax = numpy.zeros((nsim))
    for cnt in range((nsim)):
        print cnt
        
        # Generate data
        #z1 = roadnet.makedata2(G,p1bm,p1am,xM,tM,T,t,obstimes1_w,sig,vmax,alph,bet) 
        #z2 = roadnet.makedata2(G,p2bm,p2am,xM,tM,T,t,obstimes2_w,sig,vmax,alph,bet)
        
        # First stage of processing (case 1 for each object)
        obstimes1 = numpy.zeros((len(obstimes1_w)+2))
        obstimes1[0] = 0.
        for i in range(len(obstimes1_w)):
            obstimes1[i+1] = obstimes1_w[i]
        obstimes1[-1] = t
        obstimes2 = numpy.zeros((len(obstimes2_w)+2))
        obstimes2[0] = 0.
        for i in range(len(obstimes2_w)):
            obstimes2[i+1] = obstimes2_w[i]
        obstimes2[-1] = t
        [datamat1,isvalid1,pP01,v1,tcum1] = case2.arclgen(obstimes1,segls1,n,vmax,alph,bet)
        [datamat2,isvalid2,pP02,v2,tcum2] = case2.arclgen(obstimes2,segls2,n,vmax,alph,bet)
        
        # Convert to positions in 2D space
        xpts1 = numpy.zeros((n*numcandp1,2*len(obstimes1)))
        for i in range((numcandp1)):
            for k in range((n)):
                idx = k+i*n
                if isvalid1[idx]:
                    #print candp[i]
                    for j in range(len(obstimes1)):
                        pos = roadnet.getpos(G,candp1[i],datamat1[idx,j],sidelen)
                        xpts1[idx,2*j] = pos[0]
                        xpts1[idx,2*j+1] = pos[1]
        print isvalid2.sum()
        xpts2 = numpy.zeros((n*numcandp2,2*len(obstimes2)))
        for i in range((numcandp2)):
            for k in range((n)):
                idx = k+i*n
                if isvalid2[idx]:
                    #print candp[i]
                    for j in range(len(obstimes2)):
                        pos = roadnet.getpos(G,candp2[i],datamat2[idx,j],sidelen)
                        if i==14:
                            plt.plot(pos[0],pos[1],'g.')
                        xpts2[idx,2*j] = pos[0]
                        xpts2[idx,2*j+1] = pos[1]
                    #print datamat[cnt,:]
                    #print xpts[cnt,:]

        [pP1,w1,pP1w,ww1] = case2.calcpostprobs(z1,xpts1,isvalid1,pP01,sig)
        
        [pP2,w2,pP2w,ww2] = case2.calcpostprobs(z2,xpts2,isvalid2,pP02,sig)
        
        print pP1w
        print pP2w
        print 1/sum(ww1**2)
        print 1/sum(ww2**2)
                
        w1cdf = numpy.cumsum(ww1)
        w2cdf = numpy.cumsum(ww2)
        maxval = -1.
        for a in range((ntpts)):
            #print a
            random.seed(3.)
            uj = random.uniform(0,1./m1)
            i = 0
            for j in range((m1)):
                while w2cdf[i]<uj:
                    i = i+1
                # Find s
                pathidx = i/n
                [snum,seg,pos] = roadnet.getseg(G,candp2[pathidx],v2[i,:],tcum2[i,:],tax[a],sidelen)
                res2 = abs((pos-G.node[candp2[pathidx][snum]]["pos"]).sum())
                print pos
                # MC approximation of p1 and p2
                random.seed(4.)
                uk = random.uniform(0,1./m2)
                l1 = 0
                l2 = 0
                post1 = 0.
                post2 = 0.
                #print i, ww2[i]
                for k in range((m2)):
                    while w1cdf[l1]<uk:
                        l1 = l1+1
                    while w2cdf[l2]<uk:
                        l2 = l2+1
                    pathidx1 = l1/n
                    [snum1,seg1,pos1] = roadnet.getseg(G,candp1[pathidx1],v1[l1,:],tcum1[l1,:],tax[a],sidelen)
                    print pos1
                    if (seg1[0]==seg[0]) and (seg1[1]==seg[1]):
                        if snum1>0:
                            tleft = tax[a]-tcum1[l1,snum1-1]
                        else:
                            tleft = tax[a]
                        post1 = post1+scipy.stats.beta.pdf(res2/(vmax*tleft),alph,bet)/(m2*vmax*tleft)
                    elif (seg1[0]==seg[1]) and (seg1[1]==seg[0]):
                        if snum1>0:
                            tleft = tax[a]-tcum1[l1,snum1-1]
                        else:
                            tleft = tax[a]
                        #print tcum1[l1,:] 
                        #print tleft
                        #print vmax*tleft
                        #print res2
                        print (sidelen-res2)
                        print tleft
                        post1 = post1+scipy.stats.beta.pdf((sidelen-res2)/(vmax*tleft),alph,bet)/(m2*vmax*tleft)
                    pathidx2 = l2/n
                    [snum2,seg2,pos2] = roadnet.getseg(G,candp2[pathidx2],v2[l2,:],tcum2[l2,:],tax[a],sidelen)
                    #print seg, seg2
                    if (seg2[0]==seg[0]) and (seg2[1]==seg[1]):
                        if snum2>0:
                            tleft = tax[a]-tcum2[l2,snum2-1]
                        else:
                            tleft = tax[a]
                        post2 = post2+scipy.stats.beta.pdf(res2/(vmax*tleft),alph,bet)/(m2*vmax*tleft)
                    uk = uk+1./m2
                #print post1, post2
                if post2>1e-20:
                    bcoeff[cnt,a] = bcoeff[cnt,a]+math.sqrt(post1/post2)/m1
                uj = uj+1./m1
            if bcoeff[cnt,a]>maxval:
                maxval = bcoeff[cnt,a]
                tmax[cnt] = tax[a]
    
        #print bcoeff
    
    return bcoeff, tmax
        