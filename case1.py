import random
import numpy
# import densityest
import math
import roadnet
import time
import matplotlib.pyplot as plt
import mystats

from scipy.stats import norm

import testdb
import networkx as nx
from profilehooks import profile

# @profile
def calcpostprobs(G,z,t,n,K,vmean,vvar,ups,tpath):

# [pp,sp,wp] = calcpostprobs(G,z,t,n,K,vmean,vvar,ups):
# Improved approximation of the posterior
# Path candidate selection is performed once per measurement
# Inputs:
# G- road network
# z- m-length list containing position measurements as 2-length lists
# t- list or array measurement times
# n- sample size
# K- no. candidate paths
# (vmean,vvar)- statistics of vehicle movement
# ups- measurement noise variance
# Outputs:
# pp- sampled paths (list of sequences of node indices)
# sp- (estimated) distance of object along sampled paths (list of floats)
# wp- sample weights    (list of floats)

    m = len(z)
    
    Tp = roadnet.create_turn_tab(G)
     
    vmax = vmean+3*math.sqrt(vvar)
        
    logwt = numpy.zeros(n)
    psamp = []
    ssamp = numpy.zeros((n,m))
    
    # First measurement
    # Find all edges on which the target could lie
    [epos,eidx] = roadnet.find_start_paths(G,z[0],ups)
    # print "epos:",epos
    # print "eidx:",eidx
    # Sample edges (both directions are considered)

    pst = []
    pst_rev = []
    psamp = []
    for i in range(n):
        # Draw edge and distance along edge
        [pidx_samp,ssamp[i][0]] = samp_start_path(z[0],epos,ups)
        psamp.append(eidx[pidx_samp])
        if psamp[i][-2] not in pst:
            pst.append(psamp[i][-2])
        if psamp[i][-1] not in pst_rev:
            pst_rev.append(psamp[i][-1])

    # Second measurement
    # Select candidate paths (for both directions)
    [pe,ee] = roadnet.findnodes(G,z[1],ups)
    [cpath_pos,cpath_idx,turnpen] = roadnet.find_paths_by_IDs(G,Tp,pe,ee,pst,K,tpath)

    ncand = len(cpath_idx)

    [cpath_pos1,cpath_idx1,turnpen1] = roadnet.find_paths_by_IDs(G,Tp,pe,ee,pst_rev,K,tpath)

    ncand1 = len(cpath_idx1)
    pst = []
    for i in range(n):
        ncandi = 0
        cpath_idxi = []
        cpath_posi = []
        tpeni = []
        isrev = []
        s0 = []
        # Find candidate paths matching this path sample
        for a in range(ncand):
            if cpath_idx[a][0]==psamp[i][-2] and cpath_idx[a][1]==psamp[i][-1]:
                cpath_posi.append(cpath_pos[a])
                cpath_idxi.append(cpath_idx[a])
                tpeni.append(turnpen[a])
                isrev.append(0)
                s0.append(ssamp[i][0])
                ncandi = ncandi+1
        for a in range(ncand1):
            if cpath_idx1[a][0]==psamp[i][-1] and cpath_idx1[a][1]==psamp[i][-2]:
                cpath_posi.append(cpath_pos1[a])
                cpath_idxi.append(cpath_idx1[a])
                tpeni.append(turnpen1[a])
                isrev.append(1)
                ell = roadnet.mynorm(cpath_pos1[a][:,1]-cpath_pos1[a][:,0])
                s0.append(ell-ssamp[i][0])
        # Draw a sample path if there are any matching candidates
        if len(cpath_idxi)>0:
            [pidx_samp,seg_idx,ssamp[i][1],lwtj,isvalid] = samp_path(z[1],t[1]-t[0],cpath_posi,tpeni,s0,ups,vmean,vvar)
            psamp_new = cpath_idxi[pidx_samp][2:seg_idx+2]
            if isrev[pidx_samp]:
                psamp[i] = [psamp[i][1],psamp[i][0]]
            psamp[i] = psamp[i]+psamp_new
            # print "psamp:",psamp[i]
            if psamp[i][-2] not in pst:
                pst.append(psamp[i][-2])
            logwt[i] = lwtj
        # Set the weight to zero (effectively) if there are no matching candidates
        else:
            logwt[i] = -1e10
    
    # Process the remaining measurements
    for j in range(2,m):
        # Select candidate paths
        [pe,ee] = roadnet.findnodes(G,z[j],ups)
        [cpath_pos,cpath_idx,turnpen] = roadnet.find_paths_by_IDs(G,Tp,pe,ee,pst,K,[])
        ncand = len(cpath_idx)
        pst = []
        for i in range(n):
            cpath_posi = []
            cpath_idxi = []
            tpeni = []
            s0 = []
            # Find candidate paths matching this path sample
            for a in range(ncand):
                if cpath_idx[a][0]==psamp[i][-2] and cpath_idx[a][1]==psamp[i][-1]:
                    cpath_posi.append(cpath_pos[a])
                    cpath_idxi.append(cpath_idx[a])
                    tpeni.append(turnpen[a])
                    s0.append(ssamp[i][j-1])
            # Draw a sample path if there are any matching candidates
            if len(cpath_idxi)>0:
                [pidx_samp,seg_idx,ssamp[i][j],lwtj,isvalid] =samp_path(z[j],t[j]-t[j-1],cpath_posi,tpeni,s0,ups,vmean,vvar)
                psamp_new = cpath_idxi[pidx_samp][2:seg_idx+2]
                psamp[i] = psamp[i]+psamp_new
                if psamp[i][-2] not in pst:
                    pst.append(psamp[i][-2])
                logwt[i] += lwtj
            # Set the weight to zero (effectively) if there are no matching candidates
            else:
                logwt[i] = -1e10
            
    # Normalise weights
    maxwt = -1e20
    for i in range(n):
        if logwt[i]>maxwt:
            maxwt = logwt[i]
    wtilde = numpy.zeros((n))
    for i in range(n):
        wtilde[i] = math.exp(logwt[i]-maxwt)
    wt = wtilde/numpy.sum(wtilde)
    
    # Find the unique set of sampled paths
    npaths = 0
    pout = []
    seout = []
    wout = numpy.zeros((n))
    nout = 0
    for i in range(n):
        isin = 0
        if psamp[i] in pout:
            idx = pout.index(psamp[i])
            seout[idx] += wt[i]*ssamp[i][-1]
            wout[idx] += wt[i]
        else:
            pout.append(psamp[i])
            seout.append(wt[i]*ssamp[i][-1])
            wout[nout] = wt[i]
            nout += 1
            
    wout = wout[0:nout]
    seout = seout[0:nout]
    for i in range(nout):
        seout[i] /= wout[i]
        
    return pout, seout, wout


def calcpostprobs_case2(G,z,t,n,K,vmean,vvar,ups):

# [pp,ppos,wp,seg_idx,ssamp] = calcpostprobs_case2(G,z,t,n,vmean,vvar,ups):
# Improved approximation of the posterior
# Path candidate selection is performed once per measurement
# Sample positions are returned in addition to paths
# Inputs:
# G- road network
# z- m-length list containing position measurements as 2-length lists
# t- list or array of measurement times
# n- sample size
# (vmean,vvar)- statistics of vehicle movement
# ups- measurement noise variance
# Outputs:
# pp- sampled paths with node indices (list of sequences of node indices)
# ppos- samped paths with node positions (list of list of node positions)
# wp- sample weights (list of floats)
# seg_idx- index (list of lists of integers)
# ssamp- position along segment of each sampled path (m x n array of floats)

    m = len(z)
    seg_idx = [[0 for i in range(n)] for j in range(m)]
    
    Tp = roadnet.create_turn_tab(G)
     
    vmax = vmean+3*math.sqrt(vvar)
        
    logwt = numpy.zeros(n)
    validsamp = [1 for k in range(n)]
    psamp = []
    ssamp = numpy.zeros((m,n))
    
    # First measurement
    # Find all edges on which the target could lie
    [epos,eidx] = roadnet.find_start_paths(G,z[0],ups)
    # Sample edges
    pst = []
    psamp = []
    ppos = []
    for i in range(n):
        # Draw path and distance along path
        [pidx_samp,ssamp[0,i]] = samp_start_path(z[0],epos,ups)
        psamp.append(eidx[pidx_samp])
        if psamp[i][-2] not in pst:
            pst.append(psamp[i][-2])
        if psamp[i][-1] not in pst:
            pst.append(psamp[i][-1])
        seg_idx[0][i] = 0
        
    # Second measurement
    # Select candidate paths (for both directions)
    [pe,ee] = roadnet.findnodes(G,z[1],ups)
    [cpath_pos,cpath_idx,turnpen] = roadnet.find_paths_by_IDs(G,Tp,pe,ee,pst,K)
    ncand = len(cpath_idx)
    pst = []
    for i in range(n):
        ncandi = 0
        cpath_idxi = []
        cpath_posi = []
        tpeni = []
        isrev = []
        s0 = []
        # Find candidate paths matching this path sample
        for a in range(ncand):
            if cpath_idx[a][0]==psamp[i][-2] and cpath_idx[a][1]==psamp[i][-1]:
                cpath_posi.append(cpath_pos[a])
                cpath_idxi.append(cpath_idx[a])
                tpeni.append(turnpen[a])
                isrev.append(0)
                s0.append(ssamp[0,i])
            elif cpath_idx[a][1]==psamp[i][-2] and cpath_idx[a][0]==psamp[i][-1]:
                cpath_posi.append(cpath_pos[a])
                cpath_idxi.append(cpath_idx[a])
                tpeni.append(turnpen[a])
                isrev.append(1)
                ell = roadnet.mynorm(cpath_pos[a][:,1]-cpath_pos[a][:,0])
                s0.append(ell-ssamp[0,i])
        
        # Draw a sample path if there are any matching candidates
        # print "cpath_idxi:", cpath_idxi

        if len(cpath_idxi)>0:
            [pidx_samp,seg_samp,ssamp[1,i],lwtj,validsamp[i]] = samp_path(z[1],t[1]-t[0],cpath_posi,tpeni,s0,ups,vmean,vvar)
            if validsamp[i]:
                psamp_new = cpath_idxi[pidx_samp][2:seg_samp+2]
                if isrev[pidx_samp]:
                    psamp[i] = [psamp[i][1],psamp[i][0]]
                    ell = roadnet.mynorm(cpath_pos[pidx_samp][:,1]-cpath_pos[pidx_samp][:,0])
                    ssamp[0,i] = ell-ssamp[0,i]
                psamp[i] = psamp[i]+psamp_new
                seg_idx[1][i] = seg_samp
                if psamp[i][-2] not in pst:
                    pst.append(psamp[i][-2])
                logwt[i] = lwtj
            # Set the weight to zero (effectively) if there are no matching candidates
            else:
                logwt[i] = -1e10
        else:
            logwt[i] = -1e10
    # Processing remaining measurements
    for j in range(2,m):
        # Select candidate paths
        [pe,ee] = roadnet.findnodes(G,z[j],ups)
        [cpath_pos,cpath_idx,turnpen] = roadnet.find_paths_by_IDs(G,Tp,pe,ee,pst,K)
        ncand = len(cpath_idx)
        pst = []
        for i in range(n):
            if validsamp[i]:
                cpath_posi = []
                cpath_idxi = []
                tpeni = []
                s0 = []
                # Find candidate paths matching this path sample
                for a in range(ncand):
                    if cpath_idx[a][0]==psamp[i][-2] and cpath_idx[a][1]==psamp[i][-1]:
                        cpath_posi.append(cpath_pos[a])
                        cpath_idxi.append(cpath_idx[a])
                        tpeni.append(turnpen[a])
                        s0.append(ssamp[j-1,i])
                # Draw a sample path if there are any matching candidates
                if len(cpath_idxi)>0:
                    [pidx_samp,seg_samp,ssamp[j,i],lwtj,validsamp[i]] =               samp_path(z[j],t[j]-t[j-1],cpath_posi,tpeni,s0,ups,vmean,vvar)
                    if validsamp[i]:
                        seg_idx[j][i] = seg_samp+len(psamp[i])-2
                        psamp_new = cpath_idxi[pidx_samp][2:seg_samp+2]
                        psamp[i] = psamp[i]+psamp_new
                        if psamp[i][-2] not in pst:
                            pst.append(psamp[i][-2])
                        logwt[i] += lwtj
                    # Set the weight to zero (effectively) if there are no matching candidates
                    else:
                        logwt[i] = -1e10
                else:
                    logwt[i] = -1e10
    
    # Normalise weights
    maxwt = -1e20
    for i in range(n):
        if logwt[i]>maxwt:
            maxwt = logwt[i]
        
    wtilde = numpy.zeros((n))
    for i in range(n):
        wtilde[i] = math.exp(logwt[i]-maxwt)
    wt = wtilde/numpy.sum(wtilde)
    
    ppos = []
    for i in range(n):
        nseg = len(psamp[i])-1
        ppos_samp = numpy.zeros((2,nseg+1))
        for j in range(nseg+1):
            pos = G.node[psamp[i][j]]["pos"]
            ppos_samp[0,j] = pos[0]
            ppos_samp[1,j] = pos[1]
        ppos.append(ppos_samp)
        
    return psamp, ppos, wt, seg_idx, ssamp
                
def samp_start_path(z1,edge_pos,v):

# [edgeidx,z] = samp_start_path(z1,edge_pos,v)
# Samples a starting edge
# Inputs:
# z1- 2-element measurement list or array
# edge_pos- candidate edges (a list of arrays containing the positions of the edge nodes)
# v- measurement noise variance
# Outputs:
# edgeidx- index of sampled edge (integer)
# z- distance along edge (float)

    nedge = len(edge_pos)
    
    x1 = z1[0]
    y1 = z1[1]
    sqv = math.sqrt(v)
    
    mui = numpy.zeros((nedge))
    elli = numpy.zeros((nedge))
    qrho = numpy.zeros((nedge))
    # Compute weights for each edge
    for i in range(nedge):
        xi = edge_pos[i][0,0]
        yi = edge_pos[i][1,0]
        xe = edge_pos[i][0,1]
        ye = edge_pos[i][1,1]
        thi = math.atan2(ye-yi,xe-xi)
        cth = math.cos(thi)
        sth = math.sin(thi)
        elli[i] = math.sqrt((xe-xi)**2+(ye-yi)**2)
        mui[i] = (x1-xi)*cth+(y1-yi)*sth
        P = mystats.normprob(mui[i],sqv,0,elli[i])
        if abs(cth)>1e-10:
            npdf = mystats.normpdf(y1-yi,(x1-xi)*sth/cth,sqv/abs(cth))/abs(cth)
        else:
            npdf = mystats.normpdf(x1-xi,(y1-yi)*cth/sth,sqv/abs(sth))/abs(sth)
        qrho[i] = P*npdf/elli[i]
    # Normalise weights
    w = sum(qrho)
    qrho = qrho/w
    # Draw an edge
    idx = mystats.drawmultinom(qrho)
    # Draw the distance along the edge
    z = mystats.tnorm(mui[idx],v,0,elli[idx])

    return idx, z

def samp_path(z,t,cand_paths,turnpen,s0,v,vmean,vvar):

# [pathidx,segidx,z,logwt] = samp_path(z,t,cand_paths,turnpen,s0,v,vmean,vvar)
# Samples a path extension
# Inputs:
# z- 2-element measurement list or array of floats
# t- duration since last measurement 
# cand_paths- candidate paths (a list of arrays containing the positions of each node on the path)
# turnpen- turn penalty along each path
# s0- starting position in final segment of each candidate path
# v- measurement noise variance
# (vmean,vvar)- statistics of vehicle movement
# Outputs:
# pathidx- index of sampled path (integer)
# segidx- segment along path in which vehicle lies (integer)
# z- distance along segment (float)
# logwt- log of sample weight (float)

    npaths = len(cand_paths)
    sqv = math.sqrt(v)
    eps = numpy.spacing(1)
    
    # Pre-allocate
    x = z[0]
    y = z[1]
    gamma = numpy.zeros((npaths))
    C1 = numpy.zeros((npaths))
    qpath = numpy.zeros((npaths))
    pathprior = numpy.zeros((npaths))
    b = list(numpy.zeros((npaths)))
    nu = list(numpy.zeros((npaths)))
    zeta = list(numpy.zeros((npaths)))
    ell = list(numpy.zeros((npaths)))
    nsegs = [0 for i in range(npaths)]
    lb = [0 for i in range(npaths)]
    # Compute the weight of each path
    for i in range(npaths):
        sz_path = cand_paths[i].shape
        nsegs[i] = sz_path[1]-1
        if nsegs[i]==1:
            lb[i] = s0[i]
        else:
            lb[i] = 0.
        ell[i] = numpy.zeros(nsegs[i])
        # Calculate prior mean and variance
        # First segment (account for s0)
        pos1 = cand_paths[i][:,0]
        pos2 = cand_paths[i][:,1]
        ell[i][0] = numpy.linalg.norm(pos1-pos2)
        Tbarj = ell[i][0]/vmean
        kapj = ell[i][0]*ell[i][0]*vvar/(vmean**4)
        Ttilde = -s0[i]*Tbarj/ell[i][0]
        kaptilde = 0.
        c = (ell[i][0]-s0[i])/ell[i][0]
        mean_term = Tbarj
        var_term = c*c*kapj
        for j in range(1,nsegs[i]):
            pos1 = cand_paths[i][:,j]
            pos2 = cand_paths[i][:,j+1]
            ell[i][j] = numpy.linalg.norm(pos1-pos2)
            Ttilde += mean_term
            kaptilde += var_term
            Tbarj = ell[i][j]/vmean
            kapj = ell[i][j]*ell[i][j]*vvar/(vmean**4)
            mean_term = Tbarj
            var_term = kapj
        # Compute location sampling density and path weight
        xi = pos1[0]
        yi = pos1[1]
        xe = pos2[0]
        ye = pos2[1]
        thi = math.atan2(ye-yi,xe-xi)
        cth = math.cos(thi)
        sth = math.sin(thi)
        mu = (x-xi)*cth+(y-yi)*sth
        alph = ell[i][nsegs[i]-1]*(t-Ttilde)/Tbarj
        lam = ell[i][nsegs[i]-1]*ell[i][nsegs[i]-1]*(kaptilde+kapj)/(Tbarj*Tbarj)
        nu[i] = (v*alph+lam*mu)/(v+lam)
        zeta[i] = lam*v/(v+lam)
        if abs(cth)>1e-10:
            npdf1 = mystats.normpdf(y-yi,(x-xi)*sth/cth,sqv/abs(cth))/abs(cth)
        else:
            npdf1 = mystats.normpdf(x-xi,(y-yi)*cth/sth,sqv/abs(sth))/abs(sth)
        npdf2 = mystats.normpdf(mu,alph,math.sqrt(v+lam))
        sqzeta = math.sqrt(zeta[i])
        C2 = mystats.normprob(nu[i],sqzeta,lb[i],ell[i][nsegs[i]-1])
        sqlam = math.sqrt(lam)
        C1[i] = mystats.normprob(alph,sqlam,lb[i],ell[i][nsegs[i]-1])
        pathprior[i] = 1/(numpy.sum(ell[i])+turnpen[i])
        if npdf2>eps:
            gamma[i] = npdf1*npdf2*C2/C1[i]
        else:
            gamma[i] = 0.
    # Normalise path weights
    pathprior = pathprior/numpy.sum(pathprior)
    qpath = gamma*pathprior
    wt = numpy.sum(qpath)
    if wt<eps:
        # Ignore zero weighted paths
        pidx = 0
        segidx = 0
        ssamp = 0.
        logwt = -1e10
        isvalid = 0
        return pidx, segidx, ssamp, logwt, isvalid
    # If we get here we can sample a valid path
    isvalid = 1
    qpath = qpath/wt
    # Sample a path
    pidx = mystats.drawmultinom(qpath)
    segidx = nsegs[pidx]-1
    # Now sample a position
    if nsegs[pidx]==1:
        ssamp = mystats.tnorm(nu[pidx],zeta[pidx],s0[pidx],ell[pidx][nsegs[pidx]-1])
    else:
        ssamp = mystats.tnorm(nu[pidx],zeta[pidx],0,ell[pidx][nsegs[pidx]-1])
    # Calculate sample weight
    logwt = math.log(wt)+mystats.calc_samp_wt(t,ell[pidx],s0[pidx],ssamp,lb[pidx],C1[pidx],vmean,vvar)
    

    return pidx, segidx, ssamp, logwt, isvalid
            
def eg1():

# Set up for demonstration of case 1 algorithm
# This returns the path posterior probabilities for the given example

    # Environment parameters
    sidelen = 100
    
    # Measurement parameters
    sig = 0.5    # measurement standard deviation
    g = 5       # network grid size
    t = 45.     # surveillance duration
    s0 = 20.    # starting location along first edge
    se = 45.    # ending location along final edge
    obstimes = [0.,22.,t]  # observation times
        
    # Vehicle parameters
    p = [0,1,2,7,12,17,18]  # path nodes
    
    # Parameters used for path generation
    vmax = 60./3.6
    a = 5.
    b = 1.25
    
    # Algorithm parameters
    n = 100  # number of samples
    K = 100  # minimum number of candidate paths
    vmean = 10. # mean velocity
    vvar = 4. # velocity variance
    
    # Generate graph
    G = roadnet.makegridgraph(g,sidelen)
            
    # Generate data
    tact = 0
    while abs(tact-t)>0.01:
        [z,tact] = roadnet.makedata(G,p,obstimes,s0,se,sig,vmax,a,b)
        
    t0 = time.time()
    [pp,sep,wp] = calcpostprobs(G,z,obstimes,n,K,vmean,vvar,sig*sig)
    et = time.time()-t0
    print et
    
    # Plot best paths
    np = len(pp)
    idxsort = numpy.argsort(wp)            
    nplot = min([3,np])
    plt.figure(1)
    roadnet.plotpath(G,[p,pp[idxsort[np-1]],pp[idxsort[np-2]]],'-',0.)
    for j in range(len(z)):
        plt.plot(z[j][0],z[j][1],'kx',markersize=10.,markeredgewidth=2)
    plt.show()
    return pp, wp

def test_UTM():
    boundingboxbuffer_ratio = 2.0 
    z,obstimes,extent = testdb.ecourier_data(8,2,boundingboxbuffer_ratio)
    G=testdb.london_roadmap(extent)
    x0,y0 = z[0]
    dis=[]
    tmpx=G.node[315161]['pos'][0]
    tmpy=G.node[315161]['pos'][1]
    dis.append(math.sqrt((x0-tmpx)**2+(y0-tmpy)**2))
    print 'Min distance:',min(dis)

def real_eg1():
    sig = 5
    # Algorithm parameters
    n = 100  # number of samples
    K = 100  # minimum number of candidate paths
    vmean = 10. # mean velocity
    vvar = 4. # velocity variance
    
    boundingboxbuffer_ratio = 2.0 
    z,obstimes,extent, endID = testdb.ecourier_vertex_data(8,boundingboxbuffer_ratio)
    print "measurements:"
    print z
    print "observation times"
    print obstimes
    G=testdb.london_roadmap(extent)

    #plot
    # position={}
    # for node in G:
    #   position[node]=G.node[node]['pos']
    # nx.draw(G,pos=position,node_size=50)
    # for point in z:
    #     plt.plot(point[0],point[1],'kx',markersize=10.0,markeredgewidth=2)
    # plt.show()

    t0 = time.time()
    [pp,sep,wp] = calcpostprobs(G,z,obstimes,n,K,vmean,vvar,sig*sig)
    et = time.time()-t0
    print et, "secs, done algorithm"
    print endID

    return pp, sep, wp, endID

def sim1(g,nsim,t,m,n):

# [rmse,rhist] = sim1(g,nsim,t,m,n):
# Simulation for performance assessment using random paths
# Inputs:
# g- size of grid network
# nsim- number of simuations
# t- surveillance duration
# m- number of equispaced (in time) position measurements
# n- sample size for algorithm
# Outputs:
# rmse: RMS error between the true path and the MAP path (spatial)
# rhist: histogram of rank of true path
     
    # Environment parameters
    sidelen = 100
    
    # Measurement parameters
    sig = 5.    # measurement standard deviation
    s0 = 30.    # starting location along first edge
    obstimes = [0 for j in range(m)]
    for j in range(m):
        obstimes[j] = j*t/(m-1)
        
    # Vehicle parameters
    vmax = 50./3.6
    a = 5.
    b = 1.25
    vmean = vmax*a/(a+b)

    # Algorithm parameters
    K = 100  # minimum number of candidate paths
    vmean = 10. # mean velocity
    vvar = 4.   # velocity variance
    
    # Generate graph
    G = roadnet.makegridgraph(g,sidelen)
    
    rmse = [0. for j in range(len(n))]
    rhist = [[0 for i in range(n[j])] for j in range(len(n))]

    for i in range(nsim):
        # Generate random path and data
        [p,se,z] = mystats.gen_rand_path(G,t,s0,vmax,a,b,obstimes,sig)
            
        for j in range(len(n)):
            print i, j
            t0 = time.time()
            [pp,sep,wp] = calcpostprobs(G,z,obstimes,n[j],K,vmean,vvar,sig*sig)
            et = time.time()-t0
        
            # Generate output statistics
            np = len(pp)
            idxsort = numpy.argsort(wp)
            idxmax = numpy.argmax(wp)
            rmse_samp = [0 for k in range(np)]
            for k, psamp in enumerate(pp):
                rmse_samp[k] = mystats.calc_path_rmse(G,s0,p,se,psamp,sep[k])
            idxminrmse = numpy.argmin(rmse_samp)
            rhist[j][idxsort[idxminrmse]] += 1./nsim
            rmse[j] = rmse[j]+rmse_samp[idxmax]/nsim

    return rmse, rhist
                    
    
            