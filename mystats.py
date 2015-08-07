import numpy
import math
import random
from scipy.stats import norm

def calcC(delt,Ttilde,kaptilde,Tbar,kap,lb,ell):

# C = calcC(delt,Ttilde,kaptilde,Tbar,kap,lb,ell)
# Computes the normalising constant for true position prior using Gauss-Hermite 
# quadrature
# Inputs:
# delt- elapsed time (float)
# Ttilde- mean duration (float)
# kaptilde- duration variance (float)
# Tbar- mean duration for last segment (float)
# kap- duration variance for last segment (float)
# lb- lower bound (float)
# ell- segment length (float)

    eps = numpy.spacing(1)
    
    C = 0.
    # Number of quadrature points
    n = 9
    mu = delt-Ttilde
    sig = math.sqrt(kaptilde)

    [z,w] = numpy.polynomial.hermite.hermgauss(n)
    w = w/sum(w)
    
    xp = (math.sqrt(2*kap)*z+Tbar)/ell
        
    C = 0.
    for i in range(n):
        C += w[i]*(normcdf(ell*xp[i],mu,sig)-normcdf(lb*xp[i],mu,sig))/xp[i];
        
    return C

def calc_samp_wt(t,ell,s0,ssamp,lb,C1,vmean,vvar):
  
# logwt = calc_samp_wt(t,ell,s0,ssamp,lb,C1,vmean,vvar)
# Computes the sample weight accounting for the position prior approximation
# Inputs:
# t- elapsed time (float)
# ell- segment lengths (list of float)
# s0- starting position (float)
# ssamp- sampled position (float)
# lb- lower bound (float)
# C1- normalising constant for position prior approximation (float)
# (vmean,vvar)- vehicle motion statistics
# Output:
# logwt- log of sample weight

    eps = numpy.spacing(1)
    
    # Compute statistics for the position prior and its approximation
    nsegs = len(ell)
    Tbarj = ell[0]/vmean
    kapj = ell[0]*ell[0]*vvar/(vmean**4)
    Ttilde = -s0*Tbarj/ell[0]
    kaptilde = eps
    c = (ell[0]-s0)/ell[0]
    mean_term = Tbarj
    var_term = c*c*kapj
    for j in range(1,nsegs):
        Ttilde += mean_term
        kaptilde += var_term
        Tbarj = ell[j]/vmean
        kapj = ell[j]*ell[j]*vvar/(vmean**4)
        mean_term = Tbarj
        var_term = kapj
    kapwt = kaptilde+ssamp*ssamp*kapj/(ell[-1]*ell[-1])
    kapq = kaptilde+kapj
    ept = t-(Ttilde+ssamp*Tbarj/ell[-1])
    C = calcC(t,Ttilde,kaptilde,Tbarj,kapj,lb,ell[-1])
    logwt = (math.log(kapq/kapwt)+ept*ept*(1/kapq-1/kapwt))/2-math.log(C+eps)+math.log(ell[-1]*C1/Tbarj)
    
    
    return logwt
    
def normcdf(x,mu,sigma):

# P = normcdf(x,mu,sigma)
# Computes an approximation to the noraml CDF

    a = [1.253313402,-0.9999673043,0.6262972801,-0.3316218430,0.1522723563,-5.982834993e-02,1.915649350e-02,-4.644960579e-03,7.771088713e-04,-7.830823677e-05,3.534244658e-06]
    q = len(a)-1
    z = (x-mu)/sigma
    if z>0:
        v = 1.
        P = a[0]
        for i in range(q):
            v = v*z
            P = P+v*a[i+1]
        return 1-P*normpdf(z,0,1)
    else:
        y = -z
        v = 1.
        P = a[0]
        for i in range(q):
            v = v*y
            P = P+v*a[i+1]
        return P*normpdf(z,0,1)
    
def normpdf(x,mu,sigma):

# p = normpdf(x,mu,sigma)
# Computes the normal PDF.

    ep = (x-mu)/sigma
    p = math.exp(-ep*ep/2)/(math.sqrt(2*math.pi)*sigma)
    
    return p

def normprob(mu,sigma,lb,ub):

# P = normprob(mu,sigma,lb,ub)
# Computes the normal probability between (lb,ub)

    P = normcdf(ub,mu,sigma)-normcdf(lb,mu,sigma)
    
    return P


def gen_rand_path(G,t,s0,vmax,a,b,obstimes,sig):

# [p,se,z] = gen_rand_path(G,t,s0,vmax,a,b,obstimes,sig)
# Generates a random path and measurements
# Inputs:
# G- road network
# t- path duration (float)
# s0- starting position on first edge (float)
# (vmax,a,b)- statistics for path generation (floats)
# obstimes- 0\observation times (list of floats)
# measurement noise standard deviation (float)
# Output:
# p- path (list of integers)
# se- ending position on last edge (float)
# z- position measurements (list of list of floats)

    nnodes = len(G.nodes())
    npts = len(obstimes)
    
    isvalid = 0
    while not isvalid:
        p = []
        v = []
        td = []
        # Select an initial node at random
        p.append(int(math.ceil(random.uniform(0,1)*nnodes)-1))
        ni = G.neighbors(p[0])
        # Select a neighbour at random
        idx = int(math.ceil(random.uniform(0,1)*len(ni))-1)
        p.append(ni[idx])
        ell = mynorm(G.node[p[1]]["pos"]-G.node[p[0]]["pos"])
        v.append(vmax*random.betavariate(a,b))
        td.append((ell-s0)/v[-1])
        tel = (ell-s0)/v[-1]
        # Check for end of interval
        if tel>t:
            se = v*t+s0
            td[-1] = t
            islt = 0
        else:
            islt = 1
        
        while islt: # Not at the end of the interval
            # Add a new edge
            ni = G.neighbors(p[-1])
            lni = len(ni)
            P = numpy.zeros(lni)
            pos0 = G.node[p[-2]]["pos"]
            pos1 = G.node[p[-1]]["pos"]
            # Select a new node (favour straight paths)
            for j in range(lni):
                pos2 = G.node[ni[j]]["pos"]
                ep1 = pos1-pos0
                ep2 = pos2-pos1
                th = math.acos(mydot(ep2,ep1)/(mynorm(ep1)*mynorm(ep2)))
                P[j] = 1/(1.+1*abs(th)**4)
            P = P/numpy.sum(P)
            idx = drawmultinom(P)
            p.append(ni[idx])
            v.append(vmax*random.betavariate(a,b))
            ell = mynorm(G.node[p[-1]]["pos"]-G.node[p[-2]]["pos"])
            tel += ell/v[-1]
            td.append(ell/v[-1])
            # Check for end of interval
            if tel>t:
                se = ell-v[-1]*(tel-t)
                td[-1] = se/v[-1]
                islt = 0
        # Path is valid if it doesn't contain loops
        isvalid = (len(p)==len(set(p)))
    
    tcum = numpy.cumsum(td)+obstimes[0]
    nseg = len(p)-1
    
    #Generate measurements
    z = [[] for k in range(npts)]
    
    pos0 = G.node[p[0]]["pos"]
    pos1 = G.node[p[1]]["pos"]
    ell = numpy.linalg.norm(pos1-pos0)
    pos = pos0+(pos1-pos0)*s0/ell
    z[0] = pos+sig*numpy.random.randn(2)
    dprev = s0/ell
    jprev = 0
    for k in range(npts-2):
        j = (obstimes[k+1]>tcum).sum()
        if j<nseg:
            eweight = G[p[j]][p[j+1]]['weight']
            if (j != jprev):
                df = (obstimes[k+1]-tcum[j-1])*v[j]/eweight
            else:
                df = dprev+(obstimes[k+1]-obstimes[k])*v[j]/eweight
            pos0 = G.node[p[j]]["pos"]
            pos1 = G.node[p[j+1]]["pos"]
            pos = pos0+(pos1-pos0)*df
            z[k+1] = pos+sig*numpy.random.randn(2)
            dprev = df
            jprev = j
    
    pos0 = G.node[p[nseg-1]]["pos"]
    pos1 = G.node[p[nseg]]["pos"]
    ell = numpy.linalg.norm(pos1-pos0)
    pos = pos0+(pos1-pos0)*se/ell
    z[npts-1] = pos+sig*numpy.random.randn(2)
    
    return p, se, z

def gen_rand_path2(G,edel,t,s0,vmax,a,b,obstimes,sig):

# [p,se,z] = gen_rand_path(G,edel,t,s0,vmax,a,b,obstimes,sig)
# Generates a random path and measurements
# Paths must avoid certain edges
# Inputs:
# G- road network
# edel- edges which are to be avoided (list of integers)
# t- path duration (float)
# s0- starting position on first edge (float)
# (vmax,a,b)- statistics for path generation (floats)
# obstimes- 0\observation times (list of floats)
# measurement noise standard deviation (float)
# Output:
# p- path (list of integers)
# se- ending position on last edge (float)
# z- position measurements (list of list of floats)

    nnodes = len(G.nodes())
    npts = len(obstimes)
    nedel = len(edel)
    elist = G.edges()
    eps = numpy.spacing(1)
    
    isvalid = 0
    while not isvalid:
        p = []
        v = []
        td = []
        isvalid_start = 0
        while not isvalid_start:
            p = []
            p.append(int(math.ceil(random.uniform(0,1)*nnodes)-1))
            ni = G.neighbors(p[0])
            lni = len(ni)
            P = numpy.zeros(lni)
            for j in range(lni):
                P[j] = 1.
                ej = tuple([p[0],ni[j]])
                if ej in elist:
                    if elist.index(ej) in edel:
                        P[j] = 0.
                else:
                    ej = tuple([ni[j],p[0]])
                    if elist.index(ej) in edel:
                        P[j] = 0.
            sumP = numpy.sum(P)
            if sumP>eps:
                P = P/sumP
                isvalid_start = 1
            else:
                isvalid_start = 0
        idx = drawmultinom(P)
        p.append(ni[idx])
        ell = mynorm(G.node[p[1]]["pos"]-G.node[p[0]]["pos"])
        v.append(vmax*random.betavariate(a,b))
        td.append((ell-s0)/v[-1])
        tel = (ell-s0)/v[-1]
        if tel>t:
            se = v[-1]*t+s0
            td[-1] = t
            islt = 0
        else:
            islt = 1
    
    
        while islt:
            # Add a new edge
            ni = G.neighbors(p[-1])
            lni = len(ni)
            P = numpy.zeros(lni)
            pos0 = G.node[p[-2]]["pos"]
            pos1 = G.node[p[-1]]["pos"]
            for j in range(lni):
                pos2 = G.node[ni[j]]["pos"]
                ep1 = pos1-pos0
                ep2 = pos2-pos1
                th = math.acos(mydot(ep2,ep1)/(mynorm(ep1)*mynorm(ep2)))
                P[j] = 1/(1.+1*abs(th)**4)
                # Set to zero if the edge is missing
                ej = tuple([p[-1],ni[j]])
                if ej in elist:
                    if elist.index(ej) in edel:
                        P[j] = 0.
                else:
                    ej = tuple([ni[j],p[-1]])
                    if elist.index(ej) in edel:
                        P[j] = 0.
            P = P/numpy.sum(P)
            idx = drawmultinom(P)
            p.append(ni[idx])
            v.append(vmax*random.betavariate(a,b))
            ell = mynorm(G.node[p[-1]]["pos"]-G.node[p[-2]]["pos"])
            tel += ell/v[-1]
            td.append(ell/v[-1])
            if tel>t:
                se = ell-v[-1]*(tel-t)
                td[-1] = se/v[-1]
                islt = 0
        # Path is valid if it doesn't contain loops
        isvalid = (len(p)==len(set(p)))
        
    tcum = numpy.cumsum(td)+obstimes[0]
    nseg = len(p)-1
    
    # Generate measurements
    z = [[] for k in range(npts)]
    
    pos0 = G.node[p[0]]["pos"]
    pos1 = G.node[p[1]]["pos"]
    ell = numpy.linalg.norm(pos1-pos0)
    pos = pos0+(pos1-pos0)*s0/ell
    z[0] = pos+sig*numpy.random.randn(2)
    dprev = s0/ell
    jprev = 0
    for k in range(npts-2):
        j = (obstimes[k+1]>tcum).sum()
        if j<nseg:
            eweight = G[p[j]][p[j+1]]['weight']
            if (j != jprev):
                df = (obstimes[k+1]-tcum[j-1])*v[j]/eweight
            else:
                df = dprev+(obstimes[k+1]-obstimes[k])*v[j]/eweight
            pos0 = G.node[p[j]]["pos"]
            pos1 = G.node[p[j+1]]["pos"]
            pos = pos0+(pos1-pos0)*df
            z[k+1] = pos+sig*numpy.random.randn(2)
            dprev = df
            jprev = j
    
    pos0 = G.node[p[nseg-1]]["pos"]
    pos1 = G.node[p[nseg]]["pos"]
    ell = numpy.linalg.norm(pos1-pos0)
    pos = pos0+(pos1-pos0)*se/ell
    z[npts-1] = pos+sig*numpy.random.randn(2)
    
    return p, se, z
        
def calc_path_rmse(G,s0,p,se,phat,sehat):

# rmse = calc_path_rmse(G,s0,p,se,phat,sehat)
# Calculates the (spatial) rmse along true and estimated paths
# Inputs:
# G- road network
# s0- starting position (float)
# p- true path (list of integers)
# se- ending position (float)
# phat- estimated path (list of integers)
# sehat- estimated ending position (float)

    # Number of points along path
    npts = 101
    
    nseg = len(p)-1
    ell = [0 for j in range(nseg)]
    for j in range(nseg-1):
        ell[j] = mynorm(G.node[p[j+1]]["pos"]-G.node[p[j]]["pos"])
    ell[nseg-1] = se
    gsize = (sum(ell)-s0)/(npts-1)
    
    edge_idxs = 0
    ps = p[0]
    ss = s0
    z = numpy.zeros((2,npts))
    pos0 = G.node[ps]["pos"]
    z[:,0] = pos0
    # Find position at each point along true path
    for k in range(npts-1):
        s = ss+gsize
        while s>(1e-5+ell[edge_idxs]):
            s -= ell[edge_idxs]
            edge_idxs += 1
        ps = p[edge_idxs]
        pos0 = G.node[ps]["pos"]
        pos1 = G.node[p[edge_idxs+1]]["pos"]
        z[:,k+1] = pos0+s*(pos1-pos0)/ell[edge_idxs]
        ss = s
        
    nseg = len(phat)-1
    ell = [0 for j in range(nseg)]
    for j in range(nseg-1):
        ell[j] = mynorm(G.node[phat[j+1]]["pos"]-G.node[phat[j]]["pos"])
    ell[nseg-1] = sehat
    gsize = (sum(ell)-s0)/(npts-1)
    
    edge_idxs = 0
    ps = phat[0]
    ss = s0
    zhat = numpy.zeros((2,npts))
    pos0 = G.node[ps]["pos"]
    zhat[:,0] = pos0
    # Find position at each point along estimated path
    for k in range(npts-1):
        s = ss+gsize
        while s>(1e-5+ell[edge_idxs]):
            s -= ell[edge_idxs]
            edge_idxs += 1
        ps = phat[edge_idxs]
        pos0 = G.node[ps]["pos"]
        pos1 = G.node[phat[edge_idxs+1]]["pos"]
        zhat[:,k+1] = pos0+s*(pos1-pos0)/ell[edge_idxs]
        ss = s
        
    rmse = 0.
    for k in range(npts):
        epz = z[:,k]-zhat[:,k]
        rmse += mynorm(epz)/npts
    
    return rmse

def drawmultinom(wts):

# idx = drawmultinom(wts)
# Draw accoridng to given probabilities
# Input:
# wts- sampling probabilities (list of floats)
# Output:
# idx- sampled integer

    wcdf = numpy.cumsum(wts)
    idx = 0
    u = random.uniform(0,1)
    while (wcdf[idx]<u):
        idx += 1
        
    return idx

def tnorm(mu,kap,a,b):
    
# z = tnorm(mu,kap,a,b)
# Draw from a truncated normal in the interval (a,b)

    sqkap = math.sqrt(kap)
    lb = normcdf((a-mu)/sqkap,0,1)
    ub = normcdf((b-mu)/sqkap,0,1)
    u = lb+(ub-lb)*random.uniform(0,1)
    z = mu+sqkap*norm.ppf(u)
    
    return z

def get_rank_quarts(rhist):
    
# q = get_rank_quarts(rhist)
# Find the quartiles for a given histogram

    n = len(rhist)
    P = [0.25,0.50,0.75]
    q = [0,0,0]
    for j in range(len(P)):
        q[j] = n - sum(numpy.cumsum(rhist)>=P[j])+1
    
    return q
    
def mydot(a,b):

# d = mydot(a,b)
# Dot product

    n = len(a)
    d = 0.
    for i in range((n)):
        d = d+a[i]*b[i]
        
    return d
    
def mynorm(a):

# n = mynorm(a)
# Euclidean norm

    n = math.sqrt(mydot(a,a))
    
    return n

