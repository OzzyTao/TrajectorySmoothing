import random
import numpy
import math
import networkx as nx
import matplotlib.pyplot as plt
from scipy.stats import norm
import sys
sys.path.append("F:/London/K_shortest_paths/ba6eccab9a83e54e95dd")
from k_shortest_paths import k_shortest_paths

def makeemptygrid(n,sdlen):

# G = makeemptygrid(n,sdlen)
# Constructs a grid network with no edges

    nnodes = n*n
    
    G = nx.Graph()
    # Add nodes
    for i in range(nnodes):
        G.add_node(i,pos=numpy.array([i%n,i/n]))
                
    return G


def makegridgraph(n,sdlen):

# G = makegridgraph(n,sdlen)
# Constructs a n x n grid network with length sdlen between adjacent nodes

    nnodes = n*n
    
    G = nx.Graph()
    # Add nodes
    for i in range(nnodes):
        G.add_node(i,pos=numpy.array([i%n,i/n]))
        
    # Add edges
    for i in range(nnodes):
        pi = G.node[i]["pos"]
        for j in numpy.arange(i+1,nnodes):
            pj = G.node[j]["pos"]
            if abs(pi-pj).sum()==1:
                G.add_edge(i,j,weight=sdlen)
        G.node[i]["pos"] = sdlen*G.node[i]["pos"]
        
    return G
    
def create_turn_tab(G):

# Tpen = create_turn_tab(G)
# Creates a turn penalty table
# Inputs:
# G- road network
# Outputs:
# Tpen- 2D array of penalties for transitioning from one edge to another
# (size is 2*nedges x 2*nedges, only the top half is populated).

    e = G.edges()
    nedges = len(e)
    Tpen = numpy.zeros((2*nedges,2*nedges)) # Need to allow for both directions
    
    
    for i in range(nedges):
        for j in range(i+1,nedges):
            if e[i][1]==e[j][0]: # Connected edges
                # Determine the nature of the junction
                nn = len(G.neighbors(e[j][0]))
                # Determine the angle of the turn
                pos0 = G.node[e[i][0]]["pos"]
                pos1= G.node[e[j][0]]["pos"]
                pos2 = G.node[e[j][1]]["pos"]
                ep1 = pos1-pos0
                ep2 = pos2-pos1
                th = myacos(mydot(ep2,ep1)/(mynorm(ep1)*mynorm(ep2)))
                if nn==1:
                    Tpen[i,j] = 0.
                elif nn==2:
                    Tpen[i,j] = 5*(abs(th)**4)
                else:
                    Tpen[i,j] = 10*(abs(th)**4)
        for j in range(nedges):
            if e[i][1]==e[j][1]: # Connected edges
                # Determine the nature of the junction
                nn = len(G.neighbors(e[j][1]))
                # Determine the angle of the turn
                pos0 = G.node[e[i][0]]["pos"]
                pos1= G.node[e[j][1]]["pos"]
                pos2 = G.node[e[j][0]]["pos"]
                ep1 = pos1-pos0
                ep2 = pos2-pos1
                th = myacos(mydot(ep2,ep1)/(mynorm(ep1)*mynorm(ep2)))
                if nn==1:
                    Tpen[i,j+nedges] = 0.
                elif nn==2:
                    Tpen[i,j+nedges] = 10*(abs(th)**4)
                else:
                    Tpen[i,j+nedges] = 20*(abs(th)**4)
    for i in range(nedges):
        for j in range(i+1,nedges):
            if e[i][0]==e[j][1]: # Connected edges
                # Determine the nature of the junction
                nn = len(G.neighbors(e[j][0]))
                # Determine the angle of the turn
                pos0 = G.node[e[i][1]]["pos"]
                pos1= G.node[e[j][1]]["pos"]
                pos2 = G.node[e[j][0]]["pos"]
                ep1 = pos1-pos0
                ep2 = pos2-pos1
                th = myacos(mydot(ep2,ep1)/(mynorm(ep1)*mynorm(ep2)))
                if nn==1:
                    Tpen[i+nedges,j+nedges] = 0.
                elif nn==2:
                    Tpen[i+nedges,j+nedges] = 10*(abs(th)**4)
                else:
                    Tpen[i+nedges,j+nedges] = 20*(abs(th)**4)
    
    return Tpen

    
def makedata(G,p,obstimes,s0,se,sig,vmax,a,b):

# [z,t] = makedata(G,p,obstimes,s0,se,sig,vmax,a,b)
# Generate position measurements
# Inputs:
# G- road network
# p- path (list of integers)
# obstimes- observation times (list of floats)
# s0- starting position on first edge (float)
# se- ending position on last edge (float)
# sig- measurement noise standard deviation (float)
# (vmax,a,b)- vehicle motion statistics (floats)
# Outputs:
# z- position measurements (list of list of floats)
# t- total duration (float)

    nseg = len(p)-1
    npts = len(obstimes)
    
    z = [[] for k in range(npts)]#numpy.zeros((2,npts))
    
    # Generate velocities and calculate durations along each segment
    v = numpy.zeros((nseg))
    td = numpy.zeros((nseg))
    pos0 = G.node[p[0]]["pos"]
    pos1 = G.node[p[1]]["pos"]
    ell = numpy.linalg.norm(pos1-pos0)
    v[0] = vmax*random.betavariate(a,b)
    td[0] = (ell-s0)/v[0]
    for j in range(1,nseg-1):
        pos0 = G.node[p[j]]["pos"]
        pos1 = G.node[p[j+1]]["pos"]
        ell = numpy.linalg.norm(pos1-pos0)
        v[j] = vmax*random.betavariate(a,b)
        td[j] = ell/v[j]
    j = nseg-1
    ell = se
    v[j] = vmax*random.betavariate(a,b)
    td[j] = ell/v[j]
    tcum = td.cumsum()+obstimes[0]
    t = td.sum()+obstimes[0]
    
    # Generate measurements
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
    
    return z, t

def makedata2(G,pb,pa,xM,tM,T,obstimes,s0,se,sig,vmax,a,b):

# z = makedata2(G,pb,pa,xM,tM,T,obstimes,s0,se,sig,vmax,a,b)
# Generates position measurements for travel over a specified time with a specified 
# meeting point
# Inputs:
# G- road network
# pb- path before meeting (list of integers)
# pa- path after meeting (list of integers)
# xM- meeting position (array of floats)
# tM- meeting time (float)
# T- duration of meeting (float)
# obstimes- observation times (list of floats)
# s0- starting position on first edge (float)
# se- ending position on last edge (float)
# sig- measurement noise standard deviation (float)
# (vmax,a,b)- vehicle motion statistics (floats)
# Outputs:
# z- position measurements (list of list of floats)

    nsegb = len(pb)-1
    nsega = len(pa)-1
    npts = len(obstimes)
    tol = 0.05
    
    nobsa = 0
    nobsb = 0
    nobsd = 0
    for i in range(npts):
        if obstimes[i]<=tM:
            nobsb = nobsb+1
        elif obstimes[i]<=(tM+T):
            nobsd = nobsd = nobsd+1
        else:
            nobsa = nobsa+1

    z = [[] for k in range(npts)]
    
    
    # Before meeting
    isok = 0
    vb = numpy.zeros((nsegb))
    td = numpy.zeros((nsegb))
    while not isok: # Specified meeting time not yet satisfied
        # Generate velocities and calculate durations for each segment
        pos0 = G.node[pb[0]]["pos"]
        pos1 = G.node[pb[1]]["pos"]
        ell = numpy.linalg.norm(pos1-pos0)
        vb[0] = vmax*random.betavariate(a,b)
        td[0] = (ell-s0)/vb[0]
        for j in range(1,nsegb-1):
            pos0 = G.node[pb[j]]["pos"]
            pos1 = G.node[pb[j+1]]["pos"]
            ell = numpy.linalg.norm(pos1-pos0)
            vb[j] = vmax*random.betavariate(a,b)
            td[j] = ell/vb[j]
        j = nsegb-1
        pos0 = G.node[pb[j]]["pos"]
        sM = mynorm(xM-pos0)
        vb[j] = vmax*random.betavariate(a,b)
        td[j] = sM/vb[j]
        tcum = td.cumsum()+obstimes[0]
        tb = td.sum()+obstimes[0]
        # Accept path if its close to specified time
        if abs(tb-tM)<tol:
            isok = 1
    
    # Generate data points
    dprev = s0
    jprev = 0
    for k in range(nobsb):
        j = 0
        while tcum[j]<obstimes[k]:
            j += 1
        eweight = G[pb[j]][pb[j+1]]['weight']
        if j != jprev:
            df = (obstimes[k]-tcum[j-1])*vb[j]/eweight
        else:
            df = (dprev+obstimes[k]*vb[j])/eweight
        pos0 = G.node[pb[j]]["pos"]
        pos1 = G.node[pb[j+1]]["pos"]
        pos = pos0+(pos1-pos0)*df
        z[k] = pos+sig*numpy.random.randn(2)
        dprev = df
        jprev = j
    
    # During meeting
    for k in range(nobsd):
        z[nobsb+k] = xM+sig*numpy.random.randn(2)
        
    # After meeting
    isok = 0
    va = numpy.zeros((nsega))
    td = numpy.zeros((nsega))
    dM = numpy.linalg.norm(xM-G.node[pa[0]]["pos"])
    while not isok: # Specified ending time not yet satisfied
        j = 0
        # Generate velcoities and calculate durations along each segment
        eweight = G[pa[0]][pa[1]]['weight'] 
        va[j] = vmax*random.betavariate(a,b)
        td[j] = (eweight-dM)/va[j]
        for j in range(nsega-2):
            eweight = G[pa[j+1]][pa[j+2]]['weight'] 
            va[j+1] = vmax*random.betavariate(a,b)
            td[j+1] = eweight/va[j+1]
        va[nsega-1] = vmax*random.betavariate(a,b)
        td[nsega-1] = se/va[nsega-1]
        tcum = td.cumsum()
        ta = td.sum()
        # Accept path if its close to specified time
        if abs(tb+T+ta-obstimes[-1])<tol:
            isok = 1

    # Generate data points
    dprev = dM
    jprev = 0
    for k in range(nobsa-1):
        j = 0
        while (tb+T+tcum[j])<obstimes[k+nobsb+nobsd]:
            j += 1
        eweight = G[pa[j]][pa[j+1]]['weight']
        if j != jprev:
            df = (obstimes[k+nobsb+nobsd]-(tb+T+tcum[j-1]))*va[j]/eweight
        else:
            df = (dprev+(obstimes[k+nobsb+nobsd]-(tb+T+tcum[j-1]))*va[j])/eweight
        pos0 = G.node[pa[j]]["pos"]
        pos1 = G.node[pa[j+1]]["pos"]
        pos = pos0+(pos1-pos0)*df
        z[k+nobsd+nobsb] = pos+sig*numpy.random.randn(2)
        dprev = df
        jprev = j
    pos0 = G.node[pa[nsega-1]]["pos"]
    pos1 = G.node[pa[nsega]]["pos"]
    ell = numpy.linalg.norm(pos1-pos0)
    pos = pos0+(pos1-pos0)*se/ell
    z[nobsa-1+nobsd+nobsb] = pos+sig*numpy.random.randn(2)
        
    return z

def closest_point(p,p1,p2):
    
# [d,s] = closest_point(p,p1,p2)
# Finds the closest point from p to the line joining p1 and p2
# Outputs:
# d- distance from p to line (p1,p2)
# s- distance of closest point from p1
    x = float(p[0])
    y = float(p[1])
    x1 = float(p1[0])
    y1 = float(p1[1])
    x2 = float(p2[0])
    y2 = float(p2[1])
    
    
    d1 = math.sqrt((x1-x)*(x1-x)+(y1-y)*(y1-y))
    d2 = math.sqrt((x2-x)*(x2-x)+(y2-y)*(y2-y))
        
    ell = (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)
    u = (x2-x1)*(x-x1)+(y2-y1)*(y-y1)
    if u<0.:
        xi = x1
        yi = y1
        d = math.sqrt((x1-x)*(x1-x)+(y1-y)*(y1-y))
    elif u>ell:
        xi = x2
        yi = y2
        d = math.sqrt((x2-x)*(x2-x)+(y2-y)*(y2-y))
    elif True:
        b = u/ell
        xi = x1+b*(x2-x1)
        yi = y1+b*(y2-y1)
        d =  math.sqrt((xi-x)*(xi-x)+(yi-y)*(yi-y))
    
    s = math.sqrt((xi-x1)*(xi-x1)+(yi-y1)*(yi-y1))

    return d, s

def find_start_paths(G,z,ups):
    
# [eidx,epos] = find_start_paths(G,z,ups)
# Find candidate edges on which an object could lie
# Inputs:
# G- road network
# z- measurement (list of floats)
# ups- measurement noise variance (float)
# Outputs:
# epos- coordinates of nodes of candidate edges (list of 2D arrays) [[[x11,y11],[x12,y12]],[[x21,y21],[x22,y22]],...]
# eidx- indices of nodes of candidate edges (list of list of integers)

    [node_idx,near_edges] = findnodes(G,z,ups)
    nedges = len(near_edges)
    
    path_pos = []
    path_idx = []
    for i in range(nedges):
        j1 = near_edges[i][0]
        j2 = near_edges[i][1]
        pos1 = G.node[j1]["pos"]
        pos2 = G.node[j2]["pos"]
        path_pos.append(numpy.array([[pos1[0],pos2[0]],[pos1[1],pos2[1]]]))
        path_idx.append([j1,j2])
        
    return path_pos, path_idx
    
def find_paths(G,Tp_tab,z,p0,ups,K):

# [pathidx,pathpos,turnpen] = find_paths(G,Tp_tab,z,p0,ups,K):
# Construct candidate path set
# Inputs:
# G- road network
# Tp_tab- table of turn penalties (2D array of floats)
# z- measurement (list of floats)
# p0- current path (list of integers)
# ups- measurement noise variance (float)
# K- minimum number of candidate paths (integer)
# Outputs:
# path_pos- coordinates of nodes of candidate paths (list of arrays of floats) [[[x11,y11],[x12,y12]],[[x21,y21],[x22,y22]],...]
# path_idx- indices of nodes of candidate paths (list of list of integers)
# turnpen: turn penalties for each candidate path (list of floats)

    [node_idx,near_edges] = findnodes(G,z,ups)
    nnodes = len(node_idx)
    
    e = G.edges()
    nedges = len(e)
    for j in range(nedges):
        eadd = (e[j][1],e[j][0])
        e.append(eadd)
        
    path_pos = []
    path_idx = []
    turnpen = []
    end_edge = list(numpy.zeros((2)))
    for i in range(nnodes):
        j = node_idx[i]
        numcandp = 0
        d = 0
        if p0[-1]==j:
            candp0 = [[p0[-1]]]
        else:
            candp0 = []
        while numcandp<K:
            candp = candp0+list(nx.all_simple_paths(G,p0[-1],j,d))
            numcandp = len(candp)
            d = d+1
        # Only retain paths which include feasible end edges
        for k in range(numcandp):
            if len(candp[k])==1:
                end_edge[0] = p0[-2]
            else:
                end_edge[0] = candp[k][-2]
            end_edge[1] = candp[k][-1]
            if end_edge in near_edges:
                pp = candp[k]
                pp.insert(0,p0[-2])
                path_idx.append(pp)
                nnodes = len(pp)
                newpath = numpy.zeros((2,nnodes))
                for a in range(nnodes):
                    posa = G.node[pp[a]]["pos"]
                    newpath[0,a] = posa[0]
                    newpath[1,a] = posa[1]
                path_pos.append(newpath)
                tpen = 0.
                for a in range((nnodes-2)):
                    e1 = (pp[a],pp[a+1])
                    j1 = e.index(e1)
                    e2 = (pp[a+1],pp[a+2])
                    j2 = e.index(e2)
                    if (j1>j2):
                        tpen = tpen+Tp_tab[j2,j1]
                    else:
                        tpen = tpen+Tp_tab[j1,j2]
                turnpen.append(tpen)
            end_edge[0] = candp[k][-1]
            if len(candp[k])==1:
                end_edge[1] = p0[-2]
            else:
                end_edge[1] = candp[k][-2]
            if end_edge in near_edges:
                pp = candp[k]
                pp.insert(0,p0[-2])
                path_idx.append(pp)
                nnodes = len(pp)
                newpath = numpy.zeros((2,nnodes))
                for a in range(nnodes):
                    posa = G.node[pp[a]]["pos"]
                    newpath[0,a] = posa[0]
                    newpath[1,a] = posa[1]
                path_pos.append(newpath)
                tpen = 0.
                for a in range((nnodes-2)):
                    e1 = (pp[a],pp[a+1])
                    j1 = e.index(e1)
                    e2 = (pp[a+1],pp[a+2])
                    j2 = e.index(e2)
                    if (j1>j2):
                        tpen = tpen+Tp_tab[j2,j1]
                    else:
                        tpen = tpen+Tp_tab[j1,j2]
                turnpen.append(tpen)
                
    return path_pos, path_idx, turnpen
    


def find_paths_by_IDs(G,Tp_tab,end_node_idx,end_near_edges,start_node_idx,K,tpath):

# [pathidx,pathpos, turnpen] = find_paths_by_coordinates(G,Tp_tab,end_node_idx,end_near_edges,start_node_idx,K):
# Finds the candidate edges for a group of starting nodes
# Inputs:
# G- road network
# Tp_tab- table of turn penalties
# end_node_idx- indices of ending nodes (list of integers)
# end_near_edges- ending edges (list of list of integers)
# start_node_idx- starting nodes (list of integers)
# K- minimum number of candidate paths (integer)
# Outputs:
# path_pos- coordinates of nodes of candidate paths (list of arrays of floats) 
# path_idx- indices of nodes of candidate paths (list of list of integers)
# turnpen: turn penalties for each candidate path (list of floats)
    
        
    e = G.edges()
    nedges = len(e)
    for j in range(nedges):
        eadd = (e[j][1],e[j][0])
        e.append(eadd)
        
    # Insert virtual edges into the graph
    graph = insertVirtualEdges(G, start_node_idx, end_node_idx)
    
    # Compute the k candidate paths and their turning cost
    path_pos = []
    path_idx = []
    turnpen = []
    end_edge = list(numpy.zeros((2)))
    
    [candp, numcandp] = computeKCandidatepaths(graph, K, tpath)
    
    G = remove_virtual_edges(graph)
    
    # Remove single node paths
    candp = [p for p in candp if len(p)>1]
    numcandp = len(candp)

    # Insert true path into candidate list
    # if true_path and true_path not in candp:
    #     candp = [true_path]+candp
    #     numcandp +=1


    # Retain paths with the correct ending edges
    for k in range(numcandp):
        end_edge[0] = candp[k][-2]
        end_edge[1] = candp[k][-1]
        if end_edge in end_near_edges:
            pp = candp[k]
            path_idx.append(pp)
            nnodes = len(pp)
            newpath = numpy.zeros((2,nnodes))
            for a in range(nnodes):
                posa = G.node[pp[a]]["pos"]
                newpath[0,a] = posa[0]
                newpath[1,a] = posa[1]
            path_pos.append(newpath)
            tpen = 0.
            for a in range((nnodes-2)):
                e1 = (pp[a],pp[a+1])
                j1 = e.index(e1)
                e2 = (pp[a+1],pp[a+2])
                j2 = e.index(e2)
                if (j1>j2):
                    tpen = tpen+Tp_tab[j2,j1]
                else:
                    tpen = tpen+Tp_tab[j1,j2]
            turnpen.append(tpen)
        end_edge[0] = candp[k][-1]
        end_edge[1] = candp[k][-2]
        if end_edge in end_near_edges:
            pp = candp[k]
            path_idx.append(pp)
            nnodes = len(pp)
            newpath = numpy.zeros((2,nnodes))
            for a in range(nnodes):
                posa = G.node[pp[a]]["pos"]
                newpath[0,a] = posa[0]
                newpath[1,a] = posa[1]
            path_pos.append(newpath)
            tpen = 0.
            for a in range((nnodes-2)):
                e1 = (pp[a],pp[a+1])
                j1 = e.index(e1)
                e2 = (pp[a+1],pp[a+2])
                j2 = e.index(e2)
                if (j1>j2):
                    tpen = tpen+Tp_tab[j2,j1]
                else:
                    tpen = tpen+Tp_tab[j1,j2]
            turnpen.append(tpen)
            
    return path_pos, path_idx, turnpen
    


def computeKCandidatepaths(graph, K,tpath):

# [candp,numcandp] = computeKCandidatepaths(graph, K)
# Computes multi-source candidate paths
# graph- network
# K- minimum number of paths (integer)
# Outputs:
# candp- paths (list of list of integers)
# numcandp- number of paths (integer)
    graphcopy=graph.copy()
    idStart = 1000001
    idStop = 1000002
    numcandp = 0
    d = 0
    #compute k candidate paths given the maximum number of segment in each path  

    # if len(list(nx.all_simple_paths(graph,idStart,idStop)))<K:
    #     raise Exception("To simple graph")

    # while numcandp<K and d<100:
    #     candp = list(nx.all_simple_paths(graph, idStart, idStop, d))
    #     numcandp = len(candp)
    #     d += 1
    candp = k_shortest_paths(graphcopy,idStart,idStop,k=K,weight=None)[-1]
    #remove the virtual start node and end node
    for pth in candp:
        pth.pop()
        pth.pop(0)
    if tpath and tpath not in candp:
        candp.append(tpath)   
    return candp, len(candp)
    
    
def insertVirtualEdges(G, start_node_idx, end_node_idx):

# G = insertVirtualEdges(G, start_node_idx, end_node_idx):
# Inserts virtual edges into the graph

    idStart = 1000001
    idStop =1000002
    
    for nid in start_node_idx:
        G.add_edge(idStart, nid, weight = 0.01)           
        
    for nid in end_node_idx:
        G.add_edge(nid, idStop, weight = 0.01)    
    
    return G
    
def remove_virtual_edges(G):

# G = remove_virtual_edges(G)
# Removes virtual edges from the graph

    idStart = 1000001
    idStop =1000002
    G.remove_node(idStart)
    G.remove_node(idStop)
    
    
    return G
    
def remove_dups(s):

# sout = remove_dups(s)
# Removes duplicate elements from a list s

    n = len(s)
    sout = []
    for i in range(n):
        if s[i] not in sout:
            sout.append(s[i])
    
    return sout
    
def findnodes(G,z,ups):

# [node_idx,e] = findnodes(G,z,ups)
# Finds nodes and edges near a position
# Input:
# G- graph
# z- position (list of floats)
# ups- measurement noise variance (float)
# Outputs:
# node_idx- node indices (list of integers)
# e- edge indices  (list of list of integers)

    Th = 18.42 # 99.9% prob

    ed = G.edges()
    nedges = len(ed)
    node_idx = []
    e = []
    for i in range(nedges):
        j1 = ed[i][0]
        j2 = ed[i][1]
        pos1 = G.node[j1]["pos"]
        pos2 = G.node[j2]["pos"]
        [d,scl] = closest_point(z,pos1,pos2)
        dnorm = d*d/ups
        if (dnorm<Th):
            e.append(list(ed[i]))
            if j1 not in node_idx:
                node_idx.append(j1)
            if j2 not in node_idx:
                node_idx.append(j2)
    
    return node_idx, e
    
    
def getpos(G,p,d,sidelen):

# pos = getpos(G,p,d,sidelen)
# Returns the position of a distance along a path (a grid network is assumed)
# Inputs:
# G- road network
# p- path (list of integers)
# d- distance (float)
# sidelen- distance between nodes (float)
# Output:
# pos- position (array of floats)

    j = int(numpy.floor(d/sidelen))
    if j==len(p):
        j = len(p)-1
        df = 1
    else:
        df = (d-j*sidelen)/sidelen
 
    pos0 = G.node[p[j]]["pos"]
    pos1 = G.node[p[j+1]]["pos"]
    pos = pos0+(pos1-pos0)*df
    
    return pos
    
def getseg(G,p,v,tcum,t,sidelen):
    
# [j,seg,pos] = getseg(G,p,v,tcum,t,sidelen)
# Finds the segment occupied by the object at a particular time
# Inputs:
# G- road network
# p- path (list of integers)
# v- velocity along each edge (list of floats)
# tcum- accumulated travel time (list of floats)
# t- time of interest (float)
# sidelen- distance between nodes (float)
# Outputs:
# j- edge number (integer)
# seg- edge (list of two integers)
# pos- 2D position (list of two floats)

    j = 0
    while tcum[j]<t:
        j += 1
 
    seg = [p[j],p[j+1]]
    pos0 = G.node[p[j]]["pos"]
    pos1 = G.node[p[j+1]]["pos"]
    if j>0:
        dseg = v[j]*(t-tcum[j-1])
    else:
        dseg = v[j]*t
    pos = pos0+(pos1-pos0)*dseg/sidelen
    
    return j, seg, pos
   

def plotpath(G,p,lines,offset):
    
# plotpath(G,p,lines,offset)
# Plots the road network and paths
# Inputs:
# G- road network
# p- path (list of integers)
# lines- line style
# offset- offset of path from network (float)

    #plt.clf()
    nnodes = len(G.node)
    for j in range(nnodes):
        pos0 = G.node[j]["pos"]
        plt.plot(pos0[0],pos0[1],'k.',markersize=7.)
        
    elist = G.edges()
    nedges = len(elist)
    for j in range(nedges):
        pos0 = G.node[elist[j][0]]["pos"]
        pos1 = G.node[elist[j][1]]["pos"]
        plt.plot([pos0[0],pos1[0]],[pos0[1],pos1[1]],'k-',linewidth=0.75)
        
    cols = [[0,0,1],[1,0,0],[0,0.7,0],[0.5,0,0.5],[0.8,0.5,0],[0,0.5,0.8]]

    npaths = len(p)
    for i in range(npaths):
        nsegs = len(p[i])-1
        j = 0
        pos0 = G.node[p[i][j]]["pos"]
        plt.plot(pos0[0]+offset,pos0[1]+offset,'x', color=cols[i],mew=2.,markersize=10.)
        for j in range((nsegs)):
            pos0 = G.node[p[i][j]]["pos"]
            pos1 = G.node[p[i][j+1]]["pos"]
            xpos = numpy.array((pos0[0]+offset,pos1[0]+offset))
            ypos = numpy.array((pos0[1]+offset,pos1[1]+offset))
            plt.plot(xpos,ypos,color=cols[i],linewidth=2.5,linestyle=lines)
    
    # Axis
    min_xpos = 10000
    max_xpos = -10000
    min_ypos = 10000
    max_ypos = -10000
    for k in range(len(G.nodes())):
        pos = G.node[k]["pos"]
        if pos[0]>max_xpos:
            max_xpos = pos[0]
        if pos[0]<min_xpos:
            min_xpos = pos[0]
        if pos[1]>max_ypos:
            max_ypos = pos[1]
        if pos[1]<min_ypos:
            min_ypos = pos[1]
            
    plt.axis([min_xpos-50,max_xpos+50,min_ypos-50,max_ypos+50])
    plt.xlabel('x-position (m)',fontsize=16)
    plt.ylabel('y-position (m)',fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=14)
    
def mydot(a,b):

    n = len(a)
    d = 0.
    for i in range((n)):
        d = d+a[i]*b[i]
        
    return d

def mynorm(a):

    n = math.sqrt(mydot(a,a))
    
    return n
    
def myacos(value):
    value=-1.0 if value<-1.0 else 1.0 if value>1.0 else value
    return math.acos(value)      
    
            