import case1 
import psycopg2
import math
import testdb
import time
import numpy
import networkx as nx
import matplotlib.pyplot as plt

import csv
import sys

def reasonable_distance(start,end):
	if math.sqrt((start[0]-end[0])**2+(start[1]-end[1])**2)> 500:
		return False
	else:
		return True

def get_edges(route,cur):
	edgenum=len(route)-1
	edge=[]
	for i in range(edgenum):
		start=route[i]
		end = route[i+1]
		cur.execute('select getedgeid(%s,%s);' % (start,end))
		edge.append(str(cur.fetchone()[0]))
	return edge

def get_allmeasurements(measurements,vehicleid,cur):
	allmeasurements=[]
	cur.execute('select * from gettrajectory(%s, %s, %s);' % (vehicleid, measurements[0], measurements[-1]))
	for record in cur:
		allmeasurements.append(record)
	return allmeasurements

def RSME(edge, allmeasurements, cur):
	distance=[]
	for measurement in allmeasurements:
		cur.execute('select st_distance(geom_way,geom,true) as dis from hh_2po_4pgr a, ecourier_oneday b where b.id='+str(measurement[0])+' and a.id in ('+','.join(edge)+') order by dis limit 1;')
		distance.append(cur.fetchone()[0])
	total = 0.0
	for d in distance:
		total+=d*d
	return math.sqrt(total/len(distance))

def temporal_RSME(edge,allmeasurements,startvertex,cur):
	starttime = allmeasurements[0][1]
	total_time = (allmeasurements[-1][1]-starttime).total_seconds()
	distance=[]
	for measurement in allmeasurements:
		fraction = (measurement[1]-starttime).total_seconds()/total_time
		edgestr = "'{"+','.join(edge)+"}'"
		tempstr = 'select st_distance(ST_LineInterpolatePoint(mergeroute('+edgestr+' , '+str(startvertex)+'),'+str(fraction)+'),geom,true) as dis from ecourier_oneday b where b.id='+str(measurement[0])+';'
		print tempstr
		cur.execute(tempstr)
		distance.append(cur.fetchone()[0])
	total = 0.0
	for d in distance:
		total+=d*d
	return math.sqrt(total/len(distance))

def nearest_vertex(G,measurement):
	minimum_distance = 1000000
	for node in nx.nodes_iter(G):
		pos = G.node[node]['pos']
		distance = math.sqrt((measurement[0]-pos[0])**2+(measurement[1]-pos[1])**2)
		if distance< minimum_distance:
			minimum_distance=distance
			result = node
	return result

def case1_vehicle(vehicleid):
	boundingboxratio=3.0
	sig = 5
	# Algorithm parameters
	n = 100  # number of samples
	K = 100  # minimum number of candidate paths
	vmean = 10. # mean velocity
	vvar = 4. # velocity variance

	conn = psycopg2.connect("dbname=gis user=postgres password=peach")
	cur = conn.cursor()
	pointsqlx ="ST_X(ST_transform(ST_setsrid(ST_makepoint(longitude,latitude),4326),32630))"
	pointsqly ="ST_Y(ST_transform(ST_setsrid(ST_makepoint(longitude,latitude),4326),32630))"
	cur.execute("select id,"+pointsqlx+","+pointsqly+", timestamp from ecourier_oneday a where vehicleid=%s and exists (select * from hh_2po_4pgr where st_dwithin(source_geom,a.geom,0.00014) or st_dwithin(target_geom,a.geom,0.00014)) order by timestamp",(vehicleid,))
	rows = cur.fetchall()
	mesnum = len(rows)
	results = []
	for i in range(mesnum-1):
		pid,x,y,timex=rows[i]
		startid = pid
		startpoint = numpy.array([x,y])
		starttime = timex
		pid,x,y,timex=rows[i+1]
		endid = pid
		endpoint = numpy.array([x,y])
		endtime = timex
		if 20<(endtime-starttime).total_seconds()<100:
			IDlist = [startid,endid]
			z=[startpoint,endpoint]
			t=[0.0,(endtime-starttime).total_seconds()]

			minx = startpoint[0] if startpoint[0]<endpoint[0] else endpoint[0]
			maxx = startpoint[0] if startpoint[0]>endpoint[0] else endpoint[0]
			miny = startpoint[1] if startpoint[1]<endpoint[1] else endpoint[1]
			maxy = startpoint[1] if startpoint[1]>endpoint[1] else endpoint[1]

			width=maxx-minx
			height=maxy-miny
			minx=minx-width*boundingboxratio/2.0
			maxx=maxx+width*boundingboxratio/2.0
			miny=miny-height*boundingboxratio/2.0
			maxy=maxy+height*boundingboxratio/2.0
			G=testdb.london_roadmap((minx,miny,maxx,maxy))
			if G.nodes():
				startvertex = nearest_vertex(G,startpoint)
				endvertex = nearest_vertex(G,endpoint)

			print z
			print t
			if len(G.nodes())>10 and nx.has_path(G,startvertex,endvertex) and startvertex!=endvertex and reasonable_distance(z[0],z[1]):

				# position={}
				# for node in G:
				#   position[node]=G.node[node]['pos']
				# nx.draw(G,pos=position,node_size=50)
				# for point in z:
				#     plt.plot(point[0],point[1],'kx',markersize=10.0,markeredgewidth=2)
				# plt.show()

				# try:
				t0 = time.time()
				[pp,sep,wp] = case1.calcpostprobs(G,z,t,n,K,vmean,vvar,sig*sig)
				et = time.time()-t0
				print et, "secs, done algorithm"
				possibleroutes = len(pp)
				result=[]
				for j in range(possibleroutes):
					edge = get_edges(pp[j],cur)
					allmeasurements = get_allmeasurements(IDlist,vehicleid,cur)
					result.append({'route':pp[j],'possibility':wp[j],'rmse':RSME(edge,allmeasurements,cur),'rmse_temporal':temporal_RSME(edge,allmeasurements,pp[j][0],cur),'period':t[-1]})
				results.append(result)
				# except:
				# 	pass

	return results



results=case1_vehicle(16)
f = open('temp'+str(len(results))+'test.csv','wb')
try:
	writer = csv.writer(f)
	writer.writerow(('TestID','Interval','routeID','possibility','RSME','temporal_RSME'))
	tests = len(results)
	for i in range(tests):
		sorted_records=sorted(results[i], key=lambda x: x['possibility'], reverse=True)
		routes = len(sorted_records)
		for j in range(routes):
			writer.writerow((i+1,sorted_records[j]['period'],j+1,sorted_records[j]['possibility'],sorted_records[j]['rmse'],sorted_records[j]['rmse_temporal']))
finally:
	f.close()
print 'complete'

# pp, sep, wp, z=real_eg1()
# for route in pp:
# 	print route
# 	print RSME(route,z)
# print sep
# print wp