import case1 
import psycopg2
import math
import testdb
import time
import numpy
import networkx as nx
import matplotlib.pyplot as plt
from sta import Test
import pickle
import cProfile

import csv
import sys

# def reasonable_distance(start,end):
# 	if math.sqrt((start[0]-end[0])**2+(start[1]-end[1])**2)> 500:
# 		return False
# 	else:
# 		return True
class MeasureList:
	def __init__(self, records,interval):
		self.records=records
		self.length=len(records)
		self.start=0
		self.interval=interval
	def reasonable_speed(self,startrecord,endrecord):
		dist=math.sqrt((startrecord[2]-endrecord[2])**2+(startrecord[3]-endrecord[3])**2)
		duration = (endrecord[1]-startrecord[1]).total_seconds()
		if 5<dist/duration<40:
			return True
		else:
			return False
	def duration(self,startrecord,endrecord):
		low=self.interval-5
		high=self.interval+5
		if (endrecord[1]-startrecord[1]).total_seconds()<low:
			return -1
		elif (endrecord[1]-startrecord[1]).total_seconds()>high:
			return 1
		else:
			return 0
	def next(self):
		start=self.start
		end=start+1
		while end<self.length:
			startrecord=self.records[start]
			endrecord=self.records[end]
			t = self.duration(startrecord,endrecord)
			if t==-1:
				end+=1
			elif t==1:
				start+=1
				end=start+1
			else:
				if self.reasonable_speed(startrecord,endrecord):
					locations=[numpy.array([startrecord[2],startrecord[3]]),numpy.array([endrecord[2],endrecord[3]])]
					times=[startrecord[1],endrecord[1]]
					ids = [startrecord[0],endrecord[0]]
					self.start=start+1
					return locations,times,ids
				else:
					start+=1
					end=start+1
		self.start = self.length-1
		return [],[],[]



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
	cur.execute("select * from gettrajectory(%s, '%s', %s, %s);" % (vehicleid[0],vehicleid[1], measurements[0], measurements[-1]))
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
		# print tempstr
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

def surrounding_area(modelpath,modelstartvertex,truepath,startvertex,cur):
	modelpathstr = "'{"+','.join(modelpath)+"}'"+' , '+str(modelstartvertex)
	truepathstr = "'{"+','.join(truepath)+"}'"+' , '+str(startvertex)
	components = 'ARRAY[m,t,st_makeline(st_startpoint(m),st_startpoint(t)),st_makeline(st_endpoint(m),st_endpoint(t))]'
	query = "select st_area(st_makepolygon(st_linemerge(st_collect("+components+"))),true) from mergeroute("+modelpathstr+") m, mergeroute("+truepathstr+") t ;"
	# print query
	cur.execute(query)
	return cur.fetchone()[0]

def share_distance(modelpath,modelstartvertex,truepath,startvertex,cur):
	modelpathstr = "'{"+','.join(modelpath)+"}'"+' , '+str(modelstartvertex)
	truepathstr = "'{"+','.join(truepath)+"}'"+' , '+str(startvertex)
	query = "select st_length(st_intersection(m,t))/st_length(t) from mergeroute("+modelpathstr+") m, mergeroute("+truepathstr+") t ;"
	cur.execute(query)
	return cur.fetchone()[0]

def real_path(G,allmeasurements,startvertex,endvertex,cur):
	pointonedge=[]
	if len(allmeasurements)>2:
		i = 570000
		middlepoint = [measurement[0] for measurement in allmeasurements[1:-1]]
		query = "select b.source as source, b.target as target, st_linelocatepoint(b.geom_way,a.geom) as fraction, st_x(st_transform(st_lineinterpolatepoint(b.geom_way,st_linelocatepoint(b.geom_way,a.geom)),32630)),st_y(st_transform(st_lineinterpolatepoint(b.geom_way,st_linelocatepoint(b.geom_way,a.geom)),32630)) from ecourier_oneday a, hh_2po_4pgr b where a.id=%s and st_dwithin(a.geom,b.geom_way,0.00015) order by st_distance(a.geom,b.geom_way) limit 1;"
		for point in middlepoint:
			cur.execute(query % (point,))
			measurement_is_valid=cur.fetchone()
			if measurement_is_valid:
				source, target, fraction, x, y = measurement_is_valid
				if G.has_edge(source,target):
					weight= G[source][target]['weight']
					G.add_node(i,pos=numpy.array([x,y]))
					G.add_edge(source,i,weight=fraction*weight)
					G.add_edge(i,target,weight=(1-fraction)*weight)
					G.remove_edge(source,target)
					pointonedge.append(i)
					i+=1
	pointonedge= [startvertex]+pointonedge+[endvertex] 
	# shortest path on every two points and concatenate 
	realpath=[]
	for j in range(len(pointonedge)-1):
		realpath+=nx.dijkstra_path(G,pointonedge[j],pointonedge[j+1])
	# remove duplicate points
	simple_realpath=[]
	for point in realpath:
		if simple_realpath:
			if point not in pointonedge[1:-1]:
				if len(simple_realpath)>1 and point==simple_realpath[-2]:
					simple_realpath.pop()
				elif point!=simple_realpath[-1]:
					simple_realpath.append(point)
		else:
			simple_realpath.append(point)
	return simple_realpath

def cal_min_bbox(allmeasurements,cur):
	allcoordinates=[]
	pointsqlx ="ST_X(ST_transform(ST_setsrid(ST_makepoint(longitude,latitude),4326),32630))"
	pointsqly ="ST_Y(ST_transform(ST_setsrid(ST_makepoint(longitude,latitude),4326),32630))"
	for measurement in allmeasurements:
		cur.execute("select "+pointsqlx+","+pointsqly+" from ecourier_oneday where id="+str(measurement[0])+"; ")
		allcoordinates.append(cur.fetchone())
	minx = min([p[0] for p in allcoordinates])
	maxx = max([p[0] for p in allcoordinates])
	miny = min([p[1] for p in allcoordinates])
	maxy = max([p[1] for p in allcoordinates])
	return minx, maxx, miny, maxy

def case1_vehicle(vehicleid,interval):
	boundingboxratio=0.6
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
	cur.execute("select id, timestamp,"+pointsqlx+","+pointsqly+" from near_vertex_measurement where vehicleid=%s and type=%s order by timestamp",vehicleid)
	vehicle_records = MeasureList([measurement for measurement in cur],interval)
	zl,tl,idl=vehicle_records.next()
	# mesnum = len(rows)
	results = []
	while zl:
		starttime=tl[0]
		t=[(item-starttime).total_seconds() for item in tl]
		allmeasurements = get_allmeasurements(idl,vehicleid,cur)

		minx, maxx, miny, maxy = cal_min_bbox(allmeasurements,cur)
		width=maxx-minx
		height=maxy-miny
		minx=minx-width*boundingboxratio/2.0
		maxx=maxx+width*boundingboxratio/2.0
		miny=miny-height*boundingboxratio/2.0
		maxy=maxy+height*boundingboxratio/2.0

		G=testdb.london_roadmap((minx,miny,maxx,maxy))
		if G.nodes():
			startvertex = nearest_vertex(G,zl[0])
			endvertex = nearest_vertex(G,zl[-1])

		print zl
		print t
		# large graph take too long time to process
		if 5000>len(G.nodes())>10 and nx.has_path(G,startvertex,endvertex) and startvertex!=endvertex:
			try:
				t0 = time.time()
				[pp,sep,wp] = case1.calcpostprobs(G,zl,t,n,K,vmean,vvar,sig*sig)
				et = time.time()-t0
				print et, "secs, done algorithm"
				possibleroutes = len(pp)
				
				if allmeasurements and allmeasurements[0][0]==idl[0] and allmeasurements[-1][0]==idl[-1]:
					tpath = real_path(G,allmeasurements,startvertex,endvertex,cur)
					# print "TRUE PATH", tpath
					tedge = get_edges(tpath,cur)
					result=[]
					for j in range(possibleroutes):
						# print "modeled path:",pp[j]
						edge = get_edges(pp[j],cur)
						try:
							area = surrounding_area(edge,pp[j][0],tedge,startvertex,cur)
						except:
							conn.rollback()
							area = surrounding_area(edge,pp[j][-1],tedge,startvertex,cur)
						result.append({'route':pp[j],
							'possibility':wp[j],
							'RMSE':RSME(edge,allmeasurements,cur),
							'temporal_RMSE':temporal_RSME(edge,allmeasurements,pp[j][0],cur),
							'area':area,
							'common_length':share_distance(edge,pp[j][0],tedge,startvertex,cur),
							'edgeIDList':edge
							})
					newtest = Test(result)
					newtest.test_para = {'period':t[-1],'measurements':[x[0] for x in allmeasurements],"tpath":tpath,"tedge":tedge,"vehicleid":vehicleid}
					results.append(newtest)
					if len(results)>30:
						return results
					else:
						print "Success.************************", len(results)
				else:
					print "vertx edge table error",idl,allmeasurements
			except:
				conn.rollback()
				print "FAIL............................."
		zl,tl,idl=vehicle_records.next()
	return results

def profile_main():
	interval = 30
	conn = psycopg2.connect("dbname=gis user=postgres password=peach")
	cur = conn.cursor()
	cur.execute("select distinct vehicleid, type from near_vertex_measurement;")
	results=[]
	for vehicleid in cur:
		results+=case1_vehicle(vehicleid,interval)
		if len(results)>=100:
			break
	# results=case1_vehicle(8)
	f = open("f:/London/statistics3/"+str(interval)+'seconds'+str(len(results))+'test.p','wb')
	pickle.dump(results,f)
	# try:
	# 	writer = csv.writer(f)
	# 	writer.writerow(('TestID','Interval','routeID','possibility','RSME','temporal_RSME','area','common_length'))
	# 	tests = len(results)
	# 	for i in range(tests):
	# 		sorted_records=sorted(results[i], key=lambda x: x['possibility'], reverse=True)
	# 		routes = len(sorted_records)
	# 		for j in range(routes):
	# 			writer.writerow((i+1,sorted_records[j]['period'],j+1,sorted_records[j]['possibility'],sorted_records[j]['rmse'],sorted_records[j]['rmse_temporal'],sorted_records[j]['area'],sorted_records[j]['commonlength']))
	# finally:
	# 	f.close()
	print 'complete'

# cProfile.run('profile_main()')
if __name__=='__main__':
	profile_main()
