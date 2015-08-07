import networkx as nx
import psycopg2
import numpy
import time
from matplotlib import pyplot as plt 
# import utm
import math
class MeasureList:
	def __init__(self, records):
		self.records=records
		self.length=len(records)
		self.start=0
	def reasonable_speed(self,startrecord,endrecord):
		dist=math.sqrt((startrecord[2]-endrecord[2])**2+(startrecord[3]-endrecord[3])**2)
		duration = (endrecord[1]-startrecord[1]).total_seconds()
		if 5<dist/duration<40:
			return True
		else:
			return False
	def duration(self,startrecord,endrecord):
		low=25
		high=65
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




def london_roadmap(boundingbox):
	start = time.time()

	conn = psycopg2.connect("dbname=gis user=postgres password=peach")
	cur = conn.cursor()
	envelope="ST_Transform(ST_MakeEnvelope(%s,%s,%s,%s,32630),4326)" % boundingbox
	# subgraph
	cur.execute("drop table if exists subgraph;")
	cur.execute("create table subgraph as (select * from hh_2po_4pgr where geom_way && "+envelope+");")
	# cur.execute("""create or replace function getedgeid(integer,integer) returns integer as 'select id from subgraph where ($1=source and $2=target) or ($1=target and $2=source) limit 1;' language sql immutable returns null on null input;""")


	sourcex = "ST_X(ST_transform(ST_setsrid(ST_makepoint(x1,y1),4326),32630))"
	sourcey = "ST_Y(ST_transform(ST_setsrid(ST_makepoint(x1,y1),4326),32630))"
	targetx = "ST_X(ST_transform(ST_setsrid(ST_makepoint(x2,y2),4326),32630))"
	targety = "ST_Y(ST_transform(ST_setsrid(ST_makepoint(x2,y2),4326),32630))"
	cur.execute("select source,target,"+sourcex+","+sourcey+","+targetx+","+targety+",km from subgraph;")

	roadmap = nx.Graph()
	for record in cur:
		source,target,x1,y1,x2,y2,cost=record
		# print source,target
		if not roadmap.has_node(source):
			roadmap.add_node(source,pos=numpy.array([x1,y1]))
		if not roadmap.has_node(target):
			roadmap.add_node(target,pos=numpy.array([x2,y2]))
		if source!=target:
			roadmap.add_edge(source,target,weight=cost*1000)
		# print roadmap.node[source]['pos'],roadmap.node[target]['pos'],cost*1000

	end =time.time()

	print "Graph generating ..."
	print end-start, 'secs'

	# position={}
	# for node in roadmap:
	# 	position[node]=roadmap.node[node]['pos']
	# nx.draw(roadmap,pos=position,node_size=50)
	# plt.show()

	print nx.number_of_nodes(roadmap), "nodes"
	print nx.number_of_edges(roadmap), "edges"
	return roadmap

def ecourier_data(vehicleid,halfmins,boundingboxratio):
	conn = psycopg2.connect("dbname=gis user=postgres password=peach")
	cur = conn.cursor()
	pointsqlx ="ST_X(ST_transform(ST_setsrid(ST_makepoint(longitude,latitude),4326),32630))"
	pointsqly ="ST_Y(ST_transform(ST_setsrid(ST_makepoint(longitude,latitude),4326),32630))"
	cur.execute("select "+pointsqlx+","+pointsqly+", timestamp from ecourier_oneday where vehicleid=%s order by timestamp",(vehicleid,))
	x,y,time0=cur.fetchone()
	minx,maxx,miny,maxy=x,x,y,y 

	print x,y

	z=[numpy.array([x,y])]
	t=[0.0]
	for i in range(halfmins):
		x,y,timex=cur.fetchone()

		print x,y

		z.append(numpy.array([x,y]))
		t.append((timex-time0).total_seconds())
		if x<minx: 
			minx=x
		elif x>maxx:
			maxx=x
		if y<miny:
			miny=y
		elif y>maxy:
			maxy=y 
	width=maxx-minx
	height=maxy-miny
	minx=minx-width*boundingboxratio/2.0
	maxx=maxx+width*boundingboxratio/2.0
	miny=miny-height*boundingboxratio/2.0
	maxy=maxy+height*boundingboxratio/2.0
	return (z,t,(minx,miny,maxx,maxy))

def ecourier_vertex_data(vehicleid,boundingboxratio):
	halfmins = 1
	conn = psycopg2.connect("dbname=gis user=postgres password=peach")
	cur = conn.cursor()
	pointsqlx ="ST_X(ST_transform(ST_setsrid(ST_makepoint(longitude,latitude),4326),32630))"
	pointsqly ="ST_Y(ST_transform(ST_setsrid(ST_makepoint(longitude,latitude),4326),32630))"
	cur.execute("select id,"+pointsqlx+","+pointsqly+", timestamp from near_vertex_measurement order by timestamp",(vehicleid,))
	pid,x,y,time0=cur.fetchone()

	endid=pid
	endpoint = numpy.array([x,y])
	endtime=time0
	while True:
		startid = endid
		startpoint = endpoint
		starttime = endtime
		pid,x,y,timex=cur.fetchone()
		endid = pid
		endpoint = numpy.array([x,y])
		endtime = timex
		if (endtime-starttime).total_seconds()<100:
			break
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
	return (z,t,(minx,miny,maxx,maxy),IDlist)

def reasonable_speed(z,obstime):
	dist = math.sqrt((z[0][0]-z[1][0])**2+(z[0][1]-z[1][1])**2)
	duration = obstime[-1]-obstime[0]
	if 5<dist/duration<40:
		return True
	else:
		return False

def reasonable_interval(obstime):
	smalltime=30
	bigtime=90
	if smalltime<obstime[-1]-obstime[0]<bigtime:
		return True
	else:
		return False

def ecourier_pair_data(vehicleid1,vehicleid2,boundingboxratio):
	# Fetch trajectories of two vehicles with a meeting potential, which means these two trajectories
	# has a timeral overlap and the boundingboxes of these two intersects.
	conn = psycopg2.connect("dbname=gis user=postgres password=peach")
	cur = conn.cursor()
	pointsqlx ="ST_X(ST_transform(geom,32630))"
	pointsqly ="ST_Y(ST_transform(geom,32630))"
	# cur.execute("select distinct vehicleid from near_vertex_measurement order by vehicleid;")
	# vehicles = [vehicle[0] for vehicle in cur]
	cur.execute("select id,timestamp,"+pointsqlx+","+pointsqly+" from near_vertex_measurement where vehicleid=%s and type='%s' order by timestamp;" % vehicleid1)
	# vehicle1_records=[measurement for measurement in cur]
	# vehicle1_len = len(vehicle1_records)
	vehicle1_records=MeasureList([measurement for measurement in cur])
	cur.execute("select id,timestamp,"+pointsqlx+","+pointsqly+" from near_vertex_measurement where vehicleid=%s and type='%s' order by timestamp;" % vehicleid2)
	# vehicle2_records=[measurement for measurement in cur]
	# vehicle2_len = len(vehicle2_records)
	haha=[measurement for measurement in cur]
	vehicle2_records=MeasureList([measurement for measurement in cur])
	z1,t1,id1=vehicle1_records.next()
	z2,t2,id2=vehicle2_records.next()
	while z1 and z2:
		if t1[0]<t2[1] and t2[0]<t1[1]:
			cur.execute("select trajectoryintersect(%s,%s,%s,%s);" % (id1[0],id1[1],id2[0],id2[1]))
			(intersect,) = cur.fetchone()
			if intersect:
				zerotime = min(t1+t2)
				t1 = [(t-zerotime).total_seconds() for t in t1]
				t2 = [(t-zerotime).total_seconds() for t in t2]
				minx = min([p[0] for p in z1+z2])
				maxx = max([p[0] for p in z1+z2])
				miny = min([p[1] for p in z1+z2])
				maxy = max([p[1] for p in z1+z2])
				width=maxx-minx
				height=maxy-miny
				minx=minx-width*boundingboxratio/2.0
				maxx=maxx+width*boundingboxratio/2.0
				miny=miny-height*boundingboxratio/2.0
				maxy=maxy+height*boundingboxratio/2.0
				return (z1,t1,z2,t2,(minx,miny,maxx,maxy))
			# else:
			# 	z1,t1,id1=vehicle1_records.next()
			# 	z2,t2,id2=vehicle2_records.next()
		if t1[0]<=t2[0]:
			z1,t1,id1=vehicle1_records.next()
		else:
			z2,t2,id2=vehicle2_records.next()
	return ([],[],[],[],(0,0,0,0))		
	# while i1<vehicle1_len-1 and i2<vehicle2_len-1:
	# 	v1_start=vehicle1_records[i1][1]
	# 	v1_end = vehicle1_records[i1+1][1]
	# 	v2_start=vehicle2_records[i2][1]
	# 	v2_end = vehicle2_records[i2+1][1]
	# 	if (v1_end-v1_start).total_seconds()<70 and (v2_end-v2_start).total_seconds()<70:
	# 		if v1_start<v2_end and v2_start<v1_end:
	# 			cur.execute("select trajectoryintersect(%s,%s,%s,%s);" % (vehicle1_records[i1][0],vehicle1_records[i1+1][0],vehicle2_records[i2][0],vehicle2_records[i2+1][0]))
	# 			(intersect,) = cur.fetchone()
	# 			z1 = [numpy.array([record[2],record[3]]) for record in (vehicle1_records[i1],vehicle1_records[i1+1])]
	# 			z2 = [numpy.array([record[2],record[3]]) for record in (vehicle2_records[i2],vehicle2_records[i2+1])]
	# 			zerotime = min([record[1] for record in (vehicle1_records[i1],vehicle1_records[i1+1],vehicle2_records[i2],vehicle2_records[i2+1])])
	# 			t1 = [(record[1]-zerotime).total_seconds() for record in (vehicle1_records[i1],vehicle1_records[i1+1])]
	# 			t2 = [(record[1]-zerotime).total_seconds() for record in (vehicle2_records[i2],vehicle2_records[i2+1])]
	# 			if intersect and reasonable_speed(z1,t1) and reasonable_speed(z2,t2):
	# 				minx = min([record[2] for record in (vehicle1_records[i1],vehicle1_records[i1+1],vehicle2_records[i2],vehicle2_records[i2+1])])
	# 				maxx = max([record[2] for record in (vehicle1_records[i1],vehicle1_records[i1+1],vehicle2_records[i2],vehicle2_records[i2+1])])
	# 				miny = min([record[3] for record in (vehicle1_records[i1],vehicle1_records[i1+1],vehicle2_records[i2],vehicle2_records[i2+1])])
	# 				maxy = max([record[3] for record in (vehicle1_records[i1],vehicle1_records[i1+1],vehicle2_records[i2],vehicle2_records[i2+1])])
	# 				width=maxx-minx
	# 				height=maxy-miny
	# 				minx=minx-width*boundingboxratio/2.0
	# 				maxx=maxx+width*boundingboxratio/2.0
	# 				miny=miny-height*boundingboxratio/2.0
	# 				maxy=maxy+height*boundingboxratio/2.0
	# 				return (z1,t1,z2,t2,(minx,miny,maxx,maxy))
	# 			else:
	# 				i1+=1
	# 				i2+=1
	# 		elif v1_end<=v2_start:
	# 			i1+=1
	# 		else:
	# 			i2+=1
	# 	else:
	# 		if (v1_end-v1_start).total_seconds()>=70:
	# 			i1+=1
	# 		if (v2_end-v2_start).total_seconds()>=70:
	# 			i2+=1
	# return ([],[],[],[],(0,0,0,0))

