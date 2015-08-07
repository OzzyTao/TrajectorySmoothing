import sys
import pickle
import math
import numpy
import testdb
import time
import case1
import csv
from sta import Test
import networkx as nx
from case1_main import cal_min_bbox, nearest_vertex, get_edges
sys.path.append('F:/London/map-matching/SnoT_Matching')
import psycopg2
from map_matching import cal_route
sys.path.append('F:/London/map-matching/SnoT_Matching/mm_test')
from mm_test import RSME, temporal_RSME, surrounding_area, share_distance, speed

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

def get_allmeasurements(measurements,vehicleid,cur):
	allmeasurements=[]
	cur.execute("select * from gettrajectory(%s, '%s', %s, %s);" % (vehicleid[0],vehicleid[1], measurements[0], measurements[-1]))
	for record in cur:
		allmeasurements.append(record)
	return allmeasurements

# def NormalisedTemporal(allmeasurements,testid,cur):

def common_edge(realedges,edges):
	common = 0
	if realedges:
		for edge in edges:
			if edge in realedges:
				common+=1
		return common*1.0/len(realedges)
	return common

def case1_vehicle(vehicleid,interval):
	boundingboxratio=0.5
	sig = 5
	# Algorithm parameters
	n = 200  # number of samples
	K = 50  # minimum number of candidate paths
	vmean = 6. # mean velocity
	vvar = 4. # velocity variance
	includetruepath=False# without inputting true path

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
		allmeasurements_id = [measurement[0] for measurement in allmeasurements]
		addrow_query = "INSERT INTO testcases3(time_slot, start_num, end_num, measurements) VALUES (%s, %s, %s, %s) RETURNING test_id;"
		measurementsstr = "{%s}" % (','.join([str(m_id) for m_id in allmeasurements_id]),)
		cur.execute(addrow_query,(t[-1],idl[0],idl[-1],measurementsstr))
		testid = cur.fetchone()[0]
		conn.commit()

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
		if nx.has_path(G,startvertex,endvertex) and startvertex!=endvertex:
			try:
				t0 = time.time()
				tedge,tpath = cal_route(allmeasurements_id,testid,cur)
				updateroute = "{%s}" % (','.join(tedge),)
				updatevertex = "{%s}" % (','.join([str(x) for x in tpath]),)
				cur.execute('''update testcases3 set gt_mm_edges=%s where test_id = %s''',(updateroute,testid))
				cur.execute('''update testcases3 set gt_mm_vertex=%s where test_id = %s''',(updatevertex,testid))
				conn.commit()
				[pp,sep,wp] = case1.calcpostprobs(G,zl,t,n,K,vmean,vvar,sig*sig,tpath if includetruepath else None)
				et = time.time()-t0
				print et, "secs, done algorithm"
				possibleroutes = len(pp)
				# print pp
				if allmeasurements and allmeasurements[0][0]==idl[0] and allmeasurements[-1][0]==idl[-1]:
					result=[]
					for j in range(possibleroutes):
						# print "modeled path:",pp[j]
						edge = get_edges(pp[j],cur)
						sarea = surrounding_area(edge,pp[j][0],cur,testid)

						temproute = {'route':pp[j],
							'possibility':wp[j],
							'RMSE':RSME(edge,allmeasurements,cur),
							'temporal_RMSE':temporal_RSME(edge,allmeasurements,pp[j][0],cur,allmeasurements_id,testid),
							'common_length':share_distance(edge,pp[j][0],cur,testid),
							'common_edges':common_edge(tedge,edge),
							'area':sarea,
							'edgeIDList':edge
							}
						result.append(temproute)
					newtest = Test(result)
					newtest.test_para = {'period':t[-1],
						'measurements':allmeasurements_id,
						"tpath":tpath,
						"tedge":tedge,
						"vehicleid":vehicleid,
						"test_id":testid,
						'v':len(G.nodes()),
						'e':len(G.edges())}
					results.append(newtest)
					# conn.commit()
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

def profile_main(interval):
	conn = psycopg2.connect("dbname=gis user=postgres password=peach")
	cur = conn.cursor()
	cur.execute("select distinct vehicleid, type from near_vertex_measurement;")
	results=[]
	for vehicleid in cur:
		results+=case1_vehicle(vehicleid,interval)
		if len(results)>=100:
			break
	# results=case1_vehicle(8)
	f = open("f:/London/statistics17/"+str(interval)+'seconds'+str(len(results))+'test.p','wb')
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
for timegroup in [30,60,90,120,150,180]:
	profile_main(timegroup)
# profile_main(180)