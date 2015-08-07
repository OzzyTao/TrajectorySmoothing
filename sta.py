import csv
import pickle
import operator
from scipy.stats import pearsonr

def compare(col):
	def temp(x,y):
		if x[col]==y[col]:
			return x['possibility_Rank']-y['possibility_Rank']
		else:
			if col=='CLength' or col=='CEdges':
				return 1 if y[col]-x[col]>0 else -1
			else:
				return 1 if x[col]-y[col]>0 else -1
	def possi(x,y):
		return 1 if y[col]-x[col]>0 else -1
	if col=='possibility':
		return possi
	else:
		return temp
class Test:
	def __init__(self,records):
		self.ranked = False
		self.routes = []
		for record in records:
			self.routes.append({
				'route':record['route'],
				'edgeIDList':record['edgeIDList'],
				'possibility':float(record['possibility']),
				'RMSE':float(record['RMSE']),
				'TRMSE':float(record['temporal_RMSE']),
				'area':float(record['area']),
				'CLength':float(record['common_length']),
				'CEdges':float(record['common_edges'])
					})
		self.ranking()
		self.corr()
	def best_ranking(self,key):
		if self.ranked:
			self.routes.sort(key=operator.itemgetter(key+'_Rank'))
		else:
			self.routes.sort(cmp=compare(key))
		return self.routes[0]
	def ranking(self):
		keys = ['possibility','RMSE','TRMSE','area','CLength','CEdges']
		for key in keys:
			self.best_ranking(key)
			rank = 1
			for route in self.routes:
				route[key+'_Rank']=rank
				rank+=1
		self.ranked = True
	def corr(self):
		index=['RMSE','TRMSE','area','CLength','CEdges']
		index_Rank = [e+'_Rank' for e in index]
		for i in index:
			setattr(self,i+'_corr',pearsonr([route[i] for route in self.routes],[route['possibility'] for route in self.routes]))
		for i in index_Rank:
			setattr(self,i+'_corr',pearsonr([route[i] for route in self.routes],[route['possibility_Rank'] for route in self.routes]))




if __name__ == "__main__":
	filename="60s101test.csv"
	with open(filename,'rb') as csvfile:
		csv_reader = csv.DictReader(csvfile,fieldnames=[],restkey='undefined-fieldnames',delimiter=',')
		current_row = 0
		tests = []
		singletest_records = []
		test_id = '1'
		for row in csv_reader:
			current_row+=1
			if current_row == 1:
				csv_reader.fieldnames = row['undefined-fieldnames']
			else:
				# print current_row,row
				if row['TestID'] == test_id:
					singletest_records.append(row)
				else:
					tests.append(Test(singletest_records))
					print "Finish test", test_id
					test_id = row['TestID']
					singletest_records=[row]
		tests.append(Test(singletest_records))
		with open("f:/London/statistics/60s.p",'wb') as binaryfile:
			pickle.dump(tests,binaryfile)

