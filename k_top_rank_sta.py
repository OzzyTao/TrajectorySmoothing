from sta import Test
import pickle
import csv
path = 'f:/London/statistics/'
names = ['RMSE','TRMSE','area','CLength','RMSE_Rank','TRMSE_Rank','area_Rank','CLength_Rank']
# K values
fields = ['10','20','30','50','100','150'] 
suffix = 'econds101test.p'
binaryfiles = ['60seconds115test10k.p','60seconds101test20k.p','60seconds115test30k.p','60seconds115test50k.p','60seconds101test.p','60seconds115test150k.p']
typetest = []
for file in binaryfiles:
	with open(path+file,'rb') as binary:
		typetest.append(pickle.load(binary))

for name in names:
	with open(path+'k/'+name+'.csv','wb') as csvfile:
		mywriter = csv.writer(csvfile,delimiter=',')
		mywriter.writerow(fields)
		for testid in range(100):
			tmpvalues = []
			for tests in typetest:
				top_rank = tests[testid].best_ranking('possibility')
				tmpvalues.append(top_rank[name])
			mywriter.writerow(tmpvalues)