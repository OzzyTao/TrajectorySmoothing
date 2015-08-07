from sta import Test
import pickle
import csv
path = 'f:/London/statistics/'
names = ['RMSE','TRMSE','area']
# normalized indices: top ranked value against the worest_value 
fields = ['30s','40s','50s','60s','90s','120s','150s'] 
suffix = 'econds101test.p'
binaryfiles = [x+suffix for x in fields[:-1]]+['150seconds111test.p']
# binaryfiles = ['60seconds115test10k.p','60seconds101test20k.p','60seconds115test30k.p','60seconds115test50k.p','60seconds101test.p','60seconds115test150k.p']
typetest = []
for file in binaryfiles:
	with open(path+file,'rb') as binary:
		typetest.append(pickle.load(binary))

for name in names:
	with open(path+'top_rank/'+name+'_norm.csv','wb') as csvfile:
		mywriter = csv.writer(csvfile,delimiter=',')
		mywriter.writerow(fields)
		for testid in range(100):
			tmpvalues = []
			for tests in typetest:
				top_rank = tests[testid].best_ranking('possibility')
				tests[testid].best_ranking(name);
				worest_val = tests[testid].routes[-1][name]
				tmpvalues.append(top_rank[name]/worest_val)
			mywriter.writerow(tmpvalues)