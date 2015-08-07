from sta import Test
import pickle
import csv
path = 'f:/London/statistics/'
names = ['RMSE_corr','TRMSE_corr','area_corr','CLength_corr','RMSE_Rank_corr','TRMSE_Rank_corr','area_Rank_corr','CLength_Rank_corr']
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
			mywriter.writerow([getattr(tests[testid],name)[0] for tests in typetest])

