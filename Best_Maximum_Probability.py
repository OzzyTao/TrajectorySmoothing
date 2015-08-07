from sta import Test
import pickle
import csv
path = 'f:/London/statistics/'
names = ['RMSE','TRMSE','area','CLength','RMSE_Rank','TRMSE_Rank','area_Rank','CLength_Rank']
# K values
fields = ['30s','40s','50s','60s','90s','120s','150s'] 
seconds = [30,40,50,60,90,120,150]
suffix = 'econds101test.p'
binaryfiles=[field+suffix for field in fields[:-1]]+['150seconds111test.p']
# binaryfiles = ['60seconds115test10k.p','60seconds101test20k.p','60seconds115test30k.p','60seconds115test50k.p','60seconds101test.p','60seconds115test150k.p']
typetest = []
for file in binaryfiles:
	with open(path+file,'rb') as binary:
		typetest.append(pickle.load(binary))

for name in names[:4]:
	with open(path+'best_rank/'+name+'_max.csv','wb') as csvfile:
		mywriter = csv.writer(csvfile,delimiter=',')
		mywriter.writerow(['seconds','rank','probability'])
		for i in range(7):
			second = seconds[i]
			tests = typetest[i]
			for testid in range(100):
				best_rank=tests[testid].best_ranking(name)
				top_rank=tests[testid].best_ranking('possibility')
				mywriter.writerow((second,best_rank['possibility_Rank'],top_rank['possibility']))
