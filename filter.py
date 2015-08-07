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
thres = [15,30,45]
for th in thres:
	with open(path+'filtered/bestroute'+str(th)+'.csv','wb') as csvfile:
		mywriter = csv.writer(csvfile,delimiter=',')
		mywriter.writerow(['seconds','RMSE','TRMSE','area','CLength'])
		for i in range(7):
			tests=typetest[i]
			second = seconds[i]
			for testid in range(100):
				top_rank = tests[testid].best_ranking('possibility')
				if top_rank['possibility']>th/100.0:
					tmprow = [second]
					for name in names[:4]:
						best_rank = tests[testid].best_ranking(name)
						tmprow.append(best_rank['possibility_Rank'])
					mywriter.writerow(tmprow)
