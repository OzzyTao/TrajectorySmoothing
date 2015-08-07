import csv
import pickle
from sta import Test
file = "f:/London/statistics/60seconds101test.p"
outputfile = "verify.csv"
with open(file,'rb') as binaryfile:
	tests = pickle.load(binaryfile)
	test = tests[0]
	keys = test.routes[0].keys()

with open(outputfile,'wb') as output:
	mywriter = csv.DictWriter(output, delimiter = ',', fieldnames = keys)
	mywriter.writerow(dict((fn,fn) for fn in keys))
	for row in test.routes:
		mywriter.writerow(row)