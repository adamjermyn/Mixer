import numpy as np
import sys
import table
import glob
import os

def parse(fname):
	fi = open(fname,'r')

	params = []
	data = []
	err = []

	for i,line in enumerate(fi):
		line = line.rstrip()

		if i%18 == 0:
			data.append([])
			err.append([])
			params.append(map(float,line.split(' ')))		
		if i%18 > 1 and line != '' and i%18 < 9:
			data[-1].append(map(float,line.split(' ')))
		if i%18 > 9 and line != '' and i%18 < 17:
			err[-1].append(map(float,line.split(' ')))

	params = np.array(params)
	data = np.array(data)
	err = np.array(err)

	return params, data, err

def parseAll():
	os.chdir('/data/vault/asj42')
	params = np.zeros((0,6))
	data = np.zeros((0,7,7))
	err = np.zeros((0,7,7))
	for file in glob.glob("*.out"):
		p, d, e = parse(file)
		params = np.append(params, p, axis=0)
		data = np.append(data, d, axis=0)
		err = np.append(err, e, axis=0)
		print(params.shape, p.shape, data.shape, d.shape)
	return params, data, err

params, data, err = parseAll()
import pickle
pickle.dump([params, data, err],open('dumpP.dat','w+'))


