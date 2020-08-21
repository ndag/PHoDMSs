import numpy as np
from scipy.spatial.distance import pdist, squareform
from functools import partial
from itertools import combinations
import multiprocessing as mp
import sys
import os
import itertools as it
import dionysus as d

def get_dist(points):
	n = len(points)
	xpts = np.zeros((n,2))
	ypts = np.zeros((n,2))
	points = np.array(points)

	xpts[:,0] = points[:,0]
	ypts[:,0] = points[:,1]
	xdists = pdist(xpts)
	ydists = pdist(ypts)
	
	dists = np.sqrt(np.square(xdists) + np.square(ydists))
	return dists


def get_points(file, numpoints):
	points = []
	for i in range(numpoints):
		line = file.readline()
		sline = line.split(' ')
		fline = list(map(float,sline))
		points.append(fline)
	return points


def get_rips(ripsfilt, threshold):
	simplices = []
	for simplex in ripsfilt:
		if simplex.data <= threshold:
			simplices.append(tuple(simplex))
	return simplices

def compute_rank(simpliceslist, dim, s1, t1, t2, s2, t3, t4):
	f = d.Filtration()
	first = simpliceslist[s1][t1][t2]
	second = simpliceslist[s2][t3][t4]
	if first is None or second is None:
		return 0
	second = set(second)-set(first)
	if len(second) > 0:
		second = list(second)
	for simplex in first:
		f.append(d.Simplex(simplex, 0))
	for simplex in second:
		f.append(d.Simplex(simplex, 1))
	m = d.homology_persistence(f)
	dgms = d.init_diagrams(m,f)
	if len(dgms) > dim:
		return len(dgms[dim])
	else:
		return 0


def rank_indices(simpliceslist, dim, m, n, s1, t1, t2, s2, t3, t4):
	if t3<=t1<=t2<=t4 and s1<=s2:
		result = compute_rank(simpliceslist, dim, s1, t1, t2, s2, t3, t4)
		return result
	else:
		r = min([t1, n-t2-1, m-s1-1, n-t3-1, t4, s2])
		for i in range(r):
			if t3+i <=t1-i <=t2+i <=t4-i and s1+i <=s2-i:
				return 0
		return -1
			


if __name__ == '__main__':
	if len(sys.argv) == 8:
		position_file = str(sys.argv[1])
		rank_file = str(sys.argv[2])
		start_thresh = int(sys.argv[3])
		end_thresh = int(sys.argv[4])
		spacing = int(sys.argv[5])
		time_samples = int(sys.argv[6])
		dim = int(sys.argv[7])
	else:
		print('Please pass arguments:')
		print('position filename, output filename, starting threshold, ending threshold, spacing, time saples, dimension')
		exit()

	points = []
	num_times = -1
	for line in open(position_file).readlines():
		num_times += 1

	myfile = open(position_file)
	line1 = myfile.readline()
	num_points = int(line1)
	num_times = int(num_times / num_points)
	for i in range(num_times):
		points.append(get_points(myfile, num_points))
	myfile.close()
	
	thresholds = [i for i in range(start_thresh, end_thresh + 1, spacing)]
	

	useddists = [get_dist(points[i]) for i in range(0, num_times, int(num_times / (time_samples - 1))-1)]
	useddists = np.array(useddists)
	
	n = len(useddists)
	m = len(thresholds)
	
	simpliceslist = [[[None for i in range(n)] for j in range(n)] for k in range(m)]
	
	for i in range(n):
		f = d.fill_rips(useddists[i], dim + 1, thresholds[-1])
		for j in range(m):
			simpliceslist[j][i][i] = get_rips(f, thresholds[j])	

	for diag in range(n):
		for row in range(n-diag-1):
			newuseddists = np.minimum(useddists[row], useddists[row + 1])
			useddists[row] = newuseddists
			f = d.fill_rips(newuseddists, dim+1, thresholds[-1])
			for j in range(m):
				simpliceslist[j][row][row+diag+1] = get_rips(f, thresholds[j])
						

	result2 = np.zeros((m,n,n,m,n,n))
	
	indiceslist = list(it.product(range(m), range(n), range(n), range(m), range(n), range(n)))

	mp_rank = partial(rank_indices, simpliceslist, dim, m, n)
	num_cores = os.cpu_count()
	
	
	pool = mp.Pool(num_cores)
	result = pool.starmap(mp_rank, indiceslist)
	
	for i in range(len(indiceslist)):
		result2[indiceslist[i]] = result[i]	

	result2 = np.reshape(result2, (m*m*n*n*n, n))

	infostring = str(m) + " " + str(n) + " " + str(dim) + " " + str(spacing)
	with open(rank_file, 'w') as outfile:
		outfile.write(infostring +'\n')
		np.savetxt(outfile, result2, fmt='%d')
	

	
