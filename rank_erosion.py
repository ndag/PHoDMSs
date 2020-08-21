import numpy as np
import sys

def erosion(f, g, spacing, m, n):
	fshifts = np.zeros((m,n,n,m,n,n), dtype=int)
	gshifts = np.zeros((m,n,n,m,n,n), dtype=int)

	for s1 in range(m):
		for t1 in range(n):
			for t2 in range(n):
				for s2 in range(m):
					for t3 in range(n):
						for t4 in range(n):
							fshifts[s1,t1,t2,s2,t3,t4] = elementshift(f,g,s1,t1,t2,s2,t3,t4)
							gshifts[s1,t1,t2,s2,t3,t4] = elementshift(g,f,s1,t1,t2,s2,t3,t4)	
	
	d = np.amax(np.maximum(fshifts,gshifts)) * spacing

	return d


def elementshift(f,g,s1,t1,t2,s2,t3,t4):
	[m1,n1,n1,m1,n1,n1] = np.shape(f)
	
	r = min([s1, m1-s2-1, n1-t1-1, t2, t3, n1-t4-1]) 
	r = int(max([r, 0]))
	
	l = 0
	current = int(np.floor(r+l)/2)
	while r-l > 1:
		if f[s1,t1,t2,s2,t3,t4] >= g[s1-current, t1+current, t2-current, s2+current, t3-current, t4+current]:
			r = current
			current = int(np.floor((r+l)/2))
		else:
			l = current
			current = int(np.floor((r+l)/2))

	if f[s1,t1,t2,s2,t3,t4] >= g[s1-l, t1+l, t2-l, s2+l, t3-l, t4+l]:
		return l
	elif r==0:
		return 1
	else:
		return r


if __name__ == "__main__":
	if len(sys.argv) == 3:
		filename1 = sys.argv[1]
		filename2 = sys.argv[2]
	else:
		print('please input two filenames of rank functions')
		exit()

	with open(filename1, 'r') as f1:
		firstline1 = f1.readline()
		info1 = firstline1.split()
		m1 = int(info1[0])
		n1 = int(info1[1])
		dim1 = int(info1[2])
		spacing1 = int(info1[3])

	with open(filename2, 'r') as f2:
		firstline2 = f2.readline()
		info2 = firstline2.split()
		m2 = int(info2[0])
		n2 = int(info2[1])
		dim2 = int(info2[2])
		spacing2 = int(info2[3])

	if m1 != m2: 
		print('Two rank functions are incompatible for comparison')
		exit()
	if n1 != n2:
		print('Two rank functions are incompatible for comparison, differing number of time samples')
		exit()
	if dim1 != dim2:
		print('Error: Two rank functions are for homologies in different dimension')
		exit()
	if spacing1 != spacing2:
		print('Two rank functions are incompatible for comparison, differing spacing of thresholds')
		exit()

	F1 = np.loadtxt(filename1, skiprows=1)
	F2 = np.loadtxt(filename2, skiprows=1)
	m = m1
	n = n1

	F1 = np.reshape(F1, (m,n,n,m,n,n))
	F2 = np.reshape(F2, (m,n,n,m,n,n))
	
	F1[F1 == -1] = np.inf
	F2[F2 == -1] = np.inf
	F1[F1 == -2] = 0
	F2[F2 == -2] = 0
	
	
	result = erosion(F1, F2, spacing1, m, n)
	print(result)
	
