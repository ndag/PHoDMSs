import numpy as np
import time
import multiprocessing as mp
import sys

def erosion(f,tfx,tfy,tfz,g,tgx,tgy,tgz,spacing):
	nt = tfx + tfy + tfz + tgx + tgy + tgz
	newtimes = np.array(nt)
	newtimes = np.unique(newtimes)
	newtimes = np.sort(newtimes)
	n = len(newtimes)	
	
	if len(tfx) != len(tgx) or len(tfy) != len(tgy) or len(tfz) != len(tgz):
		newf = expand3D(f,tfx,tfy,tfz,newtimes)
		newg = expand3D(g,tgx,tgy,tgz,newtimes)
	else:
		newf = f
		newg = g

	fshifts = np.zeros((n,n,n), dtype=int)
	gshifts = np.zeros((n,n,n), dtype=int)
	
	for i in range(0,len(tfx)):
		for j in range(0,len(tfy)):
			for k in range(0,len(tfz)):
				fshifts[i,j,k]=elementshift(newf,newg,i,j,k)
				gshifts[i,j,k]=elementshift(newg,newf,i,j,k)
		
	needshifts = np.maximum(fshifts,gshifts)
	d = np.amax(needshifts) * spacing

	return d


def expand3D(f, tfx, tfy, tfz, nt):
	n = len(nt)

	mx = len(tfx)
	my = len(tfy)
	mz = len(tfz)
	
	newf = np.zeros((n,n,n), dtype=int)
	currentz = -1
	for k in range(n):
		if currentz < mz-1:
			if tfz[currentz + 1] == nt[k]:
				currentz += 1
						
		currenty = -1
		for i in range(n):
			if currenty < my-1:
				if tfy[currenty + 1] == nt[i]:
					currenty += 1

			currentx = -1
			for j in range(n):
				if currentx >= mx-1:
					newf[j,i,k] = f[mx-1,currenty,currentz]
				elif tfx[currentx+1] == nt[j]:
					currentx +=1
					newf[j,i,k] = f[currentx,currenty,currentz]
				else:
					newf[j,i,k] = f[currentx,currenty,currentz]
					

	return newf


def elementshift(f,g,i,j,k):
	[i1, j1, k1] = np.shape(f)
	r = min([i1-i-1, j-1, k-1])
	r = max([r, 0])
	if f[i,j,k] >= g[i,j,k]:
		return 0
	if f[i,j,k] < g[i+r,j-r,k-r]:
		return r+1
	
	l = 0
	current = int(np.floor((r+l)/2))
	while r-l > 1:
		if f[i,j,k] >= g[i+current,j-current,k-current]:
			r = current
			current = int(np.floor((r+l)/2))
		else:
			l = current
			current = int(np.floor((r+l)/2))

	if f[i,j,k] >= g[i+l,j-l,k-l]:
		return l
	else:
		return r

if __name__ == "__main__":
	if len(sys.argv) == 3:
		filename1 = sys.argv[1]
		filename2 = sys.argv[2]
	else:
		print('please input two filenames of betti-0 functions')
		exit()

	with open(filename1, 'r') as f1:
		firstline1 = f1.readline()
		info1 = firstline1.split()
		num_times1 = int(info1[0])
		start_thresh1 = int(info1[1])
		end_thresh1 = int(info1[2])
		spacing1 = int(info1[3])
		thresholds1 = [i for i in range(start_thresh1,end_thresh1+1,spacing1)]
		times1 = [i for i in range(start_thresh1, start_thresh1 + spacing1*num_times1, spacing1)]

	with open(filename2, 'r') as f2:
		firstline2 = f2.readline()
		info2 = firstline2.split()
		num_times2 = int(info2[0])
		start_thresh2 = int(info2[1])
		end_thresh2 = int(info2[2])
		spacing2 = int(info2[3])
		thresholds2 = [i for i in range(start_thresh2,end_thresh2+1,spacing2)]
		times2 = [i for i in range(start_thresh2, start_thresh2 + spacing2*num_times2, spacing2)]
	
	if spacing1 != spacing2:
		print("The current algorithm only works when the threshold parameters for the two betti-0 functions have the same spacing. Please re-generate the betti-0 functions so this is the case.")
		exit()
	
	F1 = np.loadtxt(filename1, skiprows=1)
	F1 = F1.reshape((len(thresholds1),len(times1),len(times1)))
	F2 = np.loadtxt(filename2, skiprows=1)
	F2 = F2.reshape((len(thresholds2),len(times2),len(times2)))
	
	print(erosion(F1,thresholds1,times1,times1,F2,thresholds2,times2,times2,spacing1))

