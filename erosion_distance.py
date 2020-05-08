import numpy as np
import sys

def erosion(f,tfx,tfy,tfz,g,tgx,tgy,tgz,spacing):
	nt = tfx + tfy + tfz + tgx + tgy + tgz
	newtimes = np.array(nt)
	newtimes = np.unique(newtimes)
	newtimes = np.sort(newtimes)
	ntx = np.sort(np.unique(np.array(tfx + tgx)))
	nty = np.sort(np.unique(np.array(tfy + tgy)))
	ntz = np.sort(np.unique(np.array(tfz + tgz)))

	if tfx != tgx:
		newf = expand3D(f,tfx,tfy,tfz,ntx,nty,ntz)
		newg = expand3D(g,tgx,tgy,tgz,ntx,nty,ntz)
	else:
		newf = f
		newg = g

	fshifts = np.zeros((len(ntx),len(nty),len(ntz)), dtype=int)
	gshifts = np.zeros((len(ntx),len(nty),len(ntz)), dtype=int)
	
	for i in range(len(ntx)):
		for j in range(len(nty)):
			for k in range(len(ntz)):
				fshifts[i,j,k]=elementshift(newf,newg,i,j,k)
				gshifts[i,j,k]=elementshift(newg,newf,i,j,k)
		
	needshifts = np.maximum(fshifts,gshifts)
	d = np.amax(needshifts) * spacing

	return d


def expand3D(f, tfx, tfy, tfz, ntx, nty, ntz):
	nx = len(ntx)
	ny = len(nty)
	nz = len(ntz)

	mx = len(tfx)
	my = len(tfy)
	mz = len(tfz)
	newf = np.zeros((nx, ny, nz), dtype=int)
	currentx = 0

	for k in range(nx):
		if currentx < mx - 1:
			while tfx[currentx + 1] <= ntx[k] and currentx < mx - 1:
				currentx += 1
				if currentx == mx - 1:
					break

		currenty = 0
		for i in range(ny):
			if currenty < my - 1:
				while tfy[currenty + 1] <= nty[i] and currenty < my - 1:
					currenty += 1
					if currenty == my - 1:
						break

			currentz = 0
			for j in range(nz):
				if currentz >= mz - 1:
					newf[k, i, j] = f[currentx, currenty, mz - 1]
				else:
					while tfz[currentz + 1] <= ntz[j] and currentz < mz - 1:
						currentz += 1
						if currentz == mz - 1:
							break
					newf[k, i, j] = f[currentx, currenty, currentz]

	return newf


def elementshift(f,g,i,j,k):
	[i1, j1, k1] = np.shape(f)
	r = min([i1-i-1, j-1, k-1])
	r = max([r, 0])
	if f[i,j,k] >= g[i,j,k] or f[i,j,k] == 0:
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
	if len(times1) != len(times2):
		print("The current algorithm only works when the same number of time samples are used for the two betti-0 functions. Please re-generate the betti-0 functions so this is the case.")

	F1 = np.loadtxt(filename1, skiprows=1)
	F1 = F1.reshape((len(thresholds1),len(times1),len(times1)))
	F2 = np.loadtxt(filename2, skiprows=1)
	F2 = F2.reshape((len(thresholds2),len(times2),len(times2)))
	
	print(erosion(F1,thresholds1,times1,times1,F2,thresholds2,times1,times1,spacing1))

