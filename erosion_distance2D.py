import numpy as np
import sys


def erosion(f, tfx, tfy, g, tgx, tgy, spacing, xw, yw):
    nt = tfx + tfy + tgx + tgy
    newtimes = np.array(nt)
    newtimes = np.unique(newtimes)
    newtimes = np.sort(newtimes)
    ntx = np.sort(np.unique(np.array(tfx + tgx)))
    nty = np.sort(np.unique(np.array(tfy + tgy)))

    f[f == -1] = np.inf
    g[g == -1] = np.inf

    a = xw
    xw = -yw
    yw = a

    for i in range(len(ntx)-1):
        if ntx[i+1]-ntx[i] != spacing:
            print('Please generate modules whose discretization grids are not offset by anything other than a multiple of the spacing.')
            exit()
    for i in range(len(nty)-1):
        if nty[i+1]-nty[i] != spacing:
            print('Please generate modules whose discretization grids are not offset by anything other than a multiple of the spacing.')

    if tfx != tgx or tfy != tgy:
        newf = expand2D(f, tfx, tfy, ntx, nty)
        newg = expand2D(g, tgx, tgy, ntx, nty)
    else:
        newf = f
        newg = g

    fshifts = np.zeros((len(ntx), len(nty)), dtype=int)
    gshifts = np.zeros((len(ntx), len(nty)), dtype=int)

    for i in range(len(ntx)):
        for j in range(len(nty)):
            fshifts[i, j] = elementshift(newf, newg, i, j, xw, yw)
            gshifts[i, j] = elementshift(newg, newf, i, j, xw, yw)

    needshifts = np.maximum(fshifts, gshifts)
    d = np.amax(needshifts) * spacing

    return d


def expand2D(f, tfx, tfy, ntx, nty):
    nx = len(ntx)
    ny = len(nty)

    mx = len(tfx)
    my = len(tfy)
    newf = np.zeros((nx, ny))
    currentx = 0

    for k in range(nx):
        if currentx < mx-1:
            while tfx[currentx + 1] <= ntx[k] and currentx < mx-1:
                currentx += 1
                if currentx == mx - 1:
                    break

        currenty = 0
        for i in range(ny):
            if currenty >= my-1:
                newf[k,i] = f[currentx, currenty]
            else:
                while tfy[currenty + 1] <= nty[i] and currenty < my - 1:
                    currenty += 1
                    if currenty == my - 1:
                        break
                if f[currentx,currenty] == np.inf:
                    newf[k,i] = np.inf
                else:
                    newf[k,i] = f[currentx, currenty]

    return newf


def elementshift(f, g, i2, j2, xw, yw):
    [i1, j1] = np.shape(f)

    if xw<0:
        xr = i2
    elif xw>=0:
        xr = i1-i2-1
    if yw<0:
        yr = j2
    elif yw>=0:
        yr = j1-j2-1

    r = min([np.floor(xr/abs(xw)), np.floor(yr/abs(yw))])
    r = max([r, 0])
    r = int(r)

    if f[i2, j2] >= g[i2, j2]:
        return 0
    if f[i2, j2] < g[i2 + (r * xw), j2 + (r * yw)]:
        return r + 1

    l = 0
    current = int(np.floor((r + l) / 2))
    while r - l > 1:
        if f[i2, j2] >= g[i2 + (current * xw), j2 + (current * yw)]:
            r = current
            current = int(np.floor((r + l) / 2))
        else:
            l = current
            current = int(np.floor((r + l) / 2))

    if f[i2, j2] >= g[i2 + (l * xw), j2 + (l * yw)]:
        return l
    else:
        return r


if __name__ == "__main__":
    if len(sys.argv) != 3 and len(sys.argv) != 5:
        print('Please input two filenames of 2D modules, or two filenames and integer weights for x and y values of search direction. (Default is x=-1, y=1)')
        exit()
    elif len(sys.argv) == 3:
        filename1 = sys.argv[1]
        filename2 = sys.argv[2]
        x_weight = -1
        y_weight = 1
    elif len(sys.argv) == 5:
        filename1 = sys.argv[1]
        filename2 = sys.argv[2]
        x_weight = int(sys.argv[3])
        y_weight = int(sys.argv[4])



    with open(filename1, 'r') as f1:
        firstline1 = f1.readline()
        info1 = firstline1.split()
        start_thresh1 = int(info1[0])
        end_thresh1 = int(info1[1])
        spacing1 = int(info1[2])

    with open(filename2, 'r') as f2:
        firstline2 = f2.readline()
        info2 = firstline2.split()
        start_thresh2 = int(info2[0])
        end_thresh2 = int(info2[1])
        spacing2 = int(info2[2])

    if spacing1 != spacing2:
        print("The current algorithm only works when the spacing in the discretization of the modules is constant. Please re-sample the modules so this is the case.")
        exit()

    F1 = np.loadtxt(filename1, skiprows=1)
    F2 = np.loadtxt(filename2, skiprows=1)

    samples1 = [value for value in range(start_thresh1,end_thresh1+1,spacing1)]
    samples2 = [value for value in range(start_thresh2,end_thresh2+1,spacing2)]

    if samples1[-1]>end_thresh1:
        samples1 = samples1[0:-1]
    if samples2[-1]>end_thresh2:
        samples2 = samples2[0:-1]

    print(erosion(F1, samples1, samples1, F2, samples2, samples2, spacing1, x_weight, y_weight))

