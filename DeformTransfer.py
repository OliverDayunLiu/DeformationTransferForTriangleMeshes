import numpy as np
import scipy.sparse
import scipy.sparse.linalg
import os
import re
import sys
import datetime

# given three vertices, computes the fourth.
# assume each input vertex has size (3,)
def compute4thvertex(v1, v2, v3):
    vec1 = v2 - v1
    vec2 = v3 - v1
    perp = np.cross(vec1, vec2)
    perplength = np.sqrt(np.sum(np.power(perp, 2)))
    v4 = np.add(v1, perp/float(perplength))
    return v4

def readObj(file):
    data_list = [l.strip('\n') for l in open(file).readlines()]
    for line in data_list:
        if line[7:13] == "format":
            nums = map(int, re.findall(r'\d+', line))
            V, M = nums[0], nums[1]
            break
    vertices = np.zeros((V, 3))
    indices = np.zeros((M, 3)) #only contains vertex index
    indices = indices.astype(int)
    triangles2 = np.zeros((M, 4, 3)) #M triangles, each has 4 vertices (one will be constructed), each vertex has xyz values
    vindex, findex = 0, 0
    for line in data_list:
        if line[0] == "v" and line[1] != "n":
            linearr = line.split(" ")
            secondindex = 0
            for token in linearr:
                if token != "v" and token != "":
                    vertices[vindex, secondindex] = float(token)
                    secondindex += 1
            vindex += 1
        elif line[0] == "f":
            value = map(int, re.findall(r'\d+', line))
            indices[findex, 0] = value[0]
            indices[findex, 1] = value[2]
            indices[findex, 2] = value[4]
            findex += 1

    for m in range(0, M):
        triangles2[m, 0, :] = vertices[indices[m, 0]-1, :]
        triangles2[m, 1, :] = vertices[indices[m, 1]-1, :]
        triangles2[m, 2, :] = vertices[indices[m, 2]-1, :]
        triangles2[m, 3, :] = compute4thvertex(triangles2[m, 0, :], triangles2[m, 1, :], triangles2[m, 2, :])

    return V, M, vertices, triangles2, indices

def matrixIndexToXindex(indices, m):
    #for mth triangle, what's the indices (base 0) of its first 3 vertices?
    ret = []
    v1index = indices[m, 0] - 1
    v2index = indices[m, 1] - 1
    v3index = indices[m, 2] - 1
    ret.append(v1index)
    ret.append(v2index)
    ret.append(v3index)
    return ret

# F is a vector of size 3M x 3. F[m] is the 3x3 matrix which is the transpose of transformation Sm performed on triangle m.
# Blist is a list of 2x3 matrices (there shall be M of them in total).
# Each matrix in B is abcdef mentioned in page 65 of Doctor Sumner's thesis
# suppose that target transformation = Wj~ * Rj.inverse * Qja.transpose. Then B is precisely Rj.inverse * Qja.transpose.It is
# computed by taking the qr factorization of Wj as described in the thesis.
# Solves Ax = F where A,F are known and x is the unknown.
def Transfer(F, Blist, N, indices):
    M = len(F)/3 # number of corresponding triangles
    print "Constructing A"
    A = scipy.sparse.lil_matrix((3 * M, N), dtype=np.float64)
    print "Allocation done"
    for m in range(0, M):
        indexlist = matrixIndexToXindex(indices, m)
        a = Blist[m][0, 0]
        b = Blist[m][0, 1]
        c = Blist[m][0, 2]
        d = Blist[m][1, 0]
        e = Blist[m][1, 1]
        f = Blist[m][1, 2]
        A[3 * m, indexlist[1]] = a
        A[3 * m, indexlist[2]] = d
        A[3 * m, indexlist[0]] = -a - d
        A[3 * m + 1, indexlist[1]] = b
        A[3 * m + 1, indexlist[2]] = e
        A[3 * m + 1, indexlist[0]] = -b - e
        A[3 * m + 2, indexlist[1]] = c
        A[3 * m + 2, indexlist[2]] = f
        A[3 * m + 2, indexlist[0]] = -c - f
    print "Matrix A constructed. Constructing ATA and ATF", datetime.datetime.now()
    ATA = A.transpose().dot(A)
    ATF = A.transpose().dot(F)
    print "Solving xhead"
    xhead = scipy.sparse.linalg.spsolve(ATA, ATF, use_umfpack=False)
    print "solve done"
    return xhead

def main():
    print "reading source obj file"
    sourceV, M, sourcevertices, sourcetrianglevertices, sourceindices = readObj("example/face-reference.obj")
    print "reading source transform obj file"
    transV, M, transvertices, transtrianglevertices, transindices = readObj("example/face-01-anger.obj")
    print "reading target obj file"
    targetV, M, targetvertices, targettrianglevertices, targetindices = readObj("example/face-05-laugh.obj")

    N = targetV

    print "constructing F"
    F = np.zeros((3*M, 3))
    for m in range(0,M):
        #compute source transformation matrix S = V~ * V-1
        #compute V~, the 3x3 matrix of transformed positions
        Vhead = np.zeros((3, 3))
        Vhead[:, 0] = transtrianglevertices[m, 1,:] - transtrianglevertices[m, 0,:]
        Vhead[:, 1] = transtrianglevertices[m, 2,:] - transtrianglevertices[m, 0,:]
        Vhead[:, 2] = transtrianglevertices[m, 3,:] - transtrianglevertices[m, 0,:]

        #compute V
        V = np.zeros((3, 3))
        V[:, 0] = sourcetrianglevertices[m, 1,:] - sourcetrianglevertices[m, 0,:]
        V[:, 1] = sourcetrianglevertices[m, 2,:] - sourcetrianglevertices[m, 0,:]
        V[:, 2] = sourcetrianglevertices[m, 3,:] - sourcetrianglevertices[m, 0,:]

        #compute inverse of V
        Vinverse = np.linalg.inv(V)

        #compute S
        S = np.dot(Vhead, Vinverse)

        #fill F
        F[3*m:3*m+3,:] = S.transpose()

    print "constructing Blist"
    Blist = []
    for m in range(0,M):
        wj = np.zeros((3, 2))
        wj[:, 0] = targettrianglevertices[m, 1] - targettrianglevertices[m, 0]
        wj[:, 1] = targettrianglevertices[m, 2] - targettrianglevertices[m, 0]
        q, r = np.linalg.qr(wj)
        q = np.array(q)
        r = np.array(r)
        B = np.dot(np.linalg.inv(r), np.transpose(q))
        Blist.append(B)

    #Transfering
    print "Transfering: "
    xhead = Transfer(F, Blist, N, targetindices)

    #save results
    print "Writing to file"
    resultfile = "example/result.obj"
    f = open(resultfile, "w")
    for i in range(0, N):
        f.write("v")
        if xhead[i, 0] < 0:
            f.write("   " + str(round(xhead[i, 0], 9)))
        else:
            f.write("    " + str(round(xhead[i, 0], 9)))
        if xhead[i, 1] < 0:
            f.write("   " + str(round(xhead[i, 1], 9)))
        else:
            f.write("    " + str(round(xhead[i, 1], 9)))
        if xhead[i, 2] < 0:
            f.write("   " + str(round(xhead[i, 2], 9)))
        else:
            f.write("    " + str(round(xhead[i, 2], 9)))
        f.write("\n")
    f.close()

main()



















