import os
import numpy as np
import scipy.sparse as sps
from knightRuiz import knightRuizAlg

def stchMatrix(picklePath, percentOfSparseToRemove, graphPath=None, biasValues=False, matrixFilePath=None, outputNormalizedMatrixFile=False, fithicOutputType=False, bedOutputType=False):
    import cPickle as pickle


    fname = os.path.join(picklePath, "Raw.Mtx")
    with open(fname, 'rb') as f:
        rawMatrix = pickle.load(f)
    f.close()

    R = rawMatrix.sum()

    mtxAndRemoved = removeZeroDiagonalCSR(rawMatrix, percentOfSparseToRemove)
    initialSize = rawMatrix.shape[0]
    rawMatrix = mtxAndRemoved[0]
    removed = mtxAndRemoved[1]
    newSize = rawMatrix.shape[0]


    print "Generating Normalized Matrix"
    result = knightRuizAlg(rawMatrix)
    colVec = result[0]
    if np.isnan(np.sum(colVec)):
        print "Too few rows/columns removed... try again"
        return None
    x = sps.diags(colVec.flatten(), 0, format='csr')
    
    if biasValues:
        bias = computeBiasVector(colVec)
        biasWZeros = addZeroBiases(removed, bias)
        
        biasFileName = os.path.join(picklePath, "Bias.Values")
        with open(biasFileName, 'wb') as f:
            pickle.dump(biasWZeros, f)
        f.close()

    del(colVec)

    normalizedMatrix = x.dot(rawMatrix.dot(x))
    n = normalizedMatrix.shape[0]
    print "Normalized Matrix Generated"
    
    #calculate normalization difference
    difference = abs(rawMatrix.shape[0] - sps.csr_matrix.sum(normalizedMatrix))

    scalar = R/n

    if graphPath is not None:
        from plot import plotMatrix
        print "Generating heatmap for Normalized Matrix"
        plotMatrix((normalizedMatrix*scalar), graphPath, "Normalized.Mtx")



    if outputNormalizedMatrixFile:
        from spsIO import outputMatrixFile
        print "Outputting Normalized Matrix"
#        if fithicOutputType:
#            from spsIO import outputNobleMatrix 
#            outputNobleMatrix(normalizedMatrix, "Normalized.mtx", matrixFilePath)
        if bedOutputType or fithicOutputType:
            from spsIO import outputBedMatrix
            outputBedMatrix((normalizedMatrix*scalar), "Normalized.matrix", matrixFilePath)

    fileName = os.path.join(picklePath, "Normalized.Mtx")
    with open(fileName, 'wb') as f:
        pickle.dump(normalizedMatrix, f)
    f.close()

    fileName = os.path.join(picklePath, "removed")
    with open(fileName, 'wb') as f:
        pickle.dump(removed, f)
    f.close()
    print "Normalized Matrix pickled"

def computeBiasVector(x):
#    print x
    one = np.ones((x.shape[0],1))
#    print one
    x = one/x
#    print x
    sums = np.sum(x)
#    print sums
    avg = (1.0*sums)/x.shape[0]
#    print avg
    bias = np.divide(x,avg)
    return bias

def addZeroBiases(lst, vctr):
    for values in lst:
        vctr = np.insert(vctr,values,-1,axis=0)        
    return vctr

def dropcols_coo(M, idx_to_drop):
    idx_to_drop = np.unique(idx_to_drop)
    C = M.tocoo()
    keep = ~np.in1d(C.col, idx_to_drop)
    C.data, C.row, C.col = C.data[keep], C.row[keep], C.col[keep]
    C.col -= idx_to_drop.searchsorted(C.col) # decrement column indices
    C._shape = (C.shape[0], C.shape[1] - len(idx_to_drop))
    return C.tocsr()

def removeRowCSR(mat, i):
    if not isinstance(mat, sps.csr_matrix):
        raise ValueError("works only for CSR format -- use .tocsr() first")
    n = mat.indptr[i+1] - mat.indptr[i]
    if n > 0:
        mat.data[mat.indptr[i]:-n] = mat.data[mat.indptr[i+1]:]
        mat.data = mat.data[:-n]
        mat.indices[mat.indptr[i]:-n] = mat.indices[mat.indptr[i+1]:]
        mat.indices = mat.indices[:-n]
    mat.indptr[i:-1] = mat.indptr[i+1:]
    mat.indptr[i:] -= n
    mat.indptr = mat.indptr[:-1]
    mat._shape = (mat._shape[0]-1, mat._shape[1])

def removeZeroDiagonalCSR(mtx, i=0, toRemovePre=None):
    iteration = 0
    toRemove = []
    ctr = 0
    
    if toRemovePre is not None:
        for items in toRemovePre:
            toRemove.append(items)   

    if i == 0:
        diagonal = mtx.diagonal()
#        print diagonal
        for values in diagonal:
            if values == 0:
                toRemove.append(ctr)
            ctr += 1
    
    else:
        rowSums = mtx.sum(axis=0)
        rowSums = list(np.array(rowSums).reshape(-1,)) 
        rowSums = list(enumerate(rowSums))
        for value in rowSums:
            if int(value[1]) == 0:
                toRemove.append(value[0])
                rowSums.remove(value) 
        rowSums.sort(key=lambda tup: tup[1])
        size = len(rowSums)
        perc = i/100.0
        rem = int(perc * size)
        while ctr < rem:
            toRemove.append(rowSums[ctr][0])
            ctr += 1
    list(set(toRemove))
    toRemove.sort()
    #print toRemove
    mtx = dropcols_coo(mtx, toRemove)
    for num in toRemove:
        if iteration != 0:
            num -= iteration
        removeRowCSR(mtx,num)
        iteration +=1
    return [mtx, toRemove]

def addZeroes(mtx, toAdd):
    indices = np.array(toAdd)
    i = indices - np.arange(len(indices))
    mtx = np.insert(np.insert(mtx,i,0,axis=1),i,0,axis=0)
    return mtx

def main():
    bias = np.array([0, 4, 5])
    print bias

    diag = sps.diags(bias, 0, format='lil')
    print diag.toarray()

if __name__=="__main__":
    main()
