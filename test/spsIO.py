import os,gzip
import scipy.sparse as sps
import numpy as np
import time
import scipy.io as sio
import csv


def loadBed(bedPath, matrixPath, graphPath=None):
    print "Loading..."
    print "The whole BED file will be read"
    startTime = time.time()

    with open(bedPath, 'r') as bedFile:
        data = np.loadtxt(bedFile, dtype="str")
        lengths = [(data[:, 0] == i).sum() for i in np.unique(data[:, 0])]
        lengths = np.array(lengths)
    bedFile.close()
    n = lengths.sum()

    with open(matrixPath, 'r') as matrixFile:
        data = np.loadtxt(matrixFile, dtype=int)
        ones = np.ones((data.shape[0], 2))
        zero = np.zeros((data.shape[0],1))
        toSub = np.hstack((ones,zero))
        data = data - toSub
        dirty_mtx = sps.coo_matrix((data[:, 2], (data[:, 0], data[:, 1])), shape = (n, n))
    matrixFile.close()


    #convert to csr!
    dirty_mtx = sps.csr_matrix(dirty_mtx)


    #return transpose (csc) and add to dirty_mtx(csr) to make symmetric
    transp = dirty_mtx.transpose()
    dirty_mtx = dirty_mtx + transp
    del(transp)
    R = sps.csr_matrix.sum(dirty_mtx)

    if graphPath is not None:
        from plot import plotMatrix
        plotMatrix(dirty_mtx, graphPath, "Raw.Mtx")

    #done loading!
    endTime = time.time()
    print("Loading took %f" % (endTime - startTime))
    return R, dirty_mtx


def outputBedMatrix(matrix, filename, matrixFilePath):
    saveLocation = os.path.join(matrixFilePath, filename)
    with open(saveLocation, 'w') as matrixFile:
        values = matrix.data
        colIndex = matrix.indices
        rowPtrs = matrix.indptr
        for n in range(matrix.shape[0]):
            matrixFile.write(("%s\t%s\t%.4f\n") % ((rowPtrs[n]+1), (colIndex[n]+1), (values[n])))

#TODO test
def outputNobleMatrixForNoble(matrix, filename, outputFilePath, revFragsDic, lengths, resolution=40000, chrNum = "whole"):
    saveLocation = os.path.join(outputFilePath, filename)
    
    with open(saveLocation, 'w') as matrixFile:
        values = matrix.data
        colIndex = matrix.indices
        rowPtrs = matrix.indptr
        
        if chrNum == 'whole':
            ctrI = 0
            ctrJ = 0
            for i in range(matrix.shape[0]):
                noLowChLen = True
                ch1 = "chrNULL" 
                for chromosomes in lengths:
                    currChLen = lengths[chromosomes]
                    if (ctrI < currChLen and noLowChLen is True):
                        lowChLen = currChLen
                        noLowChLen = False
                        ch1 = chromosomes
                   
                    if (ctrI < currChLen and currChLen < lowChLen):
                        lowChLen = currChLen 
                        ch1 = chromosomes
                mid1 = revFragDic[ch1][ctrI]
                ctrI+=1
                
                for j in range(matrix.shape[1]):
                    noLowChLen = True
                    ch1 = "chrNULL" 
                    for chromosomes in lengths:
                        currChLen = lengths[chromosomes]
                        if (ctrJ < currChLen and noLowChLen is True):
                            lowChLen = currChLen
                            noLowChLen = False
                            ch2 = chromosomes
                       
                        if (ctrJ < currChLen and currChLen < lowChLen):
                            lowChLen = currChLen 
                            ch2 = chromosomes
                    mid2 = revFragDic[ch2][ctr]
                    ctrJ+=1
                     
                    matrixFile.write(str(ch1)+"\t"+str(mid1)+"\t"+str(matrix[i][j])+"\t"+str(ch2)+"\t"+str(mid2)+"\t"+str(matrix[i][j])+"\n")

        else:
            ctrI = 0
            ctrJ = 0
            for i in range(matrix.shape[0]):
                ch1 = chrNum 
                mid1 = revFragDic[ch1][ctrI]
                ctrI+=1
                
                for j in range(matrix.shape[1]):
                    ch2 = chrNum
                    mid2 = revFragDic[ch2][ctr]
                    ctrJ+=1
                     
                    matrixFile.write(str(ch1)+"\t"+str(mid1)+"\t"+str(matrix[i][j])+"\t"+str(ch2)+"\t"+str(mid2)+"\t"+str(matrix[i][j])+"\n")
    
    matrixFile.close()
             

#TODO Test
def outputNobleMatrixForBed(matrix, filename, outputFilePath, bedFilePath):
    saveLocation = os.path.join(outputFilePath, filename)
    res = 0
    fragDic = {}
    with open(bedFilePath, 'r') as bedFile:
        for lines in bedFile:
            line = lines.rstrip().split()
            chrNum = line[0]
            start = line[1]
            en = line[2]
            if res == 0: res = int(en)-int(start)
            mid = int(start)+int(res/2)
            index = int(line[3]) - 1
            fragDic[index] = [chrNum, mid]
   

    lineCount = 0
    with open(saveLocation, 'w') as matrixFile:
        values = matrix.data
        colIndex = matrix.indices
        rowPtrs = matrix.indptr
        for n in range(matrix.shape[0]):
            i = rowPtrs[n]
            j = colIndex[n]
            k = values[n]
            matrixFile.write(str(fragDic[i][0])+"\t"+str(fragDic[i][1])+"\t"+str(fragDic[j][0])+"\t"+str(fragDic[j][1])+"\t"+"%.4f\n" % (k))
            lineCount += 1
            if lineCount%1000000==0: print("%d million lines have been written" % int(lineCount/1000000))



#loadNoble converts a raw data file into a sparse matrix to be computed
##chrNum takes the format of "chrN" where N is the desired chromosome 
#  to load
def loadNoble(chrNum, resolution, noblePath, lenDic, fragDic, graphPath=None):
    import math
    hiCFile = gzip.open(noblePath, 'r')
    print "Loading..."  
    startTime = time.time()
    n = 0
    halfRes = resolution/2

    if chrNum == 'whole':
        for key in lenDic:
            n += int(math.ceil(1.0*int(lenDic[key])/resolution))

    else:
        n = int(math.ceil(1.0*int(lenDic[chrNum])/resolution))
    
    #construct a sparse matrix of max resolution
    dirty_mtx = sps.lil_matrix((n,n), dtype = np.int64) 

    #load values into the sparse matrix constructed earlier! 
    for line in hiCFile:

        if chrNum != "whole" and line.startswith(chrNum):
            fileLine = line.rstrip().split()
            i = (int(fileLine[1])-halfRes)/resolution
            j = (int(fileLine[3])-halfRes)/resolution
            k = float(fileLine[4])


            if (fileLine[0] == fileLine[2]):
                try:
                    dirty_mtx[i,j] = k
                    dirty_mtx[j,i] = k

                except:
      #              print fileLine
                    continue

        else:
            fileLine=line.rstrip().split()
            firstCh = fileLine[0]
            secCh = fileLine[2]

            mid1 = int(fileLine[1])
            mid2 = int(fileLine[3])
            k = float(fileLine[4])

            i = fragDic[firstCh][mid1]
            j = fragDic[secCh][mid2]
            dirty_mtx[i,j] = k
            dirty_mtx[j,i] = k

    if graphPath is not None:
        from plot import plotMatrix 
        plotMatrix(dirty_mtx, graphPath, "Raw.Mtx")

    #convert to csr!
    dirty_mtx = dirty_mtx.tocsr()
    R = sps.csr_matrix.sum(dirty_mtx)

    #done loading!
    endTime = time.time()
    print("Loading took %f" % (endTime - startTime))

    hiCFile.close()
    return R, dirty_mtx

def outputMatrixFile(matrix, fileName, matrixFilePath):
    fname = os.path.join(matrixFilePath, fileName) 
    sio.mmwrite(fname, matrix, field='real')

def outputBiasFileNobleForNoble(biasCol, lengths, revFragDic, outputFilePath, chrNum='whole', resolution=40000):

    bpath=os.path.join(outputFilePath, "Bias.Values")
    
    with open(bpath, 'w') as f:
        
        if chrNum == 'whole':
            ctr = 0
            for values in np.nditer(biasCol):
                noLowChLen = True
                ch = "chrNULL" 
                for chromosomes in lengths:
                    currChLen = lengths[chromosomes]
                    if (ctr < currChLen and noLowChLen is True):
                        lowChLen = currChLen
                        noLowChLen = False
                        ch = chromosomes
                   
                    if (ctr < currChLen and currChLen < lowChLen):
                        lowChLen = currChLen 
                        ch = chromosomes
                    

                mid = revFragDic[ch][ctr]
                ctr+=1
                f.write(ch + "\t" + str(mid) + "\t" + str(values) + "\n")
        
        else:
            halfRes = resolution/2
            ctr = 0
            for values in np.nditer(biasCol):
                mid = ctr*resolution + halfRes
                f.write(chrNum + "\t" + str(mid) + "\t" + str(values) + "\n") 
                ctr += 1

    f.close()


def outputBiasFileNobleForBed(biasCol, outputFilePath, bedFilePath):


    bpath=os.path.join(outputFilePath, "Bias.Values")
    
    with open(bpath, 'w') as f:
    
        with open(bedFilePath, 'r') as bedFile:
            lines = bedFile.readlines()
            #print lines
            ctr = 0
            for values in np.nditer(biasCol):        
                fileLine = lines[ctr].rstrip().split() 
                chrNum = fileLine[0]
                mid = ((int(fileLine[1])+int(fileLine[2]))/2)
                f.write(chrNum + "\t" + str(mid) + "\t" + str(values) + "\n")
                ctr += 1
        bedFile.close()
    f.close()


def outputBiasFileBedForNoble(biasCol, lengths, revFragDic, outputFilePath, chrNum='whole', resolution=40000):
    bpath=os.path.join(outputFilePath, "Bias.Values")
    
    with open(bpath, 'w') as f:
        
        if chrNum == 'whole':
            ctr = 0
            for values in np.nditer(biasCol):
                noLowChLen = True
                ch = "chrNULL" 
                for chromosomes in lengths:
                    currChLen = lengths[chromosomes]
                    if (ctr < currChLen and noLowChLen is True):
                        lowChLen = currChLen
                        noLowChLen = False
                        ch = chromosomes
                   
                    if (ctr < currChLen and currChLen < lowChLen):
                        lowChLen = currChLen 
                        ch = chromosomes
                    

                mid = revFragDic[ch][ctr]
                ctr+=1
                f.write(ch + "\t" + str(mid) + "\t" + str(ctr+1) + "\t" + str(values) + "\n")
        
        else:
            halfRes = resolution/2
            ctr = 0
            for values in np.nditer(biasCol):
                mid = ctr*resolution + halfRes
                f.write(chrNum + "\t" + str(mid) + "\t" + str(ctr+1) + "\t" + str(values) + "\n") 
                ctr += 1

    f.close()



def outputBiasFileBedForBed(biasCol, outputFilePath, bedFilePath):
    bpath=os.path.join(outputFilePath, "Bias.Values")
    
    with open(bpath, 'w') as f:
        with open(bedFilePath, 'r') as bedFile:
            lines = bedFile.readlines()
            print lines
            ctr = 0
            for values in np.nditer(biasCol):        
                fileLine = lines[ctr].rstrip().split() 
                chrNum = fileLine[0]
                mid = ((int(fileLine[1])+int(fileLine[2]))/2)
                f.write(chrNum + "\t" + str(mid) + "\t" + str(ctr+1) + "\t" + str(values) + "\n")
                ctr += 1
        bedFile.close()
    f.close()

