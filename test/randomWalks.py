import numpy as np
import os
import scipy.sparse as sps


def randomWalk(alpha, steps, normalizedMatrix):
    steps = int(steps)

    n = normalizedMatrix.shape[0]
    I = np.identity(n)
    walk = np.linalg.matrix_power(alpha*I + (1-alpha)*normalizedMatrix, steps)

    walk = sps.csr_matrix(walk)
    return [walk, n]

def smoothedContactCounts(Alpha, Steps, normalizedMtx, R, graphPath=None, matrixFilePath = None, outputSmoothedMatrixFile=False):
    for alpha in np.linspace(Alpha[0], Alpha[1], Alpha[2]).tolist():
        for steps in np.linspace(Steps[0], Steps[1], Steps[2]).tolist():
            rw = randomWalk(alpha, steps, normalizedMtx)
            walk = rw[0]
            n = rw[1]
            scalar = R/n
            walk *= scalar
            name = "SmoothedCC.Alpha=" + str(alpha) + ".t=" + str(steps)
            if graphPath is not None:
                from plot import plotMatrix
                plotMatrix(walk, graphPath, name)



            if outputSmoothedMatrixFile:
                from spsIO import outputBedMatrix 
                outputBedMatrix(walk, name, matrixFilePath)
            

def main():
    pass

if __name__=="__main__":
    main()
