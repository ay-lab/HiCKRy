import scipy.sparse as sps
import numpy as np
import sys
import resource

##FUNCTION DESCRIPTION
# knighRuizAlg is an implementation of the matrix balancing algorithm
#  developed by Knight and Ruiz. The goal is to take a matrix A and
#  find a vector x such that, diag(x)*A*diag(x) returns a doubly
#  stochastic matrix

##PARAMETERS
#A is a given numpy array
#tol is error tolerance
#f1 boolean indicating if the intermediate convergance statistics
# should also be outputted
def knightRuizAlg(A, tol=1e-6, f1 = False):
    n = A.shape[0]
    e = np.ones((n,1), dtype = np.float64)
    res = []


    Delta = 3
    delta = 0.1
    x0 = np.copy(e)
    g = 0.9

    etamax = eta = 0.1
    stop_tol = tol*0.5
    x = np.copy(x0)

    rt = tol**2.0
    v = x * (A.dot(x))
    rk = 1.0 - v
#    rho_km1 = np.dot(rk.T, rk)[0, 0]
    rho_km1 = ((rk.transpose()).dot(rk))[0,0]
    rho_km2 = rho_km1
    rout = rold = rho_km1
    
    MVP = 0 #we'll count matrix vector products
    i = 0 #outer iteration count

    if f1:
        print ("it        in. it      res\n"),

    while rout > rt: #outer iteration
        i += 1

        if i > 30:
            break


        k = 0
        y = np.copy(e)
        innertol = max(eta ** 2.0 * rout, rt)
        
        while rho_km1 > innertol: #inner iteration by CG
            k += 1
            if k == 1:
                Z = rk / v
                p = np.copy(Z)
                #rho_km1 = np.dot(rk.T, Z)
                rho_km1 = (rk.transpose()).dot(Z)
            else:
                beta = rho_km1 / rho_km2
                p = Z + beta * p

            if k > 10:
                break





            #update search direction efficiently
            w = x * A.dot(x * p) + v * p
           # alpha = rho_km1 / np.dot(p.T, w)[0,0]
            alpha = rho_km1 / (((p.transpose()).dot(w))[0,0])
            ap = alpha * p
            #test distance to boundary of cone
            ynew = y + ap
            
            if np.amin(ynew) <= delta:
                
                if delta == 0:
                    break

                ind = np.where(ap < 0.0)[0]
                gamma = np.amin((delta - y[ind]) / ap[ind])
                y += gamma * ap
                break

            if np.amax(ynew) >= Delta:
                ind = np.where(ynew > Delta)[0]
                gamma = np.amin((Delta - y[ind]) / ap[ind])
                y += gamma * ap
                break

            y = np.copy(ynew)
            rk -= alpha * w
            rho_km2 = rho_km1
            Z = rk / v
            #rho_km1 = np.dot(rk.T, Z)[0,0]
            rho_km1 = ((rk.transpose()).dot(Z))[0,0]
        x *= y
        v = x * (A.dot(x))
        rk = 1.0 - v
        #rho_km1 = np.dot(rk.T, rk)[0,0]
        rho_km1 = ((rk.transpose()).dot(rk))[0,0]
        rout = rho_km1
        MVP += k + 1
        
        #update inner iteration stopping criterion
        rat = rout/rold
        rold = rout
        res_norm = rout ** 0.5
        eta_o = eta
        eta = g * rat
        if g * eta_o ** 2.0 > 0.1:
            eta = max(eta, g * eta_o ** 2.0)
        eta = max(min(eta, etamax), stop_tol / res_norm)
        if f1:
            print ("%03i %06i %03.3f %e %e \n") % \
                (i, k, res_norm, rt, rout), 
            res.append(res_norm)
    if f1:
        print ("Matrix - vector products = %06i\n") % \
            (MVP),
    
    #X = np.diag(x[:,0])   
    #x = X.dot(A.dot(X))
    return [x,(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1000)]

def main():
    test = sps.rand(5,5,density=0.5,format='csr')
    print test.toarray()
    CC = ((test.sum())/(test.shape[0]*2))
    CCother = test.sum()/test.size 
    result = knightRuizAlg(test)
    col = result[0]
    x = sps.diags(col.flatten(), 0, format='csr')
    mtx = x.dot(test.dot(x))


    print mtx.toarray()
    CCmtx = CC * mtx
    CCothermtx = CCother * mtx
    print CCmtx.toarray()
    print()
    print CCothermtx.toarray()


if __name__=="__main__":
    main()
