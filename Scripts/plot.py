import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import os
import scipy.sparse as sps
import numpy as np

def plotMatrix(matrix, path, name):
    if sps.issparse(matrix):
        if not sps.isspmatrix_coo(matrix): 
            matrix = sps.coo_matrix(matrix)

        x = ((matrix.shape[0])*(matrix.shape[0]))
        sparse = x - matrix.data.shape[0]
        zeroes = np.zeros((sparse,))
        big = sps.hstack((zeroes, matrix.data))

        #print big

        #percentile = np.percentile(big.toarray(), 95) 
        percentile =  np.percentile(sparse, 95) 

        if percentile == 0:
            percentile =  np.percentile(sparse, 95) 

        cm = plt.cm.get_cmap('afmhot_r')
        fig = plt.figure()
        ax = fig.add_subplot(111, axisbg = 'white')  

        plot = ax.scatter(matrix.col, matrix.row, c=matrix.data, s=0.3, edgecolor='', vmax = percentile, cmap = cm) 
            
        cb = plt.colorbar(plot)
        cb.set_label("Contact Counts")
        ax.set_aspect('equal')
        plot = plt.gca().invert_yaxis()
        
        ax.set_xlim((-0.5, (matrix.shape[0]) - 0.5))
        ax.set_ylim((-0.5, (matrix.shape[0]) - 0.5))
        filename = (name + '.png')
        fullpath = os.path.join(path,filename)
        fig.savefig(fullpath,dpi=300)

    else:
        print "wefucked"
        vmaxLim = np.percentile(matrix, 95)
        fig, ax = plt.subplots()
        m = ax.matshow(counts, origin="bottom",cmap="afmhot_r", vmax=vmaxLim) 


        ax.axhline(-0.5, color="#000000", linewidth=1, linestyle="--")
        ax.axvline(-0.5, color="#000000", linewidth=1, linestyle="--")


        cb = fig.colorbar(m)
        cb.set_label("Contact counts")

        ax.set_xlim((-0.5, len(counts) - 0.5))
        ax.set_ylim((-0.5, len(counts) - 0.5))
        filename = (name + '.png')
        fullpath = os.path.join(path,filename)
        fig.savefig(fullpath,dpi=300)
                         


def plotTime(path):
    percentRemoved = []
    timeForKR = []
    outIterKR = []
    inIterKR = []

    timePath = os.path.join(path, "timeData")
    with open(timePath, 'r') as file:
        for line in file:
            line = line.rstrip().split()
            percentRemoved.append(line[0])
            timeForKR.append(line[2])
            outIterKR.append(line[3])
            inIterKR.append(line[4])
    file.close()
   
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax3 = ax1.twinx()
    axes = [ax1, ax2, ax3]
    
    fig.subplots_adjust(right=0.75)

    axes[-1].spines['right'].set_position(('axes', 1.2))
    axes[-1].set_frame_on(True)
    axes[-1].patch.set_visible(False)



    ax1.plot(percentRemoved, timeForKR, 'ro')
    ax1.set_xlabel("Percent Removed")
    ax1.set_ylabel("Time for KR", color = 'r')
    ax1.tick_params(axis='y', color = 'r')

    ax2.plot(percentRemoved, outIterKR, 'bo')
    ax2.set_ylabel("Outer Iterations for KR", color = 'b')
    ax2.tick_params(axis='y', color = 'b')


    ax3.plot(percentRemoved, inIterKR, 'go')
    ax3.set_ylabel("Inner Iterations for KR", color = 'g')
    ax3.tick_params(axis='y', color = 'g')    

    full = os.path.join(path, "timeData.png")
    plt.savefig(full, dpi=300)



def main():
    print "nothing here"

if __name__ == "__main__":
    main()
