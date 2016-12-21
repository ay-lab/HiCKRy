#! /usr/bin/env python
#@author: arya

##import statements
import sys
import numpy as np
import os

try:
    import cPickle as pickle
except:
    import pickle
    print "slower pickle being used"

import scipy.sparse as sps
from spsIO import loadBed, loadNoble 
from matrixOperations import stchMatrix 
import warnings
warnings.filterwarnings('ignore')
import argparse


###
#run using main.py filePath chrLens Resolution chrNum(s)/'whole' outputDir --balance --graphs
parser = argparse.ArgumentParser(description = "Read instructions")
parser.add_argument('filePath', help="Path to the file(s) to be read. If bed/matrix file put bed file first and matrix file second", nargs='+', metavar="FILE")
parser.add_argument('type', help="Type of file to be read. If bed/matrix combo type 'bed', if not type 'fithic'")
parser.add_argument('outputPath', help='Path to directory for outputted files')


#reqd. for fithic format
parser.add_argument('-r', '--resolution', help='Resolution i.e. 40000', type=int)
parser.add_argument('-l' ,'--chrLens', help='Input file for chromosome lengths')
parser.add_argument('-c', '--chrNum', help='Chromosome number to read in following format: chrNum If whole genome reading is required the command becomes: whole')


parser.add_argument('-b','--bias', help='Should we calculate the bias values?', action='store_true')
parser.add_argument('-tr', '--toRemove', help='Percent of sparse row/columns to remove from matrix', type = int)
parser.add_argument('-g','--graphs', help='Turn on output graphs', action='store_true')
parser.add_argument('-o', '--outputDS', help='Output doubly stochastic matrix only.', action = 'store_true')

def alphaStep(s):
    try:
        list = s.split(',')
        list = map(float, list)
    except:
        raise argparse.ArgumentTypeError("Arguments should be ints denoting 'start, stop, number of evenly spaced samples'")
    return list

parser.add_argument('-as','--alphaStep', help="Describe range of alpha and step through two arguments. The first will be alpha and the next will be step. Each are ints denoting 'start, stop, number of evenly spaced numbers within range'", nargs=2, type=alphaStep)
parser.add_argument('-osm', '--outputSmoothedMatrix', help="Option to output a smoothed matrix per the specified range of the alphaStep option.", action='store_true')
parser.add_argument('-f', '--outputFormat', help="Describe how you would like the output files (biasvalues and/or doubly stochastic matrix) to be formatted. If 'fithic' the typical fithic format will be followed, 'bed' yields a .matrix file format similar to bed/matrix files")
args = parser.parse_args()





### PARSE FILE TYPE ###

#determine if we're looking at a bed file
bedFileType = False 
if args.type == 'bed':
    bedFileType = True
    if len(args.filePath) is not 2:
        print "Requires bed file and a matrix file" 
        print args.filePath
        sys.exit(2)

#determine if we're looking at a noble file
nobleFileType = False
if args.type == 'fithic':
    nobleFileType = True
    if len(args.filePath) is not 1:
        print "Only the path to the noble file is required"
        print args.filePath
        sys.exit(2)
    if args.resolution is None:
        print "Default resolution is 40000, if this is not the resolution of the data. Retry with --resolution option" 
        resolution = 40000
    if args.chrLens is None:
        print "Need chrLens to match against for Noble file type. Retry with --chrLens option"
        sys.exit(2)
    
### PARSE FILE PATH ###
if bedFileType:
    if not os.path.exists(args.filePath[0]) or not os.path.exists(args.filePath[1]):
        print "File paths do not exist"
        sys.exit(2)

    bedPath = args.filePath[0]
    matrixPath = args.filePath[1]

if nobleFileType:
    if not os.path.exists(args.filePath[0]):
        print "File path does not exist"
        sys.exit(2)

    noblePath = args.filePath[0]

### PARSE OUTPUT PATH ###
outputPath = args.outputPath


#####NOW PARSING OPTIONAL ARGUMENTS#####

### PARSE RESOLUTION ###
if args.resolution:
    resolution = args.resolution
    if bedFileType:
        print "Resolution option will be ignored... Resolution will be determined by bed file"


### PARSE CHR LENGTHS ###
if args.chrLens:
    import math
    revFragsDic = {}
    allFragsDic = {}
    lenDic = {}
    lens=[]
    lengths = {}
    c = 0
    with open(args.chrLens, 'r') as infile: 
        for line in infile:
            ch,l=line.split()
            lenDic[str(ch)] = int(l)
            if ch not in allFragsDic:
                allFragsDic[ch]={}
                revFragsDic[ch] = {}
            for i in range(int(math.ceil(1.0*int(l)/resolution))):
                mid=int(resolution/2)+i*resolution
                allFragsDic[ch][mid]=c
                revFragsDic[ch][c]=mid
                c+=1
            lens.append(c) 
            lengths[ch] = c
    print lens 
    infile.close()

### PARSE CHRNUM ###
wholeGenome = False
if args.chrNum:
    if args.chrNum == 'whole':
        print "Whole Genome will be read"
        wholeGenome = True
        chrNum = 'whole'
    else:
        chrNum = args.chrNum


### PARSE BIAS VALUES OPTION ###
biasValues = False
if args.bias:
    biasValues = True


### PARSE PERCENT OF SPARSE ROWS/COLUMNS TO REMOVE ###
percentOfSparseToRemove = 0
if args.toRemove:
    if args.toRemove > 100 or args.toRemove < 0:
        print "Illegal value of toRemove option must be between 0-100" 
        sys.exit(2)
    percentOfSparseToRemove = args.toRemove

### PARSE GRAPH OPTION ###
graphOption = False
if args.graphs:
    graphOption = True


### PARSE OUTPUT DOUBLY STOCHASTIC MATRIX FILE OPTION ###
outputNormalizedMatrixFile = False
if args.outputDS:
    outputNormalizedMatrixFile = True


### PARSE RANDOMWALK ALPHASTEP OPTION ###
randomWalkOption = False
if args.alphaStep:
    randomWalkOption = True
    from randomWalks import smoothedContactCounts
    Alpha = args.alphaStep[0]
    Steps = args.alphaStep[1]

### PARSE OUTPUT SMOOTHED MATRIX FILE OPTION ###
outputSmoothedMatrixFile = False
if args.outputSmoothedMatrix:
    outputSmoothedMatrixFile = True


### PARSE OUTPUT FORMAT OPTION ###
fithicOutputFormat = False
bedOutputFormat = False
if args.outputFormat:
    if args.outputFormat == "fithic":
        fithicOutputFormat = True
    elif args.outputFormat == "bedOutputFormat":
        bedOutputFormat = True
    if outputNormalizedMatrixFile is False and biasValues is False:
        print "Nothing to output. -f option is being misused. Either remove it or use --outputDS or -bias option."
        sys.exit(2)

if nobleFileType:
    fithicOutputFormat = True

if bedFileType:
    bedOutputFormat = True

if biasValues and fithicOutputFormat and bedFileType:
    from spsIO import outputBiasFileNobleForBed

if biasValues and fithicOutputFormat and nobleFileType:
    from spsIO import outputBiasFileNobleForNoble

if biasValues and bedOutputFormat and bedFileType:
    from spsIO import outputBiasFileBedForBed

if biasValues and bedOutputFormat and nobleFileType:
    from spsIO import outputBiasFileBedForNoble


#### MAKE DIRECTORIES ####

#define base output folder
base = os.path.basename(args.filePath[0])
outputPath = os.path.join(outputPath, "NormedHiC_" + base)
if not os.path.exists(outputPath):
    try:
        os.makedirs(outputPath)
    except:
        if not os.path.isdir(outputPath):
            raise

#define bias path
if biasValues:
    biasPath = os.path.join(outputPath, "Bias.Values")
    if not os.path.exists(biasPath):
        try:
            os.makedirs(biasPath)
        except:
            if not os.path.isdir(biasPath):
                raise

#define pickle path
picklePath = os.path.join(outputPath, "Pickled.Files")
if not os.path.exists(picklePath):
    try:
        os.makedirs(picklePath)
    except:
        if not os.path.isdir(picklePath):
            raise

#define graph path
graphPath = None
if graphOption:
    graphPath = os.path.join(outputPath, "Graphs")
    if not os.path.exists(graphPath):
        try:
            os.makedirs(graphPath)
        except:
            if not os.path.isdir(graphPath):
                raise

#define matrixFilePath
matrixFilePath = None
if outputNormalizedMatrixFile or outputSmoothedMatrixFile:
    matrixFilePath = os.path.join(outputPath, "Matrix.Files")
    if not os.path.exists(matrixFilePath):
        try:
            os.makedirs(matrixFilePath)
        except:
            if not os.path.isdir(matrixFilePath):
                raise

def main():
    if bedFileType:
        R = loadBed(bedPath, matrixPath, picklePath, graphPath)

    if nobleFileType:
        R = loadNoble(chrNum, resolution, noblePath, picklePath, lenDic, allFragsDic, graphPath)

    stchMatrix(picklePath, percentOfSparseToRemove, graphPath, biasValues, matrixFilePath, outputNormalizedMatrixFile, fithicOutputFormat, bedOutputFormat)
   
    if biasValues:
        if fithicOutputFormat and bedFileType:
            outputBiasFileNobleForBed(picklePath, biasPath, bedPath)

        if fithicOutputFormat and nobleFileType:
            outputBiasFileNobleForNoble(picklePath, lengths, revFragsDic, biasPath, chrNum, resolution)

        if bedOutputFormat and bedFileType:
            outputBiasFileBedForBed(picklePath, biasPath, bedPath)

        if bedOutputFormat and nobleFileType:
           outputBiasFileBedForNoble(picklePath, lengths, revFragsDic, biasPath, chrNum, resolution)


    if randomWalkOption:
        smoothedContactCounts(Alpha, Steps, picklePath, R, graphPath, matrixFilePath, outputSmoothedMatrixFile)


if __name__=='__main__':
    main()
