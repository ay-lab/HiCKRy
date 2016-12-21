# HiCKRy

## What is HiCKRy?
A Python based tool developed by the Ay Lab to normalize Hi-C data using the Knight-Ruiz algorithm for matrix normalization. 

## Why use HiCKRy?
HiCKRy sets itself apart from other Hi-C normalization tools by being:
- Flexible
- Efficient
- Fast

The field of chromatin conformation capture technology has been a testimony to the impact high-throughput sequencing has had on bioinformatics. The rise of high-throughput sequencing has allowed scientists to go from analyzing the chromosome conformation at a single locus (3C), to a set of loci (5C, ChIA-PET), to the whole genome (Hi-C).

Unfortunately, this plethora of data also presents the computational challenge of dealing with it. As Hi-C data continues to increase in both size and resolution, developing tools that can handle this 'data deluge' is necessary for the continued analysis of Hi-C data and all the wonderful insights it may afford.

### How do I use HiCKRy?
First make sure that all the necessary dependencies are installed.
- Python 2.7
- Numpy/Scipy
- Matplotlib - only if using the graphing option

Using the following command details the various options of HiCKRy.
```
main.py -h
```

