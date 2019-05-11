# TDA-of-K-means-on-1st-PL
There are two files with C++ code:

1. Features_MPI.cpp: takes categorical time series as input from Categories.txt, 
obtains the Walsh-Fourier transform (WFT) and then constructs the first-order 
persistence landscape (PL) for each time series. These form the features that 
will be input into the Kmeans procedure. 

2. Kmeans_MPI.cpp: takes PLs as inputs into the K-means algorithm and outputs 
the cluster labels and total within cluster sum of squares. 

Note: Both cpp code are run in parallel on the UConn HPC cluster.
The steps are shown below.

Step 1: First, compile the code by using:
module purge
module load java/1.8.0_162
module load zlib/1.2.11
module load gcc/5.4.0-alt
module load mpi
mpic++ -std=c++11 Features_MPI.cpp -o Features_MPI
mpic++ -std=c++11 Kmeans_MPI.cpp -o Kmeans_MPI

Step 2: Submit these jobs to the cluster using the respective files: 
Features_MPI.sh and Kmeans_MPI.sh 


Input and Output File Description

Input Files:
a) The input file categories.txt must be in the same folder as the code files 
and is read into the main function in Features_MPI.cpp. 
The format of this file is a matrix. 

Row 1 consists of names for the entries of the features. 
Starting from Row 2, the first value is the name of the row and the
rest of the values are for the features. The example below saves two objects 
of size two. The names for the features are V1 and V2 and the names for these 
objects are 1 and 2.
 The file will then be:  
"V1" "V2"
"1" 0 0 
"2" 1 0

b) The input file sample.txt should also be in the same folder as the 
code files and is read into the main function in Features_MPI.cpp. 
It saves the random indexes when assigning time series to different processors of 
the Multi-Processor Interface (MPI). The format of the file is a list of integers 
denoting the indices of the categorical input time series. 

Output Files from Features_MPI.cpp:

(i) WFT1.txt to WFT100.txt: These are the WFTs from all the 100 processors.
 
Row 1 gives the names of the entries of the WFT.
 
Starting from Row 2 are the space separated values for each entry of the WFT.
For example, the file will be:
 
V1 V2

1e-4 0
0    1


(ii) land1.txt to land100.txt: These are the PLs from all 100 processors. 
Row 1 gives the names of the entries of the PLs. Starting from Row 2 are space 
separated values for each entry of the PL.
 For example, the file will be:

V1 V2 

0 0

1 1e-3

(iii) Features_MPI_time.txt: This file contains computation times for constructing 
the WFTs and the PLs.
 Row 1 is the computation time for the WFT in seconds.
 
Row 2 is the computation time for the 1st PL in seconds.


Output Files from Kmeans_MPI.cpp:

(1) tda2.txt to tda5.txt: These are the comma separated cluster labels while 
setting the number of clusters as 2, 3, 4, and 5 respectively. Importantly, 
the order of the cluster labels are based on the orders in sample.txt, which 
are the random indexes when assigning time series to different processors. 


(2) tdawithin2.txt to tdawithin5.txt: These are the total within cluster sums of 
squares for each cluster (group) when setting the number of clusters as 2, 3, 4, 
and 5 respectively.
 When the number of clusters is 2, there will be 2 entries in 
tdawithin2.txt. 

(3) mpitda2time.txt to mpitda5time.txt: These are the computation times for the 
Kmeans algorithm while setting the number of clusters as 2, 3, 4, and 5 respectively.

