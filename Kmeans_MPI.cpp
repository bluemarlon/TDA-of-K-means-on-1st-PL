#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits> 
#include <random>
#include <mpi.h>
#include <cstdlib> 
#include <algorithm>    
#include <time.h> 

using namespace std;

vector<string> Readline(string line){
    /* reads one line of the file  */
    stringstream   lineStream(line);
    vector<string> row = {};
    string         cell;
    while(getline(lineStream, cell, ',')){
        row.push_back(cell);
    }
    return row;
}

void Readtxt(
    string sstm,
    vector<double>& X, /* all 1st landscapes of one processor */
    int& dim, /* the dimension of the observations */
    int& n /* the number of objects/time series */
    ) {
        /* it reads the txt file and save into X */
        ifstream mfile (sstm);
        string line;
        if(mfile.is_open()){
            getline(mfile, line);       
            vector<string> row = Readline(line);
            dim  = row.size();
            while(getline(mfile, line))
            {
                row = Readline(line);
                n++;
                for(int i=0; i < dim; i++) {
                    X.push_back(stod(row[i].c_str()));
                }
            }
            mfile.close();
        }
    }

double EuclideanDistance(
    int dim, /* dim of the 1st PL (100) */
    int label, /* the cluster label */
    int obs, /* the index of the 1st landscapes */
    vector<double> centroid, /* the vector saving all centroid */
    vector<double>& X /* all 1st landscapes of one processor */
    ){
    /* compute the Euclidean distance for one 1st landscapes to the corresponding cluster label */
    double S=0;
    for(int i=0;i< dim;i++){
        S+=pow( centroid[i+dim*label]-X[i+dim*obs],2);
    }
    return pow(S, 0.5);
}

vector<double> FindColMinMax(
    int dim, /* dim of the 1st PL (100) */
    int n, /* the number of time series */
    vector<double>& X /* all 1st landscapes of one processor */
    ){
    vector<double> minima(dim,numeric_limits<double>::max()), maxima(dim, numeric_limits<double>::min());
    for(int i=0;i< n; i++){
        for(int j=0;j<dim;j++){
            if(X[i*dim+j] < minima[j]){ minima[j] = X[i*dim+j]; }
            if(X[i*dim+j] > maxima[j]){ maxima[j] = X[i*dim+j]; }
        }
    }
    minima.insert(minima.end(), maxima.begin(), maxima.end());
    return minima;
}

vector<double> InitializeMeans(
    int dim, /* dim of the 1st PL (100) */
    vector<double>& X, /* all 1st landscapes of one processor */
    int k, // clusters
    vector<double> range
    ){
    /*Initialize means to random numbers between
 *  *     the min and max of each column/feature*/
    vector<double> means(dim*k, 0);
    for(int i=0; i<k; i++){
        for(int j=0; j< dim; j++){
            /*Set value to a random float (adding +-1 to avoid a wide placement of a mean) */
            default_random_engine generator;
            uniform_real_distribution<double> distribution(range[j],range[dim+j]);
            means[i*dim+j] = distribution(generator);
        }
    }
    return means;
}

int Classify(
    int dim, /* dim of the 1st PL (100) */
    int k, /* # clusters */
    int obs, /* index of observations */
    vector<double> centroid, /* centroid */
    vector<double>& X /* all 1st landscapes of one processor */
    )
{
    /*find the obs-th object's cluster label */
    double imax = numeric_limits<double>::max();
    int index = -1;
    for(int i=0; i< k;i++){
        /*Find distance from item to mean*/
        double dis = EuclideanDistance(dim, i, obs, centroid, X);
        if(dis < imax){
            imax = dis;
            index = i;
        }
    }
    return index;
}

vector<int> FindClusters(
        int n, /* the number of time series */
        int k,/* # clusters */
        int dim, /* dim of the 1st PL (100) */
        vector<double> centroid,
        vector<double>& X // all 1st landscapes of one processor
    ){
    /* find the corresponding cluster label for 1st PL*/
    vector<int> clusters;
    for(int i=0;i<n;i++){
        /* Classify item into a cluster*/
        int index = Classify(dim, k, i, centroid, X);
        /* Add item to cluster*/
        clusters.push_back(index);
    }
    return clusters;
}

vector<double> CalculateMeans(
    int n, /* the number of time series */
    int k,/* # clusters */
    int dim, /* dim of the 1st PL (100) */
    vector<double> means, /*centoirds*/
    vector<double>& X, /*the dataset*/
    int maxIterations=1000,
    double threshold = 0.05
    )
{ /* K-means algorithm */

    /*Initialize clusters, the array to hold the number of items in a class*/
    vector<double> clustersize(k, 0);
    /*An array to hold the cluster an item is in*/
    vector<short> labels(n, 0);
    /*calculate means*/
    for(int obs=0;obs<n;obs++){
            /*Classify item into a cluster and update the corresponding means*/
            int index = Classify(dim, k, obs, means, X);
            clustersize[index]++;
            labels[obs] = index;
    }
    means.resize(0);//initialize means all 0;
    means.resize(dim*k);
    /*update the centroids the first time*/
    for(int obs=0;obs<n;obs++){
        for(int ii=0;ii<dim; ii++){
            means[labels[obs]*dim+ii] += X[obs*dim+ii]/clustersize[labels[obs]];
        }
    }
    /*use a flag to exit output*/
    int flagout=0;
    bool flagin = true;
    for(int i=0; i< maxIterations; i++){
        for(int obs=0;obs<n;obs++){
            /*Classify item into a cluster and update the corresponding means*/
            int index = Classify(dim, k, obs, means, X);
            /*Item changed cluster*/
            if(index != labels[obs]){
                flagin = false;
                if(i==0)
                    flagout++;
		        /* UpdateMean(clusterSizes[index],means[index],item);*/
                for(int ii=0;ii<dim;ii++){
                    /*update per observations*/
                    /*update the old cluster index mean*/
                    means[labels[obs]*dim + ii] -= (X[obs*dim + ii] - means[labels[obs]*dim + ii]) / (clustersize[labels[obs]]-1);
                    /*update the new cluster index mean*/
                    means[index*dim + ii] += (X[obs*dim + ii] - means[index*dim + ii]) / (clustersize[index]+1);
                }
                clustersize[labels[obs]]--;
                clustersize[index]++;
                labels[obs] = index;
            }
        }
        /*Nothing changed, return*/
        if(flagin){ /*excending K-means loop*/
           break;
        }
        flagin = true;
    }
    if((double)flagout/n <= threshold) {means.insert(means.end(), 1);
	}else{
        means.insert(means.end(), 0);
    }
    return means;
}
vector<double> Calculatewithin(
        int n, /* the number of time series */
        int k,/* # clusters */
        int dim, /* dim of the 1st PL (100) */
        vector<double> centroid, /* all centroids */
        vector<int> labels, /* cluster labels */
        vector<double>& X /*the dataset*/
    ){ /* compute the total within sum of square */
        vector<double> totwithin(k, 0);
        for(int i=0; i<n; i++){
            totwithin[labels[i]] += pow(EuclideanDistance(dim, labels[i], i, centroid, X), 2);
        }
        return totwithin;
}
int main(
	int argc, char * argv[]
    )
{  
    /*takes categorical time series as input from Categories.txt, 
      obtains the Walsh-Fourier transform (WFT) and then constructs the first-order 
      persistence landscape (PL) for each time series. These form the features that 
      will be input into the Kmeans procedure. 
    */
    /*
    i-th processor takes the corrsponding land(i).txt;
    processor 0 write out the within sum of square of each cluster group in tda(k)totwithin.txt
    and the cluster labels in tda(k).txt, 'k' for the number of cluster in the K-means
    */
    int world_rank, world_size;
    string kint = argv[1];
    string readfile = argv[2];
    string writefile = argv[3];
    /*Initialize the parallel section*/
    MPI_Init( &argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    int k=stoi(kint); /*cluster k*/
    int n=0; /*# of observations*/
    string sstmrd; /*data of reading the file*/
    string sstmwt; /*data of writing the file*/
    /*initialize X, the dataset*/
    vector<double> X;
    sstmrd += "../clust/data/" + readfile + to_string(world_rank+1) + ".txt";
    int dim;
    Readtxt(sstmrd, X, dim, n);
    clock_t begin_time = clock();
    /*Find the minima and maxima for each entry of the 1st landscapes*/
    vector<double> range = FindColMinMax(dim, n, X);

    /*Initialize means at random points*/
    vector<double> means = InitializeMeans(dim, X,k,range);
    /* K-means */
    vector<double> centroid = CalculateMeans(n,k,dim,means,X);
    int flag = 0;//flag=1 for true, no changes
    centroid.resize(centroid.size()-1);
    MPI_Barrier(MPI_COMM_WORLD);
    int outer=0;
    for(; outer<100;outer++){
    /* state two  */
    // master processor update the centroid and send the centroid back to other processors
        double *val = nullptr;
        if(world_rank==0) //only master receive
        {
            val = (double*) calloc(dim*k*world_size, sizeof(double));
        }
        MPI_Gather(&centroid[0],dim*k, MPI_DOUBLE, val, dim*k, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
        if(world_rank==0) //master node update centroid
        {
            vector<double> centroids(val, val + dim*k*world_size);
            vector<double> centroid1 = CalculateMeans(k*world_size, k, dim, centroid, centroids);
            centroid1.resize(centroid1.size()-1);
            centroid.resize(0);
            centroid.insert(centroid.end(), centroid1.begin(), centroid1.end());
        }
        free(val);
        MPI_Bcast(&centroid[0], dim*k, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        /*state three*/
        //each node do K-means;
        vector<double> centroid1 = CalculateMeans(n,k,dim,centroid,X); 
        flag = centroid1[centroid1.size()-1];
        centroid1.resize(centroid1.size()-1);
        centroid = centroid1;
        int flag1;
        MPI_Allreduce(&flag, &flag1, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        if(flag1 == 1) break;
    }
    vector<int> labels = FindClusters(n,k, dim, centroid, X);
    vector<double> totwithin = Calculatewithin(n, k, dim, centroid, labels, X);
    /*state final; make all the output into one processor and then master-write */
    int *val = nullptr;
    if(world_rank>0){
        MPI_Gather(&n, 1, MPI_INT, nullptr, 0, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gatherv(&labels[0], n, MPI_INT, nullptr, 0, nullptr,MPI_INT, 0, MPI_COMM_WORLD);
    }
    if(world_rank==0){
        vector<int> n1(world_size);
        MPI_Gather(&n, 1, MPI_INT, &n1[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
        int valsum = 0;//record cusum of obs on each processors
        vector<int> displs(world_size); //locations of starting point on each processor
        for(int i=0; i < n1.size(); i++){
            displs[i] += valsum;
            valsum += n1[i];   
        }
        val = (int*) calloc(valsum, sizeof(int));
        MPI_Gatherv(&labels[0], n, MPI_INT, val, &n1[0], &displs[0], MPI_INT, 0, MPI_COMM_WORLD);
        labels.resize(0);
        labels.insert(labels.end(), val, val + valsum);
    }
    free(val);
    vector<double> totgobal(k,0);
    MPI_Reduce(&totwithin[0], &totgobal[0], k, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    double times = double( clock () - begin_time ) /  CLOCKS_PER_SEC;//no read & write
    double timecost;
    MPI_Reduce(&times, &timecost, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    //master write
    if(world_rank==0){
        string ss = "./" + writefile;
        sstmwt = ss+to_string(k)+".txt";
        ofstream myfile;
        myfile.open(sstmwt);
        for(auto i: labels){
            myfile << i << ",";
        }
            
        myfile<<"\n";
        myfile.close();
        
        sstmwt = ss+"within"+to_string(k)+".txt";
        myfile.open(sstmwt);

        for(int i=0; i<totgobal.size();i++){
            myfile << totgobal[i] << ",";
        }

        myfile<<"\n";
        myfile.close();
    }
    if(world_rank==0){
        cout << "total time cost: " << timecost << endl;
    }
    /*close the parallel section*/
    MPI_Finalize();
    return 0;
}

