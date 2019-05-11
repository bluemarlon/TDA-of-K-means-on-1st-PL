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


//read one line of the file
vector<string> Readline(string line){
    stringstream   lineStream(line);
    vector<string> row = {};
    string         cell;
    while(getline(lineStream, cell, ' ')){
        row.push_back(cell);
    }
    return row;
}

//read the whole txt file
void Readtxt(
    string sstm,
    vector<double>& X,
    int& dim, //the dimension of the observations
    int n,
    //int valsum //the cumsum of previous counts = the starting point
    vector<int>& vet //the index of obs; the index should be in order 
) {
    ifstream mfile (sstm);
    string line;
    if(mfile.is_open()){
        //cout << "file is open" << endl;
        //the .txt file has colnames
        getline(mfile, line);       
        vector<string> row = Readline(line);
        dim  = row.size();//txt files
        int count = 0;
        for(auto ncount: vet){
            while( count!=ncount ){ count++;getline(mfile, line);}
            row = Readline(line);
            for(int i=0; i < dim; i++) {
                //row[i+1] as .txt file also has rownames;
                X.push_back(stod(row[i+1].c_str()));//stod to float
            }
        }
        mfile.close();
    }
}
void Readindex(
    string sstm,
    vector<int>& index,
    int n,
    int valsum //the cumsum of previous counts = the starting point
){
    ifstream mfile (sstm);
    string line;
    if(mfile.is_open()){
        for(int count = 0; count < valsum; count ++ ){
            getline(mfile, line);
        }
        int ncount = 1;
        while(getline(mfile, line) && ncount <= n)
        {   
            index.push_back(stoi(line));
            ncount++;
        }
        mfile.close();
    }
}

vector<double> stoffer(
    vector<double> x, //the categorical t.s. with 0 padding
    int num //the length after 0 padding
    ){
        // it outputs the WFT
    vector<short> Ipower(15, 0); //the size of ts no larger than 2^15
    vector<double> y(num, 0);
    for(int i=0; i<num; i++){
        int IB = i, IL = 0;
        while(true){
            int IBD = IB/2;
            if(IB == IBD * 2){
                Ipower[IL] = 0;
            }else{
                Ipower[IL] = 1;
            }
            if(IB==0 || IBD == 0){
                break;
            }
            IB = IBD;
            IL++;
        }
        int IP = 0;
        int IFAC = num;
        for(int t1=0; t1 < IL + 1; t1++){
            IFAC = IFAC/2;
            IP += IFAC * Ipower[t1];
        }
        y[IP] = x[i];
    }
    x.resize(0);
    x.insert(x.end(), y.begin(), y.end());
    int Iter((int)(log(num)/log(2)));
    for(int M=0; M < Iter; M++){
        int nump;
        if(M==0){
            nump = 1;
        }else{
            nump *= 2;
        }
        int Mnum = num/nump;
        int Mnum2 = Mnum/2;
        int alph = 1;
        for(int MP = 0; MP < nump; MP++){
            int IB = MP * Mnum;
            for(int MP2 = 0; MP2 < Mnum2; MP2++){
                int mnum21 = Mnum2 + MP2 + IB;
                int IBA = IB + MP2;
                y[IBA] = x[IBA] + alph * x[mnum21];
                y[mnum21] = x[IBA] - alph * x[mnum21];
            }
            alph = -alph;
        }
        double r = 1/pow(num, 0.5);
        for(int i=0; i<num; i++){
            x[i] = y[i] * r;
        }
    }
    return y;
}


int main(
	int argc, char * argv[]
    )
{  
    /*takes PLs as inputs into the K-means algorithm and outputs the cluster labels and total within cluster sum of squares. */
    // reads categories.txt of the raw data. categories.txt contains both rownames and columnnames of 250883 * 1441 matrix
    // reads sample.txt containing random indexes when assigning time series
    // write out: WFT#.txt and land#.txt for each processor and processor 0 output (cout) the computing time for WFT and 1st landscapes
    int world_rank, world_size;
    MPI_Init( &argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    string sstmrd; /*data of reading the file*/
    string sstmwt; /*data of writing the file*/
    vector<double> X; /*saving the raw time series based on the random order*/
    sstmrd += "./categories.txt";
    int dim;
    int n = 2508;
    if(world_rank==99){
        n = 2590;
    }
    // random indexes when assigning time series to processors
    string sstindex = "./sample.txt"; 
    vector<int> index; /*it saves the random indexes*/
    Readindex(sstindex, index, n, world_rank*2508);
    Readtxt(sstmrd, X, dim, n, index);
    clock_t begin_time = clock();
    //start to make features;
    //1. WFT
    vector<double> WFTs;
    int next2 = (int) pow(2,(int)(log(dim)/log(2))+1);
    for(int i=0;i< n; i++){
        vector<double> x(X.begin()+i*dim, X.begin()+(i+1)*dim);
        vector<double> x1(next2-dim, 0);
        //0 padding;
        x.insert(x.end(),x1.begin(), x1.end());
        vector<double> WFT = stoffer(x, next2);
        WFTs.insert(WFTs.end(), WFT.begin(), WFT.end());
    }
    // computing time without read & write
    double time1 = double( clock () - begin_time ) /  CLOCKS_PER_SEC;
    double WFTcost;
    MPI_Reduce(&time1, &WFTcost, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if(world_rank==0){
        cout << "WFTcost: " << WFTcost << endl;
    }
    sstmwt = "./WFT"+to_string(world_rank+1)+".txt";
    ofstream myfile;
    myfile.open(sstmwt);
    for(int j =0;j< dim-1;j++){
        myfile << "V" << (j+1) << ",";
    }
    myfile << "V" << (dim-1) <<"\n";
    for(int i=0;i<n;i++){
        for(int j =0;j< dim-1;j++){
            myfile << WFTs[i*dim+j] << ",";
        }
        myfile << WFTs[i*dim+dim-1]<<"\n";
    }
    myfile.close();
    MPI_Barrier(MPI_COMM_WORLD);
    //2. 1st landscapes;
    begin_time = clock();
    //give L=100 as the size of the 1st landscape;
    int L=100;
    //find the max and min of each WFT;
    vector<double> vmax, vmin;
    for(int i=0;i<n;i++){
        double mmax = -1000, mmin = 1000;
        for(int j =0;j< dim;j++){
            mmax = max(mmax, WFTs[i*dim+j]);
            mmin = min(mmin, WFTs[i*dim+j]);
        }
        vmax.push_back(mmax);
        vmin.push_back(mmin);
    }
    //presetting the range of the landscapes
    double rmax = 1e-13, rmin = -1e-13;
    //grid size of the lanscapes
    double step = (rmax-rmin)/(L-2+1);
    vector<double> lands; /*saves the 1st landscapes*/
    for(int i=0;i<n;i++){
        vector<double> land;
        for(double j=0;j<L;j++){
            land.push_back(max(min(rmin+step*(j) - vmin[i], vmax[i] - rmin - step*(j)), (double)0));
        }
        lands.insert(lands.end(), land.begin(), land.end());
    }
    // the computation time of the 1st landscapes without read & write
    time1 = double( clock () - begin_time ) /  CLOCKS_PER_SEC;
    double Landcost;
    MPI_Reduce(&time1, &Landcost, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if(world_rank==0){
        cout << "Landcost: " << Landcost << endl;
    }
    sstmwt = "./land"+to_string(world_rank+1)+".txt";
    myfile.open(sstmwt);
    for(int j =0;j< L-1;j++){
        myfile << "V" << (j+1) << ",";
    }
    myfile << "V"<< (L)<< "\n";
    for(int i=0;i<n;i++){
        for(int j =0;j< L-1;j++){
            myfile << lands[i*L+j] << ",";
        }
        myfile << lands[i*L + L-1] << "\n";
    }
    myfile.close();
    /*close the parallel section*/
    MPI_Finalize();
    return 0;
}
