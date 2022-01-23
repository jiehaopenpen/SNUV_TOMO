//
// Created by boxiao on 22.05.18.
//

#ifndef TOMO_SLS_SPARSEMATRIX_H
#define TOMO_SLS_SPARSEMATRIX_H


#include <string>
#include <cstdio>
#include <iostream>
#include <vector>
using namespace std;

#include "../Eigen/SparseCore"
using namespace Eigen;



typedef struct SM{
    double* vals;
    int* inner;
    int* outer;
    int NNZ;
    int Nrows;
    int Ncols;

}SM;

typedef struct SM_c{
    double* vals;
    int* inner;
    int* outer;
    int NNZ;
    int Nrows;
    int Ncols;

}SM_c;

bool load_csr(SM &A,string fn) {

    FILE *f = fopen(fn.c_str(), "rb");
    if (!f){
        return false;
    }

    int xyn[3];
    fread(xyn, sizeof(int32_t), 3, f);

    A.Nrows = xyn[0];
    A.Ncols = xyn[1];
    A.NNZ = xyn[2];
    A.vals=(double*)malloc(xyn[2]*sizeof(double));
    A.inner = (int *)malloc(xyn[2]*sizeof(int32_t));
    A.outer = (int *)malloc(xyn[0]*sizeof(int32_t));

    int j = 0;
    //cout<<outer[j]<<endl;
    //cout<<vals<<endl;
    for (int k = 0; k < xyn[2]; ++k){
        int32_t row_col[2];
        fread(row_col, sizeof(int32_t), 2, f);
        double val;
        fread(&val, sizeof(double), 1, f);

        if(row_col[0]!=j){

            //if(row_col[0]==5) cout<<j<<" adsfadsf "<<k<<endl;
            while(j< row_col[0]-1){

                //cout<<"????"<<endl;
                A.outer[j] = k;
                j++;
            }

            A.outer[j] = k;
            j++;
        }

        A.inner[k] = row_col[1];
        A.vals[k] = val;
        //row_old = row_col[0];
        //cout<<outer[j]<<endl;

    }
    //outer[j] = xyn[2];
    while (j<=xyn[0]-1){
        A.outer[j] = xyn[2];
        j++;
    }
    return true;
}


bool load_csc(SM_c &A, string fn) {

    FILE *f = fopen(fn.c_str(), "rb");
    if (!f){
        return false;
    }

    int xyn[3];
    fread(xyn, sizeof(int32_t), 3, f);

    A.Nrows = xyn[0];
    A.Ncols = xyn[1];
    A.NNZ = xyn[2];
    A.vals=(double*)malloc(xyn[2]*sizeof(double));
    A.inner = (int *)malloc(xyn[2]*sizeof(int32_t));
    A.outer = (int *)malloc(xyn[1]*sizeof(int32_t));

    int j = 0;
    //cout<<outer[j]<<endl;
    //cout<<vals<<endl;
    for (int k = 0; k < xyn[2]; ++k){
        int32_t row_col[2];
        fread(row_col, sizeof(int32_t), 2, f);
        double val;
        fread(&val, sizeof(double), 1, f);

        if(row_col[1]!=j){

            //if(row_col[0]==5) cout<<j<<" adsfadsf "<<k<<endl;
            while(j< row_col[1]-1){

                //cout<<"????"<<endl;
                A.outer[j] = k;
                j++;
            }

            A.outer[j] = k;
            j++;
        }

        A.inner[k] = row_col[0];
        A.vals[k] = val;
        //row_old = row_col[0];
        //cout<<outer[j]<<endl;

    }
    //outer[j] = xyn[2];
    while (j<=xyn[1]-1){
        A.outer[j] = xyn[2];
        j++;
    }
    return true;
}


#endif //TOMO_SLS_SPARSEMATRIX_H
