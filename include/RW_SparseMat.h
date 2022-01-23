#ifndef CT_2D_RW_SPARSEMAT_H
#define CT_2D_RW_SPARSEMAT_H

#include <string>
#include <cstdio>
#include <iostream>
#include <vector>
using namespace std;

#include "../Eigen/SparseCore"
using namespace Eigen;

bool save_triplets_bin(SparseMatrix<double, RowMajor> &sparseMat, string fn) {

    FILE *f = fopen(fn.c_str(), "wb");
    if (!f){
        return false;
    }

    int xyn[3] = {sparseMat.rows(), sparseMat.cols(), sparseMat.nonZeros()};
    fwrite(xyn, sizeof(int32_t), 3, f);

    for (int k = 0; k < sparseMat.outerSize(); ++k){
        for (SparseMatrix<double, RowMajor>::InnerIterator it(sparseMat,k); it; ++it){
            int32_t row_col[2] = {it.row(), it.col()};
            fwrite(row_col, sizeof(int32_t), 2, f);
            double val = it.value();
            fwrite(&val, sizeof(double), 1, f);
        }
    }

    fclose(f);
    return true;
}

bool load_triplets_bin(SparseMatrix<double, RowMajor> &sparseMat, string fn) {

    FILE *f = fopen(fn.c_str(), "rb");
    if (!f){
        return false;
    }

    int xyn[3];
    fread(xyn, sizeof(int32_t), 3, f);
    cout<<xyn[0]<<"  "<<xyn[1]<<"   "<<xyn[2]<<endl;
    sparseMat.resize(xyn[0], xyn[1]);
    vector<Triplet<double>> trips(xyn[2]);

    for (int k = 0; k < trips.size(); ++k){
        int32_t row_col[2];
        fread(row_col, sizeof(int32_t), 2, f);
        double val;
        fread(&val, sizeof(double), 1, f);

        trips[k] = Triplet<double>(row_col[0], row_col[1], val);
    }

    fclose(f);
    sparseMat.setFromTriplets(trips.begin(), trips.end());
    return true;
}


bool load_triplets_bin_c(SparseMatrix<double, ColMajor> &sparseMat, string fn) {

    FILE *f = fopen(fn.c_str(), "rb");
    if (!f){
        return false;
    }

    int xyn[3];
    fread(xyn, sizeof(int32_t), 3, f);
    cout<<xyn[0]<<"  "<<xyn[1]<<"   "<<xyn[2]<<endl;
    sparseMat.resize(xyn[0], xyn[1]);
    vector<Triplet<double>> trips(xyn[2]);

    for (int k = 0; k < trips.size(); ++k){
        int32_t row_col[2];
        fread(row_col, sizeof(int32_t), 2, f);
        double val;
        fread(&val, sizeof(double), 1, f);

        trips[k] = Triplet<double>(row_col[0], row_col[1], val);
    }

    fclose(f);
    sparseMat.setFromTriplets(trips.begin(), trips.end());
    return true;
}
#endif //CT_2D_RW_SPARSEMAT_H
