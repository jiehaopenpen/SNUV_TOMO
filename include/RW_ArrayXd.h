#ifndef CT_2D_RW_ARRAYXD_H
#define CT_2D_RW_ARRAYXD_H

#include <string>
#include <cstdio>
#include <iostream>
using namespace std;

#include "../Eigen/Core"
using namespace Eigen;

bool saveArrayXd(ArrayXd &array, string fn){

    FILE *f = fopen(fn.c_str(), "wb");
    if (!f){
        return false;
    }

    int size = array.size();
    //cout<<size<<endl;
    fwrite(&size, sizeof(int32_t), 1, f);

    for (int k = 0; k < size; ++k){
        double val = array(k);
        fwrite(&val, sizeof(double), 1, f);
    }

    fclose(f);
    return true;
}

bool loadArrayXd(ArrayXd &array, string fn){

    FILE *f = fopen(fn.c_str(), "rb");
    if (!f){
        return false;
    }

    int size;
    fread(&size, sizeof(int32_t), 1, f);
    array.setZero(size);

    for (int k = 0; k < size; ++k){
        double val;
        fread(&val, sizeof(double), 1, f);
        array(k) = val;
    }

    fclose(f);
    return true;
}

bool loadArrayXi(ArrayXi &array, string fn){

    FILE *f = fopen(fn.c_str(), "rb");
    if (!f){
        return false;
    }

    int size;
    fread(&size, sizeof(int32_t), 1, f);
    array.setZero(size);

    for (int k = 0; k < size; ++k){
        int32_t val;
        fread(&val, sizeof(int32_t), 1, f);
        array(k) = val;
    }

    fclose(f);
    return true;
}
#endif //CT_2D_RW_ARRAYXD_H
