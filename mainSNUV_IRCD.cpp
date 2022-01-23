#include <iostream>
#include <fstream>
#include <cmath>
#include "include/SparseMatrix.h"
#include "include/RW_ArrayXd.h"
#include "include/RW_SparseMat.h"
#include "settings.h"
using namespace std;
#include "include/MP.h"



int main() {

    // Load projection matrix A
    SparseMatrix<double, ColMajor> A;
    SparseMatrix<double, ColMajor> *Ap = &A;
    if(!load_triplets_bin_c(A,"../projections/matrixProj_" + sett::objSelection + "_nProj=" + to_string(sett::nProj)))
        cout<<"Matrix A Loading Failed!!!"<<endl;

    //Load projections y
    ArrayXd proj;
    ArrayXd* projp = &proj;
    if(!loadArrayXd(proj, "../projections/proj_" + sett::objSelection + "_nProj=" + to_string(sett::nProj)))
        cout<<"Projection y Loading Failed!!!"<<endl;
    
    if(sett::visulize_original){ //generate a png file to visualize the original raw image
        ArrayXd img;
        if(!loadArrayXd(img, "../objects/"+ sett::objSelection))
            cout<<"Image X Loading Failed!!!"<<endl;
        ifstream infile("../objects/"+ sett::objSelection + ".png"); //prevent from overwriting.
        if (!infile.good()) 
            cv::imwrite("../objects/"+ sett::objSelection + ".png", EigenToCV(&img));
    }
    //Parameters
    double Sigma0 = 1*pow(10,-4);     // Variance Sigma0^2
    double SigmaZ = 5*pow(10,-4);     // Variance SigmaZ^2
    double beta = 1;                  // For Lp Norm only
    int Strategy = 1;                 // Strategy, 1 for Plain-SNUV, 2 for Smoothed-Lp Norm

    ArrayXd ReconImg = IRCD(Ap,projp,Sigma0,SigmaZ,beta,Strategy);

    // Write the reconstruction in png format for visualization
    cv::imwrite("../reconstructions/SNUV_" + sett::objSelection + "_nProj=" + to_string(sett::nProj) + "_strategy=" + to_string(Strategy) + ".png", EigenToCV(&ReconImg));

    // Write the reconstruction in raw data for better precision
    saveArrayXd(ReconImg, "../reconstructions/SNUV_" + sett::objSelection + "_nProj=" + to_string(sett::nProj) + "_strategy=" + to_string(Strategy));
    
    return 0;
}



