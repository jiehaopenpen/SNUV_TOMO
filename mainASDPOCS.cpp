#include <iostream>
#include <cmath>
#include <chrono>
using namespace std;

#include "Eigen/Core"
#include "Eigen/SparseCore"
using namespace Eigen;

#include "settings.h"
#include "include/RW_SparseMat.h"
#include "include/RW_ArrayXd.h"

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "include/MP.h"

ArrayXd ASDPOCS_algo_TV(ArrayXd proj, SparseMatrix<double, RowMajor> A, int iterNumber, int N_grad, double Epsilon,
                     double a, double a_red, double r_max, double Beta, double Beta_red, double eta = pow(10,-8));

ArrayXd ASDPOCS_algo_L1gradient(ArrayXd proj, SparseMatrix<double, RowMajor> A, int iterNumber, int N_grad, double Epsilon,
                     double a, double a_red, double r_max, double Beta, double Beta_red, double eta = pow(10,-8));

int main() {

    // Load the projection vector
    ArrayXd proj;
    loadArrayXd(proj, "../projections/proj_" + sett::objSelection + "_nProj=" + to_string(sett::nProj));

    // Load the projection matrix A
    SparseMatrix<double, RowMajor> A;
    load_triplets_bin(A, "../projections/matrixProj_" + sett::objSelection + "_nProj=" + to_string(sett::nProj));
    A.makeCompressed();

    cout << "Projection vector y and projection matrix A loaded.\n";

    // Parameters for the reconstruction
    int iterNumber = 8000;
    int N_grad = 20;
    double Epsilon = pow(10, -1);
    double a = 0.2;
    double a_red = 0.999;
    double Beta = 2;
    double Beta_red = 0.9999;
    double r_max = 0.4;

    // Reconstruction with ASD-POCS (constrained TV minimization)
    auto start = chrono::high_resolution_clock::now();
    //ArrayXd reconImg = ASDPOCS_algo_TV(proj, A, iterNumber, N_grad, Epsilon, a, a_red, r_max, Beta, Beta_red);
    ArrayXd reconImg = ASDPOCS_algo_L1gradient(proj, A, iterNumber, N_grad, Epsilon, a, a_red, r_max, Beta, Beta_red);
    auto end = chrono::high_resolution_clock::now();
    cout << "Timer: "<< chrono::duration_cast<chrono::seconds>(end - start).count() <<'\n';

    ArrayXd img;
    loadArrayXd(img, "../objects/"+ sett::objSelection);
    cout << "Reconstruction RMSE = " << (reconImg - img).matrix().norm()/sqrt(A.innerSize()) << endl;
    // Save the reconstruction
    saveArrayXd(reconImg, "../reconstructions/ASDPOCS_" + sett::objSelection + "_nProj=" + to_string(sett::nProj));
    cv::imwrite("../reconstructions/ASDPOCS_" + sett::objSelection + "_nProj=" + to_string(sett::nProj) + ".png", EigenToCV(&reconImg));

    return 0;
}



/* ASD-POCS algorithm to minimize the image TV under positivity constraint and data-fidelity constraint (norm(y - Ax) < Epsilon)
 *
 * Input:   proj        <-  projection vector (observations)
 *          A           <-  projection matrix (obtained by distance-driven projection)
 *          iterNumber  <-  number of iterations
 *          N_grad      <-  number of iterations in the ASD step
 *          Epsilon     <-  bound in the constraint norm(y - Ax) < Epsilon
 *          a           <-  quantity used to initialize the gradient descent step size (dtvg)
 *          a_red, r_max    <-  parameters controlling the evolution of dtvg
 *          Beta, Beta_red  <-  Beta is the ART-relaxation parameter, Beta_red is the factor used to reduce Beta
 *                              after each ASDPOCS iteration
 *          eta         <-  small quantity added to the gradient of the image TV to avoid zero denominators
 *
 * Output:  f_res   <-  reconstructed image obtained after the POCS step of the last iteration
 *
 */

ArrayXd ASDPOCS_algo_TV(ArrayXd proj, SparseMatrix<double, RowMajor> A, int iterNumber, int N_grad, double Epsilon,
                        double a, double a_red, double r_max, double Beta, double Beta_red, double eta){

    // Some useful constants
    int L = A.innerSize();        // Number of pixels
    int N = A.outerSize();  // Number of observations

    // Compute the diagonal coefficients of the matrix A * A^T
    ArrayXd AAt;
    AAt.setZero(N);
    for(int n = 0; n < N; ++n){
        AAt.coeffRef(n) = A.row(n).squaredNorm();
    }

    ArrayXd f;
    f.setZero(L);
    ArrayXd f_res;
    double dtvg;
    double dg;
    int idx;

    for(int i = 0; i < iterNumber; ++i){
        ArrayXd f0 = f;

        // >> POCS step
        for(int n = 0; n < N; ++n){     // Consistency to the data (ART with relaxation parameter Beta)
            double error = proj(n);
            for (SparseMatrix<double, RowMajor>::InnerIterator it(A,n); it; ++it) {
                error -= it.value() * f(it.index());
            }
            for (SparseMatrix<double, RowMajor>::InnerIterator it(A,n); it; ++it) {
                f(it.index()) += Beta * it.value() * error/AAt(n);
            }
        }

        f = (f <= 0).select(0, f);      // Positivity constraint
        f = (f >= 1).select(1, f);
        // << POCS step

        f_res = f;

        ArrayXd g = ArrayXd::Zero(N);
        for(int n = 0; n < N; ++n){
            for (SparseMatrix<double, RowMajor>::InnerIterator it(A,n); it; ++it) {
                g(n) += it.value() * f(it.index());
            }
        }

        double dd = (g - proj).matrix().norm();
        double dp = (f - f0).matrix().norm();

        if(i == 0){
            dtvg = a * dp;
        }

        // >> ASD step
        ArrayXd df;
        for(int n_grad = 0; n_grad < N_grad; ++n_grad){
            // >> Compute the approximated gradient of the TV norm and store it in df
            df.setZero(L);
            idx = 0;

            // First column
            df.coeffRef(idx) = - (f(idx+1)-f(idx))/sqrt(pow(f(idx+1)-f(idx),2) + eta)
                               - (f(idx+sett::Lx)-f(idx))/sqrt(pow(f(idx+sett::Lx)-f(idx),2) + eta);
            ++idx;

            for(int lx = 1; lx < sett::Lx-1; ++lx){
                df.coeffRef(idx) = (f(idx) - f(idx-1))/sqrt(pow(f(idx) - f(idx-1),2) + eta )
                                   - (f(idx+1)-f(idx))/sqrt(pow(f(idx+1)-f(idx),2) + eta)
                                   - (f(idx+sett::Lx)-f(idx))/sqrt(pow(f(idx+sett::Lx)-f(idx),2) + pow(f(idx+sett::Lx)-f(idx-1+sett::Lx),2) + eta);
                ++idx;
            }

            df.coeffRef(idx) = (f(idx) - f(idx-1))/sqrt(pow(f(idx) - f(idx-1),2) + eta )
                               - (f(idx+sett::Lx)-f(idx))/sqrt(pow(f(idx+sett::Lx)-f(idx),2) + pow(f(idx+sett::Lx)-f(idx-1+sett::Lx),2) + eta);
            ++idx;

            // Second to second-to-last columns
            for(int ly = 1; ly < sett::Ly-1; ++ly){
                df.coeffRef(idx) = (f(idx) - f(idx-sett::Lx))/sqrt(pow(f(idx) - f(idx-sett::Lx),2) + eta)
                                   - (f(idx+1)-f(idx))/sqrt(pow(f(idx+1)-f(idx),2) + pow(f(idx+1)-f(idx+1-sett::Lx),2) + eta)
                                   - (f(idx+sett::Lx)-f(idx))/sqrt(pow(f(idx+sett::Lx)-f(idx),2) + eta);
                ++idx;

                for(int lx = 1; lx < sett::Lx-1; ++lx){
                    df.coeffRef(idx) = (2*f(idx) - f(idx-1) - f(idx-sett::Lx))/sqrt(pow(f(idx) - f(idx-1),2) + pow(f(idx) - f(idx-sett::Lx),2) + eta)
                                       - (f(idx+1)-f(idx))/sqrt(pow(f(idx+1)-f(idx),2) + pow(f(idx+1)-f(idx+1-sett::Lx),2) + eta)
                                       - (f(idx+sett::Lx)-f(idx))/sqrt(pow(f(idx+sett::Lx)-f(idx),2) + pow(f(idx+sett::Lx)-f(idx-1+sett::Lx),2) + eta);
                    ++idx;
                }

                df.coeffRef(idx) = (2*f(idx) - f(idx-1) - f(idx-sett::Lx))/sqrt(pow(f(idx) - f(idx-1),2) + pow(f(idx) - f(idx-sett::Lx),2) + eta)
                                   - (f(idx+sett::Lx)-f(idx))/sqrt(pow(f(idx+sett::Lx)-f(idx),2) + pow(f(idx+sett::Lx)-f(idx-1+sett::Lx),2) + eta);
                ++idx;
            }

            // Last column
            df.coeffRef(idx) = (f(idx) - f(idx-sett::Lx))/sqrt(pow(f(idx) - f(idx-sett::Lx),2) + eta)
                               - (f(idx+1)-f(idx))/sqrt(pow(f(idx+1)-f(idx),2) + pow(f(idx+1)-f(idx+1-sett::Lx),2) + eta);
            ++idx;

            for(int lx = 1; lx < sett::Lx-1; ++lx){
                df.coeffRef(idx) = (2*f(idx) - f(idx-1) - f(idx-sett::Lx))/sqrt(pow(f(idx) - f(idx-1),2) + pow(f(idx) - f(idx-sett::Lx),2) + eta )
                                   - (f(idx+1)-f(idx))/sqrt(pow(f(idx+1)-f(idx),2) + pow(f(idx+1)-f(idx+1-sett::Lx),2) + eta);
                ++idx;
            }

            df.coeffRef(idx) = (2*f(idx) - f(idx-1) - f(idx-sett::Lx))/sqrt(pow(f(idx) - f(idx-1),2) + pow(f(idx) - f(idx-sett::Lx),2) + eta );
            // << Compute the approximated gradient of the TV norm and store it in df

            f -= dtvg * (df.matrix().normalized().array());
        }
        // << ASD step

        dg = (f - f_res).matrix().norm();
        if( (dg > r_max*dp) && (dd > Epsilon)){
            dtvg *= a_red;      // Reduce gradient descent step size
        }
        Beta *= Beta_red;       // Reduce ART relaxation parameter

        if(i%100 == 0 || i == iterNumber-1) cout << "Iteration " << (i + 1) << endl;

    }

    return f_res;

}


// Modified ASD-POCS algorithm (the cost minimized under constraint is the sum of the absolute jumps along the lines and the columns)
inline double h_derivative(double x, double eta){
    double res;
    if ( x < -eta){
        res = -1;
    } else if(eta < x){
        res = 1;
    } else{
        res = x/eta;
    }

    return res;
}

ArrayXd ASDPOCS_algo_L1gradient(ArrayXd proj, SparseMatrix<double, RowMajor> A, int iterNumber, int N_grad, double Epsilon,
                                double a, double a_red, double r_max, double Beta, double Beta_red, double eta){

    int L = A.innerSize();      // Number of pixels
    int N = A.outerSize();      // Number of observations

    // Compute the diagonal coefficients of the matrix A * A^T
    ArrayXd AAt;
    AAt.setZero(N);
    for(int n = 0; n < N; ++n){
        AAt.coeffRef(n) = A.row(n).squaredNorm();
    }

    ArrayXd f;
    f.setZero(L);
    ArrayXd f_res;
    double dtvg;
    double dg;
    int idx;

    for(int i = 0; i < iterNumber; ++i){
        ArrayXd f0 = f;

        // >> POCS step
        for(int n = 0; n < N; ++n){     // Consistency to the data (ART with relaxation parameter Beta)
            double error = proj(n);
            for (SparseMatrix<double, RowMajor>::InnerIterator it(A,n); it; ++it) {
                error -= it.value() * f(it.index());
            }
            for (SparseMatrix<double, RowMajor>::InnerIterator it(A,n); it; ++it) {
                f(it.index()) += Beta * it.value() * error/AAt(n);
            }
        }

        f = (f <= 0).select(0, f);      // Positivity constraint
        // << POCS step

        f_res = f;

        ArrayXd g = ArrayXd::Zero(N);
        for(int n = 0; n < N; ++n){
            for (SparseMatrix<double, RowMajor>::InnerIterator it(A,n); it; ++it) {
                g(n) += it.value() * f(it.index());
            }
        }

        double dd = (g - proj).matrix().norm();
        double dp = (f - f0).matrix().norm();

        if(i == 0){
            dtvg = a * dp;
        }

        // >> ASD step
        ArrayXd df;
        for(int n_grad = 0; n_grad < N_grad; ++n_grad){
            // >> Compute the approximated gradient of the TV norm and store it in df
            df.setZero(L);
            idx = 0;

            // First column
            df.coeffRef(idx) = - h_derivative(f(idx+1) - f(idx), eta) - h_derivative(f(idx+sett::Lx) - f(idx), eta);
            ++idx;

            for(int lx = 1; lx < sett::Lx-1; ++lx){
                df.coeffRef(idx) = h_derivative(f(idx) - f(idx-1),eta) - h_derivative(f(idx+1) - f(idx), eta)
                                   - h_derivative(f(idx+sett::Lx) - f(idx), eta);
                ++idx;
            }

            df.coeffRef(idx) = h_derivative(f(idx) - f(idx-1), eta) - h_derivative(f(idx+sett::Lx) - f(idx), eta);
            ++idx;

            // Second to second-to-last columns
            for(int ly = 1; ly < sett::Ly-1; ++ly){
                df.coeffRef(idx) = h_derivative(f(idx) - f(idx-sett::Lx), eta) - h_derivative(f(idx+1) - f(idx), eta)
                                   - h_derivative(f(idx+sett::Lx) - f(idx), eta);
                ++idx;

                for(int lx = 1; lx < sett::Lx-1; ++lx){
                    df.coeffRef(idx) = h_derivative(f(idx) - f(idx-1), eta) + h_derivative(f(idx) - f(idx-sett::Lx), eta)
                                       - h_derivative(f(idx+1) - f(idx), eta) - h_derivative(f(idx+sett::Lx) - f(idx), eta);
                    ++idx;
                }

                df.coeffRef(idx) = h_derivative(f(idx) - f(idx-1), eta) + h_derivative(f(idx) - f(idx-sett::Lx), eta)
                                   - h_derivative(f(idx+sett::Lx) - f(idx), eta);
                ++idx;
            }

            // Last column
            df.coeffRef(idx) = h_derivative(f(idx) - f(idx-sett::Lx), eta) - h_derivative(f(idx+1) - f(idx), eta);
            ++idx;

            for(int lx = 1; lx < sett::Lx-1; ++lx){
                df.coeffRef(idx) = h_derivative(f(idx) - f(idx-1), eta) + h_derivative(f(idx) - f(idx-sett::Lx), eta)
                                   - h_derivative(f(idx+1) - f(idx), eta);
                ++idx;
            }

            df.coeffRef(idx) = h_derivative(f(idx) - f(idx-1), eta) + h_derivative(f(idx) - f(idx-sett::Lx), eta);
            // << Compute the approximated gradient of the TV norm and store it in df

            f -= dtvg * (df.matrix().normalized().array());
        }
        // << ASD step

        dg = (f-f_res).matrix().norm();
        if( (dg > r_max*dp) && (dd > Epsilon)){
            dtvg *= a_red;      // Reduce gradient descent step size
        }
        Beta *= Beta_red;       // Reduce ART relaxation parameter

        if(i%100 == 0 || i == iterNumber-1) cout << "Iteration " << (i + 1) << endl;

    }
    for(int got=0;got<L;got++){
            if (f_res(got)<0){
                f_res(got)  = 0;
                
            }
            if (f_res(got) > 1){
               f_res(got)= 1;
            }
        }

    return f_res;

}

