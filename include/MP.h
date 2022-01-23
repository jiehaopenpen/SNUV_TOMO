#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>


ArrayXd IRCD(SparseMatrix<double, ColMajor> *A, ArrayXd *proj, double SigmaE, double SigmaZ, double beta, int Strategy);
void UpdateX_l(ArrayXd& X, ArrayXXd* Sigma, SparseMatrix<double, ColMajor>* A, ArrayXd* proj, ArrayXd& Z, double SigmaE, double SigmaZ, int index);
void UpdateSigma(ArrayXd* X, ArrayXXd& Sigma, double SigmaE, int index);
void UpdateSigma_lp(ArrayXd* X, ArrayXXd& Sigma, double SigmaE, double beta, int index);
cv::Mat EigenToCV(ArrayXd *img);

ArrayXd IRCD(SparseMatrix<double, ColMajor> *A, ArrayXd* proj, double SigmaE, double SigmaZ, double beta, int Strategy){
    int L = (*A).cols();  // Number of pixels
    int N = (*A).rows();  // Number of measurements
    int B = 2*L - sett::Lx - sett::Ly; // Number of nearest neighbor pairs
    

    ArrayXd img;
    loadArrayXd(img, "../objects/"+ sett::objSelection);
    srand((unsigned int) time(0));
    ArrayXd X = ArrayXd::Constant(L,0.2);
    //X = img + 0.01*ArrayXd::Random(L);
    cout<<"Initialization RMSE: "<<(X-img).matrix().norm()/sqrt(X.size())<<endl;

    ArrayXd *Xp = &X;
    ArrayXd X_old = ArrayXd::Zero(L);
    ArrayXd Z = ArrayXd::Zero(N);

    // Initialize the auxilliary variable Z.
    for(int i = 0;i<L;i++){
        for (SparseMatrix<double,ColMajor>::InnerIterator it(*A,i); it; ++it){
            Z(it.index()) += X(i) * it.value();
        }
    }

    ArrayXXd Sigma = ArrayXXd::Constant(L,2,pow(10,-3)); // A Lx2 array to store unknow variance, there are reduant terms but enhance the readability 
    ArrayXXd* Sigmap = &Sigma;
    int l = 0;
      
    double sharpness = 0, sparse = 0;   // Sharpness, sparsity measure
    int iter = 0;
    for(iter = 0; iter < sett::max_num; iter++ ) {
        // Old vector memorization
        X_old = X;

        // Random update sequence
        srand ( unsigned ( std::time(0) ) );
        vector<int> myvector;
        for (int i=1; i<L; ++i) myvector.push_back(i);
        std::random_shuffle ( myvector.begin(), myvector.end() );

        switch (Strategy) {
            case 1: {                                               // No prior with batch update
                for (l = 0; l < L; l++) {
                    if(sett::Ran_order) UpdateX_l(X, Sigmap, A, proj, Z, SigmaE, SigmaZ, myvector[l]);
                    else UpdateX_l(X, Sigmap, A, proj, Z, SigmaE, SigmaZ, l);
                }
                for (l = 0; l < L; l++) {                           
                    UpdateSigma(Xp, Sigma, SigmaE, l);
                }
                break;
            }


            case 2: {                                               // Smoothed-Lp Norm using IRCD
                for (l = 0; l < L; l++) {
                    UpdateX_l(X, Sigmap, A, proj, Z, SigmaE, SigmaZ, l);
                }
                for (l = 0; l < L; l++) {                           
                    UpdateSigma_lp(Xp, Sigma, SigmaE, beta, l);
                }
                break;
            }

            default: {
                cout<<"No this kind of strategy!!!"<<endl;
                break;
            }
        }
        
        //cout<<Sigma.matrix().lpNorm<0>()<<endl;
        //cout << (Sigma != 0).count()<<endl;
        sharpness = Sigma.sum()/(Sigma != 0).count();
        sparse = (Sigma != 0).count()/(double)B;
        
        double diff = (X_old-X).matrix().norm()/sqrt(X.size());
        double rmse = (X-img).matrix().norm()/sqrt(X.size());

        //Doing some intermediate display and (if intermediate) save intermediate results 
        if( (iter+1) % sett::display == 0 ){
            cout<<"Iteration No. "<< iter+1 <<" Diff: "<< diff <<", RMSE: "<< rmse <<", Sharp: "<< sharpness <<", Sparse:"<< sparse <<endl;
            if (sett::intermediate) {
                cv::imwrite("../reconstructions/intermediate_results/SNUV_" + sett::objSelection + "_nProj=" + to_string(sett::nProj) + "_strategy=" + to_string(Strategy) + "_iter=" + to_string(iter+1) + ".png", EigenToCV(Xp));
            }
            
        }
        // If the difference between iterations is really small, break!
        if(diff<1*pow(10,-10)) {
            cout<<"Iteration No. "<<iter+1<<" Diff: "<<diff<<", RMSE: "<<rmse<<", Sharp: "<<sharpness<<", Sparse:"<<sparse<<endl;
            break;    
        }
    }
    return X;
}

void UpdateX_l(ArrayXd& X, ArrayXXd* Sigma, SparseMatrix<double, ColMajor> *A, ArrayXd* proj,ArrayXd& Z, double SigmaE, double SigmaZ, int index){
    // Coordinate descent of each element of the image 

	int r = index%sett::Lx;
	int c = index / sett::Lx;
    // Find closed-form solution of the least square problems (one possible implementation)
	double xi_fw = 0, xi_bw = 0, w_fw = 0, w_bw = 0;
	if (r != 0) {
		xi_fw += X(index - 1) / (SigmaE + (*Sigma)(index - 1, 1));
		w_fw += 1 / (SigmaE + (*Sigma)(index - 1, 1));
	}
	if (c != 0) {
		xi_fw += X(index - sett::Lx) / (SigmaE + (*Sigma)(index - sett::Lx, 0));
		w_fw += 1 / (SigmaE + (*Sigma)(index - sett::Lx, 0));
	}
	if (r != sett::Lx - 1) {
		xi_fw += X(index + 1) / (SigmaE + (*Sigma)(index, 1));
		w_fw += 1 / (SigmaE + (*Sigma)(index, 1));
	}
	if (c != sett::Ly - 1) {
		xi_fw += X(index + sett::Lx) / (SigmaE + (*Sigma)(index, 0));
		w_fw += 1 / (SigmaE + (*Sigma)(index, 0));
	}

	for (SparseMatrix<double, ColMajor>::InnerIterator it(*A, index); it; ++it) {
		w_bw += pow(it.value(), 2) / SigmaZ;
		xi_bw += it.value() * ((*proj)(it.index()) - Z(it.index()) + it.value() * X(index)) / SigmaZ;
	}
	double temp = (xi_fw + xi_bw) / (w_fw + w_bw);
	if (temp>1) temp = 1;
	else if (temp<0) temp = 0;
	double delta = X(index) - temp;
	X(index) = temp;

    // Update the auxilliary variable  Z
	for (SparseMatrix<double, ColMajor>::InnerIterator it(*A, index); it; ++it) {
		Z(it.index()) -= delta * it.value();
	}

}



void UpdateSigma(ArrayXd* X, ArrayXXd& Sigma, double SigmaE, int index){
    // Unknown variance update for Plain-SNUV
	int r = index%sett::Lx;
	int c = index / sett::Lx;

	if ((r != sett::Lx - 1) && (c != sett::Ly - 1)) {
		Sigma(index, 0) = max(0., pow((*X)(index) - (*X)(index + sett::Lx), 2) - SigmaE);
		Sigma(index, 1) = max(0., pow((*X)(index) - (*X)(index + 1), 2) - SigmaE);
	}
	else if ((r == sett::Lx - 1) && (c != sett::Ly - 1)) {
		Sigma(index, 0) = max(0., (pow((*X)(index) - (*X)(index + sett::Lx), 2) - SigmaE));
	}
	else if ((r != sett::Lx - 1) && (c == sett::Ly - 1)) {
		Sigma(index, 1) = max(0., pow((*X)(index) - (*X)(index + 1), 2) - SigmaE);
	}
	else {
		
	}


}


void UpdateSigma_lp(ArrayXd* X, ArrayXXd& Sigma, double SigmaE, double beta, int index){
    // Unknown variance update for Smoothed Lp Norms
    int r = index%sett::Lx;
    int c = index/sett::Lx;
    double mi = 2 - sett::p;

    if((r!=sett::Lx-1) && (c!=sett::Ly-1)){
        Sigma(index,0) = max(0., pow(abs((*X)(index)-(*X)(index+sett::Lx))/beta,mi)-SigmaE);
        Sigma(index,1) = max(0., pow(abs((*X)(index)-(*X)(index+1))/beta,mi)-SigmaE);
    }
    else if((r==sett::Lx-1) && (c!=sett::Ly-1)){
        Sigma(index,0) = max(0., pow(abs((*X)(index)-(*X)(index+sett::Lx))/beta,mi)-SigmaE);
    }
    else if((r!=sett::Lx-1) && (c==sett::Ly-1)){
        Sigma(index,1) = max(0., pow(abs((*X)(index)-(*X)(index+1))/beta,mi)-SigmaE);
    }
    else{
        
    }

}

cv::Mat EigenToCV(ArrayXd* img){                    // Transform from eigen array to opencv matrix for the same image
    cv::Mat M(sett::Lx,sett::Ly,CV_32FC1);
    int j = 0;
    for(int r = 0; r < sett::Lx; r++){
        for(int c = 0; c < sett::Ly; c++){
            j = r + sett::Lx * c;
            M.at<float>(r,c) = (*img)(j)*255;
        }
    }
    return M;
}