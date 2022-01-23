#define _USE_MATH_DEFINES // includes mathematical constants defined in cmath


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
#include "include/convert_ArrayXd_CSV.h"
using namespace sett;


// 2D FAN-BEAM GEOMETRY
// ROTATING OBJECT, SOURCE AND LINE DETECTOR FIXED

/********************************************************************
 *                   DISTANCE-DRIVEN PROJECTION                     *
 ********************************************************************/

pair<ArrayXd, SparseMatrix<double, RowMajor>> distanceDrivenProjection(ArrayXd img, int iview);
pair<ArrayXd, SparseMatrix<double, RowMajor>> distanceDrivenProj_YSlices(ArrayXd img, double angle, int incIdxPix=1, int idxPixStart=0, int c1=0, int c2=1);
pair<ArrayXd, SparseMatrix<double, RowMajor>> distanceDrivenProj_XSlices(ArrayXd img, double angle, int incIdxPix=1, int idxPixStart=0, int c1=0, int c2=1);

pair<ArrayXd, SparseMatrix<double, RowMajor>> distanceDrivenProjection_h(ArrayXd img, int iview);
pair<ArrayXd, SparseMatrix<double, RowMajor>> distanceDrivenProj_YSlices_h(ArrayXd img, double angle, int incIdxPix=1, int idxPixStart=0, int c1=0, int c2=1);
pair<ArrayXd, SparseMatrix<double, RowMajor>> distanceDrivenProj_XSlices_h(ArrayXd img, double angle, int incIdxPix=1, int idxPixStart=0, int c1=0, int c2=1);

     double dxPix, dyPix;  
     ArrayXd xBoundaryPix, yBoundaryPix, xCenterPix, yCenterPix;

int main() {

    // Open the image (the loaded image has been flattened column-major)
    ArrayXd img(Lx*Ly*factor*factor);
    if(!loadArrayXd(img, "../objects/"+ objSelection)){
        cout<<"Raw Image Loading Failed!!! Please Run the transformer script first!!!"<<endl;
    };
    
    //cout<<img.sum()<<endl;
    // Initialize the projection vector and the sparse matrix A: proj = A*img in practice
    ArrayXd projc(Nu * nProj),projh(Nu * nProj);
    SparseMatrix<double, RowMajor> A(Nu*nProj, Lx*Ly);
    SparseMatrix<double, RowMajor> A1(Nu*nProj, Mxh*Myh);

    //cout<<"Projecting the coarse resolution grid"<<endl;
            /*ArrayXd img1;
        img1.setZero(img.size()/pow(factor,2));
        img1 +=0.5;*/
        ArrayXd img1 = ArrayXd::Zero(Lx*Ly);
    for(int c = 0; c <  Ly*factor; c++){
        for(int r = 0; r <  Lx*factor ; r++){
            int j = r +  Lx*factor * c;    
                img1((r/factor)+Lx*floor(c/factor))+=img(j);

        }
    }
    img1 /= (factor*factor);

    
    for(int i = 0; i < nProj; i++){    

        pair<ArrayXd, SparseMatrix<double, RowMajor>> res = distanceDrivenProjection(img1, i);
        //cout<<res.first<<endl;
        projc.segment<Nu>(i*Nu) = res.first;
        A.middleRows(i*Nu,Nu) = res.second;
        //cout << i << endl;
    }
    A.makeCompressed();
    //cout<<projc<<endl;
    save_triplets_bin(A, "../projections/matrixProj_" + objSelection + "_nProj=" + to_string(nProj));
    saveArrayXd(projc, "../projections/proj_" + objSelection + "_nProj=" + to_string(nProj));
    A.resize(0,0);

//     cout<<"Projecting the fine resolution grid"<<endl;
//     for(int i = 0; i < nProj; i++){    

//         //cout<<img1.size()<<endl;
//         //img.setZero(img.size());
//         //img +=0.5;
//         //cout<<img.size()<<endl;
//         pair<ArrayXd, SparseMatrix<double, RowMajor>> res = distanceDrivenProjection_h(img, i);
//         projh.segment<Nu>(i*Nu) = res.first;
//         A1.middleRows(i*Nu,Nu) = res.second;
//         //cout << i << endl;
//     }
//     A1.makeCompressed();

 
//     /*if(AddNoiseToProj == true){
    
//     unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//     default_random_engine generator (seed);
//     normal_distribution<double> distribution (0.0,SigmaZ);
//     for (int j=0; j<projh.size(); ++j)
//         projh(j)=projh(j)+distribution(generator);
//     }
//     int ini = 10000;
//     for(int j = ini;j < ini+100;j++){
//         cout<<projh(j)<<"    "<<projc(j)<<endl;
//     } */   

//     // Store the projection vector proj and the system matrix A into binary files
//     cout<<abs(projc.sum())/projc.size()<<endl;
//    // cout<<"Saving the results with error:"<<(projc-projh).matrix().abs()/sqrt(projc.size())<<endl;
//     //cout<<"Saving the results with error:"<<abs((projc-projh)).sum()/projc.size()<<endl;
//     //cout<<projc<<endl;
//     saveArrayXd(projc, "../projections/proj_" + objSelection + "_nProj=" + to_string(nProj));
//     save_triplets_bin(A1, "../projections/matrixProj_" + objSelection + "h_nProj=" + to_string(nProj));

    cout<<"Projection Done"<<endl;

    return 0;
}




/* Function distanceDrivenProjection(ArrayXd img, int iview):
 *
 *  Depending on the view angle it chooses the good way to perform the distance-driven projection
 *  i.e. it uses distanceDrivenProj_YSlices if the projection has to be done column-wise and
 *  distanceDrivenProj_XSlices if the projection has to be done row-wise.
 *
 * Inputs:  img     <-  image to project
 *          iview   <-  index (in the array deg) of the rotation angle of the image
 *
 * Outputs: projection  <-  distance-driven projection of the image for the considered view angle
 *          A           <-  corresponding projection matrix i.e. projection = A * img
 *
 */
pair<ArrayXd,SparseMatrix<double, RowMajor>> distanceDrivenProjection(ArrayXd img, int iview){
    double angle = deg(iview);
    cout<<iview<<endl;
    if((-45 <= angle) && (angle <=45)){

        return distanceDrivenProj_YSlices(img, angle);


    } else if((45 < angle) && (angle < 135)){

        return distanceDrivenProj_XSlices(img, angle, -1, Ly, 1, 0);

    } else if(((135 <= angle) && (angle<=180)) || ((-180 <= angle) && (angle <= -135))){

        return distanceDrivenProj_YSlices(img, angle, -1, Lx, 1, 0);

    } else{

        return distanceDrivenProj_XSlices(img, angle);

    }
}



/* Function distanceDrivenProj_YSlices(ArrayXd img, double angle, int incIdxPix, int idxPixStart, int c1, int c2):
 *
 *  Perform the distance-driven projection by processing the image column-wise.
 *
 *  Inputs: img     <-  image to project
 *          angle   <-  rotation angle of the image (in degrees)
 *          incIdxPix, idxPixStart, c1, c2
 *                  <-  variables to browse the boundary projections and manage the pixel index correctly
 *                      when doing the distance-driven projection
 *
 *  Outputs:    projection  <-  distance-driven projection of the image for the considered angle
 *              A           <-  corresponding projection matrix i.e. projection = A * img
 *
 */
pair<ArrayXd, SparseMatrix<double, RowMajor>> distanceDrivenProj_YSlices(ArrayXd img, double angle, int incIdxPix, int idxPixStart, int c1, int c2){
    ArrayXd projection;
    projection = ArrayXd::Zero(Nu);
    int B1,B2;
    SparseMatrix<double, RowMajor> A;
    //cout<<img<<endl;
    if (img.size()!=Lx*Ly){
        A.resize(Nu, Lx*Ly*pow(factor,2));
        A.reserve(VectorXi::Constant(Nu, 4*Ly*factor));
        B1 = Ly*factor;
        B2 = Lx*factor;
             dxPix = Sx/(Lx*factor);
     dyPix = Sy/(Ly*factor);
    
     xBoundaryPix = (ArrayXd::LinSpaced(Lx*factor + 1, 0, Lx*factor)) * dxPix + (Su - Sx)/2.;
     yBoundaryPix = (ArrayXd::LinSpaced(Ly*factor + 1, 0, Ly*factor)) * dyPix + (DSD - DSO - Sy/2.);
     xCenterPix = (ArrayXd::LinSpaced(Lx*factor, 0, Lx*factor - 1)) * dxPix + (Su - Sx + dxPix)/2.;
     yCenterPix = (ArrayXd::LinSpaced(Ly*factor, 0, Ly*factor - 1)) * dyPix + (DSD - DSO - (Sy - dyPix)/2.);
     if(idxPixStart!=0) idxPixStart *=factor;
        
    }
    else{
        A.resize(Nu, Lx*Ly);
        A.reserve(VectorXi::Constant(Nu, 4*Ly));
        B1 = Ly;
        B2 = Lx;
        dxPix = Sx/Lx;
        dyPix = Sy/Ly;
    
        xBoundaryPix = (ArrayXd::LinSpaced(Lx + 1, 0, Lx)) * dxPix + (Su - Sx)/2.;
        yBoundaryPix = (ArrayXd::LinSpaced(Ly + 1, 0, Ly)) * dyPix + (DSD - DSO - Sy/2.);
        xCenterPix = (ArrayXd::LinSpaced(Lx, 0, Lx - 1)) * dxPix + (Su - Sx + dxPix)/2.;
        yCenterPix = (ArrayXd::LinSpaced(Ly, 0, Ly - 1)) * dyPix + (DSD - DSO - (Sy - dyPix)/2.);
    }
    //cout<<B1<<"   "<<B2<<endl;
    

    double theta = angle * M_PI/180.0; // Angle in radian

    for(int indexSlice = 0; indexSlice < B1; indexSlice++){ // Loop through the image column by column
        //cout<<indexSlice<<endl;
        // Compute coordinates after rotation of the pixel boundaries of the y-slice (column) under consideration
        ArrayXd xRotImg = cos(theta) * (xBoundaryPix - xImgCenter) - sin(theta) * (yCenterPix(indexSlice) - yImgCenter) + xImgCenter;
        ArrayXd yRotImg = sin(theta) * (xBoundaryPix - xImgCenter) + cos(theta) * (yCenterPix(indexSlice) - yImgCenter) + yImgCenter;


        // Pixel boundaries of the y-slice (column) under consideration are projected onto a common axis (here the line
        // obtained by indefinitely prolonging the line segment detector)
        ArrayXd xBoundaryProj = (-ySource)*(xRotImg - xSource)/(yRotImg - ySource) + xSource;
        //cout<<"succeed"<<endl;
        /* We now turn to computing the length of overlap between each pixel projection (of the column under consideration)
         * and each pixel of the line segment detector.
         * To do that we browse the axis used to project the pixel boundaries of both the image and the detector.
         */
 
        int PixIdx_x = idxPixStart; // PixIdx_x and DetIdx_x respectively keep track of the image pixel and the detector
        int DetIdx_x = 0;           // pixel under consideration when browsing the axis used to project the boundaries.
        while (xDet(DetIdx_x) > xBoundaryProj(PixIdx_x + incIdxPix)){
            PixIdx_x = PixIdx_x + incIdxPix;
        }
  //if(indexSlice == 219) {cout<<xBoundaryProj<<endl;}
        while (xBoundaryProj(PixIdx_x) > xDet(DetIdx_x + 1)){
            DetIdx_x = DetIdx_x + 1;
            //cout<<DetIdx_x<<endl;
        }
     
        double leftBoundary_x = max(xBoundaryProj(PixIdx_x), xDet(DetIdx_x));
    
       
        // The axis onto which the pixel boundaries are projected is browsed.
        while((PixIdx_x <= (B2 -  c2)) && (PixIdx_x >= c1) && (DetIdx_x <= (Nu - 1))){
            // oX <- length of overlap between the current image and detector pixels
            double oX = (min(xBoundaryProj(PixIdx_x + incIdxPix), xDet(DetIdx_x + 1)) - leftBoundary_x)/duDet;

            // cos_alpha <- cosine of the angle between the column under consideration and
            //              the line going from the source to the current detector
            double cos_alpha = incIdxPix * (cos(theta) * ySource - sin(theta)*(xSource - (xDet(DetIdx_x) + duDet/2.)))
                               /sqrt(pow(ySource,2) + pow(xSource - (xDet(DetIdx_x) + duDet/2.),2));

            // Update the projection of the image obtained by the detector pixel under consideration
            projection.coeffRef(DetIdx_x) += oX * (dyPix/cos_alpha) * img(PixIdx_x-c1 + indexSlice*B2);
            //cout<< img(PixIdx_x-c1 + indexSlice*B2)<<endl;

            // Store the coefficient weighting the pixel value in the projection matrix A
            A.insert(DetIdx_x, (PixIdx_x - c1) + indexSlice * B2) = oX * (dyPix/cos_alpha);

            if (xDet(DetIdx_x + 1) <= xBoundaryProj(PixIdx_x + incIdxPix)){
                DetIdx_x = DetIdx_x + 1;
                leftBoundary_x = xDet(DetIdx_x);
            } else{
                PixIdx_x = PixIdx_x + incIdxPix;
                leftBoundary_x = xBoundaryProj(PixIdx_x);
            }
        }

    }
  

    A.makeCompressed();

    return make_pair(projection, A);
}



/* Function distanceDrivenProj_XSlices(ArrayXd img, double angle, int incIdxPix, int idxPixStart, int c1, int c2):
 *
 *  Perform the distance-driven projection by processing the image row-wise.
 *
 *  Inputs: img     <-  image to project
 *          angle   <-  rotation angle of the image (in degrees)
 *          incIdxPix, idxPixStart, c1, c2
 *                  <-  variables to browse the boundary projections and manage the pixel index correctly
 *                      when doing the distance-driven projection
 *
 *  Outputs:    projection  <-  distance-driven projection of the image for the considered angle
 *              A           <-  corresponding projection matrix i.e. projection = A * img
 *
 */
pair<ArrayXd, SparseMatrix<double, RowMajor>> distanceDrivenProj_XSlices(ArrayXd img, double angle, int incIdxPix, int idxPixStart, int c1, int c2){
    ArrayXd projection;
    //
    projection = ArrayXd::Zero(Nu);
    SparseMatrix<double, RowMajor> A;
    int B1,B2;
    if (img.size()!=Lx*Ly){
        A.resize(Nu, Lx*Ly*pow(factor,2));
        A.reserve(VectorXi::Constant(Nu, 4*Ly*factor));
        B1 = Ly*factor;
        B2 = Lx*factor;
             dxPix = Sx/(Lx*factor);
     dyPix = Sy/(Ly*factor);
    
     xBoundaryPix = (ArrayXd::LinSpaced(Lx*factor + 1, 0, Lx*factor)) * dxPix + (Su - Sx)/2.;
     yBoundaryPix = (ArrayXd::LinSpaced(Ly*factor + 1, 0, Ly*factor)) * dyPix + (DSD - DSO - Sy/2.);
     xCenterPix = (ArrayXd::LinSpaced(Lx*factor, 0, Lx*factor - 1)) * dxPix + (Su - Sx + dxPix)/2.;
     yCenterPix = (ArrayXd::LinSpaced(Ly*factor, 0, Ly*factor - 1)) * dyPix + (DSD - DSO - (Sy - dyPix)/2.);
        if(idxPixStart!=0) idxPixStart *=factor;
    }
    else{
        A.resize(Nu, Lx*Ly);
        A.reserve(VectorXi::Constant(Nu, 4*Ly));
        B1 = Ly;
        B2 = Lx;
        dxPix = Sx/Lx;
        dyPix = Sy/Ly;
    
        xBoundaryPix = (ArrayXd::LinSpaced(Lx + 1, 0, Lx)) * dxPix + (Su - Sx)/2.;
        yBoundaryPix = (ArrayXd::LinSpaced(Ly + 1, 0, Ly)) * dyPix + (DSD - DSO - Sy/2.);
        xCenterPix = (ArrayXd::LinSpaced(Lx, 0, Lx - 1)) * dxPix + (Su - Sx + dxPix)/2.;
        yCenterPix = (ArrayXd::LinSpaced(Ly, 0, Ly - 1)) * dyPix + (DSD - DSO - (Sy - dyPix)/2.);
    }

    double theta = angle * M_PI/180.0; //Angle in radian

    for(int indexSlice=0; indexSlice < B2; indexSlice++){ // Loop through the image row by row
        //cout<<indexSlice<<endl;
        // Compute coordinates after rotation of the pixel boundaries of the x-slice (row) under consideration
        ArrayXd xRotImg = cos(theta)*(xCenterPix(indexSlice) - xImgCenter) - sin(theta)*(yBoundaryPix-yImgCenter) + xImgCenter;
        ArrayXd yRotImg = sin(theta)*(xCenterPix(indexSlice) - xImgCenter) + cos(theta)*(yBoundaryPix-yImgCenter) + yImgCenter;

        // Pixel boundaries of the x-slice (row) under consideration are projected onto a common axis (here the line
        // obtained by indefinitely prolonging the line segment detector)
        ArrayXd yBoundaryProj = (-ySource)*(xRotImg - xSource)/(yRotImg - ySource) + xSource;

        /* We now turn to computing the length of overlap between each pixel projection (of the row under consideration)
         * and each pixel of the line segment detector.
         * To do that we browse the axis used to project the pixel boundaries of both the image and the detector.
         */

        int PixIdx_y = idxPixStart; // PixIdx_y and DetIdx_x respectively keep track of the image pixel and the detector
        int DetIdx_x = 0;           // pixel under consideration when browsing the axis used to project the boundaries.
        while (xDet(DetIdx_x) > yBoundaryProj(PixIdx_y + incIdxPix)){
            PixIdx_y = PixIdx_y + incIdxPix;
        }
        //cout<<xDet<<endl;
        //if(indexSlice == 177 && angle == sett::deg(8)) cout<<PixIdx_y<<"   "<<DetIdx_x + 1<<endl;
        while (yBoundaryProj(PixIdx_y) > xDet(DetIdx_x + 1)){
            
            DetIdx_x = DetIdx_x + 1;
           // cout<<DetIdx_x<<endl;
        }
        
        double leftBoundary_x = max(yBoundaryProj(PixIdx_y), xDet(DetIdx_x));
        
        // The axis onto which the pixel boundaries are projected is browsed.
        while((PixIdx_y <= (B1 -  c2)) && (PixIdx_y >= c1) && (DetIdx_x <= (Nu - 1))){
            // oX <- length of overlap between the current image and detector pixels
            double oX = (min(yBoundaryProj(PixIdx_y + incIdxPix), xDet(DetIdx_x + 1)) - leftBoundary_x)/duDet;

            // cos_alpha <- cosine of the angle between the row under consideration and
            //              the line going from the source to the current detector
            double cos_alpha = -incIdxPix * (cos(theta - M_PI_2)*ySource - sin(theta - M_PI_2)*(xSource - (xDet(DetIdx_x) + duDet/2.)))
                        /sqrt(pow(ySource,2) + pow(xSource - (xDet(DetIdx_x) + duDet/2.),2));

            // Update the projection of the image obtained by the detector pixel under consideration
            projection.coeffRef(DetIdx_x) += oX * (dyPix/cos_alpha) * img(indexSlice + (PixIdx_y - c1)*B2);

            // Store the coefficient weighting the pixel value in the projection matrix A
            A.insert(DetIdx_x, indexSlice + B2 * (PixIdx_y - c1)) = oX * (dyPix/cos_alpha);

            if (xDet(DetIdx_x + 1) <= yBoundaryProj(PixIdx_y + incIdxPix)){
                DetIdx_x = DetIdx_x + 1;
                leftBoundary_x = xDet(DetIdx_x);
            } else{
                PixIdx_y = PixIdx_y + incIdxPix;
                leftBoundary_x = yBoundaryProj(PixIdx_y);
            }
        }

    }

    A.makeCompressed();

    return make_pair(projection, A);
}







//////////////////////////////////////////////////////////hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh//////////////////////////////////////////////////////////////////////////////////////////////
////////////////////hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh/////////////////////
pair<ArrayXd,SparseMatrix<double, RowMajor>> distanceDrivenProjection_h(ArrayXd img, int iview){
    double angle = deg(iview);
    cout<<iview<<endl;
    if((-45 <= angle) && (angle <=45)){

        return distanceDrivenProj_YSlices_h(img, angle);


    } else if((45 < angle) && (angle < 135)){

        return distanceDrivenProj_XSlices_h(img, angle, -1, Ly*factor, 1, 0);

    } else if(((135 <= angle) && (angle<=180)) || ((-180 <= angle) && (angle <= -135))){

        return distanceDrivenProj_YSlices_h(img, angle, -1, Lx*factor, 1, 0);

    } else{

        return distanceDrivenProj_XSlices_h(img, angle);

    }
}



/* Function distanceDrivenProj_YSlices_h(ArrayXd img, double angle, int incIdxPix, int idxPixStart, int c1, int c2):
 *
 *  Perform the distance-driven projection by processing the image column-wise.
 *
 *  Inputs: img     <-  image to project
 *          angle   <-  rotation angle of the image (in degrees)
 *          incIdxPix, idxPixStart, c1, c2
 *                  <-  variables to browse the boundary projections and manage the pixel index correctly
 *                      when doing the distance-driven projection
 *
 *  Outputs:    projection  <-  distance-driven projection of the image for the considered angle
 *              A           <-  corresponding projection matrix i.e. projection = A * img
 *
 */
pair<ArrayXd, SparseMatrix<double, RowMajor>> distanceDrivenProj_YSlices_h(ArrayXd img, double angle, int incIdxPix, int idxPixStart, int c1, int c2){
    ArrayXd projection;
    projection = ArrayXd::Zero(Nu);
    //cout<<Mxh<<"   "<<Myh<<endl;
    SparseMatrix<double, RowMajor> A(Nu, Mxh*Myh);
    dxPix = Sx/(Lx*factor);
    dyPix = Sy/(Ly*factor);
    
    
     xBoundaryPix = (ArrayXd::LinSpaced(Lx*factor + 1, 0, Lx*factor)) * dxPix + (Su - Sx)/2.;
     yBoundaryPix = (ArrayXd::LinSpaced(Ly*factor + 1, 0, Ly*factor)) * dyPix + (DSD - DSO - Sy/2.);
     xCenterPix = (ArrayXd::LinSpaced(Lx*factor, 0, Lx*factor - 1)) * dxPix + (Su - Sx + dxPix)/2.;
     yCenterPix = (ArrayXd::LinSpaced(Ly*factor, 0, Ly*factor - 1)) * dyPix + (DSD - DSO - (Sy - dyPix)/2.);
    //A.reserve(VectorXi::Constant(Nu, 4*Ly));

    double theta = angle * M_PI/180.0; // Angle in radian
    int count = -1;
    
    for(int indexSlice = upperlefth[1]; indexSlice <= lowerrighth[1]; indexSlice++){ // Loop through the image column by column
        //cout<<indexSlice<<endl;
        // Compute coordinates after rotation of the pixel boundaries of the y-slice (column) under consideration
        //cout<<yCenterPix(indexSlice)<<endl;
        ArrayXd xRotImg = cos(theta) * (xBoundaryPix - xImgCenter) - sin(theta) * (yCenterPix(indexSlice) - yImgCenter) + xImgCenter;
        ArrayXd yRotImg = sin(theta) * (xBoundaryPix - xImgCenter) + cos(theta) * (yCenterPix(indexSlice) - yImgCenter) + yImgCenter;
        
        
        // Pixel boundaries of the y-slice (column) under consideration are projected onto a common axis (here the line
        // obtained by indefinitely prolonging the line segment detector)
        ArrayXd xBoundaryProj = (-ySource)*(xRotImg - xSource)/(yRotImg - ySource) + xSource;

        /* We now turn to computing the length of overlap between each pixel projection (of the column under consideration)
         * and each pixel of the line segment detector.
         * To do that we browse the axis used to project the pixel boundaries of both the image and the detector.
         */

        int PixIdx_x = idxPixStart; // PixIdx_x and DetIdx_x respectively keep track of the image pixel and the detector
        int DetIdx_x = 0;           // pixel under consideration when browsing the axis used to project the boundaries.
        while (xDet(DetIdx_x) > xBoundaryProj(PixIdx_x + incIdxPix)){
            PixIdx_x = PixIdx_x + incIdxPix;
        }

        while (xBoundaryProj(PixIdx_x) > xDet(DetIdx_x + 1)){
            DetIdx_x = DetIdx_x + 1;
        }

        double leftBoundary_x = max(xBoundaryProj(PixIdx_x), xDet(DetIdx_x));

        // The axis onto which the pixel boundaries are projected is browsed.
        int indicator = -1;
        while((PixIdx_x <= (Lx*factor -  c2)) && (PixIdx_x >= c1) && (DetIdx_x <= (Nu - 1))){
            // oX <- length of overlap between the current image and detector pixels
            //cout<<"indexSlice: "<<indexSlice<<endl;
            //cout<<"DetIdx_x: "<<DetIdx_x<<endl;
            //cout<<"PixIdx_x: "<<PixIdx_x<<endl;
            double oX = (min(xBoundaryProj(PixIdx_x + incIdxPix), xDet(DetIdx_x + 1)) - leftBoundary_x)/duDet;

            // cos_alpha <- cosine of the angle between the column under consideration and
            //              the line going from the source to the current detector
            double cos_alpha = incIdxPix * (cos(theta) * ySource - sin(theta)*(xSource - (xDet(DetIdx_x) + duDet/2.)))
                               /sqrt(pow(ySource,2) + pow(xSource - (xDet(DetIdx_x) + duDet/2.),2));

            // Update the projection of the image obtained by the detector pixel under consideration
                projection.coeffRef(DetIdx_x) += oX * (dyPix/cos_alpha) * img(PixIdx_x-c1 + indexSlice*Lx);
            if((PixIdx_x - c1) >= upperlefth[0] && (PixIdx_x - c1)<= lowerrighth[0]){
                
            
            projection.coeffRef(DetIdx_x) += oX * (dyPix/cos_alpha) * img(PixIdx_x-c1 + indexSlice*Lx);
            

            // Store the coefficient weighting the pixel value in the projection matrix A
            //cout<<DetIdx_x<<", "<<(PixIdx_x - c1) + indexSlice * Lx<<endl;
                
                if(indicator != (PixIdx_x - c1) + indexSlice * Lx*factor){count++;}
                indicator = (PixIdx_x - c1) + indexSlice * Lx*factor;
                //cout<<count<<", "<<indicator<<endl;
                //if(count == 1 && indicator ==19) cout<<match(indicator)<<endl;
                
                
                /// The reason why it encount errors here sometime; The map "match", which corresponds pixel in different resolution stage together, is incomplete if the high resolution area is too   big for the detector. Sometimes the pixels, which didn't appear in the first angel will be assigned to 0. Thus if they apear in another angle of view, the programm just doesn't know where it belongs to. Then It will give you the error: can not insert twice. Currently I think no need to fix it, since if it is that big, the whole reconstruction doesn't make sense any more.
                if(angle == 0){match(indicator) = count; A.insert(DetIdx_x, count) = oX * (dyPix/cos_alpha);}
                else{A.insert(DetIdx_x, match(indicator)) = oX * (dyPix/cos_alpha);}
                
                
            //cout<<"yige"<<endl;
            }

            if (xDet(DetIdx_x + 1) <= xBoundaryProj(PixIdx_x + incIdxPix)){
                DetIdx_x = DetIdx_x + 1;
                leftBoundary_x = xDet(DetIdx_x);
            } 
            else{
                PixIdx_x = PixIdx_x + incIdxPix;
                leftBoundary_x = xBoundaryProj(PixIdx_x);
            }
        }

    }
  

    A.makeCompressed();

    return make_pair(projection, A);
}



/* Function distanceDrivenProj_XSlices_h(ArrayXd img, double angle, int incIdxPix, int idxPixStart, int c1, int c2):
 *
 *  Perform the distance-driven projection by processing the image row-wise.
 *
 *  Inputs: img     <-  image to project
 *          angle   <-  rotation angle of the image (in degrees)
 *          incIdxPix, idxPixStart, c1, c2
 *                  <-  variables to browse the boundary projections and manage the pixel index correctly
 *                      when doing the distance-driven projection
 *
 *  Outputs:    projection  <-  distance-driven projection of the image for the considered angle
 *              A           <-  corresponding projection matrix i.e. projection = A * img
 *
 */
pair<ArrayXd, SparseMatrix<double, RowMajor>> distanceDrivenProj_XSlices_h(ArrayXd img, double angle, int incIdxPix, int idxPixStart, int c1, int c2){
    ArrayXd projection;
    projection = ArrayXd::Zero(Nu);
    SparseMatrix<double, RowMajor> A(Nu, Mxh*Myh);
    //A.reserve(VectorXi::Constant(Nu, 4*Lx));
        dxPix = Sx/(Lx*factor);
    dyPix = Sy/(Ly*factor);
    
     xBoundaryPix = (ArrayXd::LinSpaced(Lx*factor + 1, 0, Lx*factor)) * dxPix + (Su - Sx)/2.;
     yBoundaryPix = (ArrayXd::LinSpaced(Ly*factor + 1, 0, Ly*factor)) * dyPix + (DSD - DSO - Sy/2.);
     xCenterPix = (ArrayXd::LinSpaced(Lx*factor, 0, Lx*factor - 1)) * dxPix + (Su - Sx + dxPix)/2.;
     yCenterPix = (ArrayXd::LinSpaced(Ly*factor, 0, Ly*factor - 1)) * dyPix + (DSD - DSO - (Sy - dyPix)/2.);

    double theta = angle * M_PI/180.0; //Angle in radian
    int count = -1;
    for(int indexSlice = upperlefth[0]; indexSlice <= lowerrighth[0]; indexSlice++){ // Loop through the image row by row

        // Compute coordinates after rotation of the pixel boundaries of the x-slice (row) under consideration
        ArrayXd xRotImg = cos(theta)*(xCenterPix(indexSlice) - xImgCenter) - sin(theta)*(yBoundaryPix-yImgCenter) + xImgCenter;
        ArrayXd yRotImg = sin(theta)*(xCenterPix(indexSlice) - xImgCenter) + cos(theta)*(yBoundaryPix-yImgCenter) + yImgCenter;

        // Pixel boundaries of the x-slice (row) under consideration are projected onto a common axis (here the line
        // obtained by indefinitely prolonging the line segment detector)
        ArrayXd yBoundaryProj = (-ySource)*(xRotImg - xSource)/(yRotImg - ySource) + xSource;

        /* We now turn to computing the length of overlap between each pixel projection (of the row under consideration)
         * and each pixel of the line segment detector.
         * To do that we browse the axis used to project the pixel boundaries of both the image and the detector.
         */

        int PixIdx_y = idxPixStart; // PixIdx_y and DetIdx_x respectively keep track of the image pixel and the detector
        int DetIdx_x = 0;           // pixel under consideration when browsing the axis used to project the boundaries.
        while (xDet(DetIdx_x) > yBoundaryProj(PixIdx_y + incIdxPix)){
            PixIdx_y = PixIdx_y + incIdxPix;
        }
        //if(indexSlice == 177 && angle == deg(8)) cout<<PixIdx_y<<"   "<<DetIdx_x + 1<<endl;
        while (yBoundaryProj(PixIdx_y) > xDet(DetIdx_x + 1)){
            DetIdx_x = DetIdx_x + 1;
        }

        double leftBoundary_x = max(yBoundaryProj(PixIdx_y), xDet(DetIdx_x));

        // The axis onto which the pixel boundaries are projected is browsed.
        int indicator = -1;
        while((PixIdx_y <= (Ly*factor -  c2)) && (PixIdx_y >= c1) && (DetIdx_x <= (Nu - 1))){
            // oX <- length of overlap between the current image and detector pixels
            double oX = (min(yBoundaryProj(PixIdx_y + incIdxPix), xDet(DetIdx_x + 1)) - leftBoundary_x)/duDet;

            // cos_alpha <- cosine of the angle between the row under consideration and
            //              the line going from the source to the current detector
            double cos_alpha = -incIdxPix * (cos(theta - M_PI_2)*ySource - sin(theta - M_PI_2)*(xSource - (xDet(DetIdx_x) + duDet/2.)))
                        /sqrt(pow(ySource,2) + pow(xSource - (xDet(DetIdx_x) + duDet/2.),2));

            // Update the projection of the image obtained by the detector pixel under consideration
            projection.coeffRef(DetIdx_x) += oX * (dyPix/cos_alpha) * img(indexSlice + (PixIdx_y - c1)*Lx);
            //cout<<"indexSlice: "<<indexSlice<<endl;
            //cout<<"DetIdx_x: "<<DetIdx_x<<endl;
            //cout<<"PixIdx_y: "<<PixIdx_y<<endl;

            // Store the coefficient weighting the pixel value in the projection matrix A
            
            if((PixIdx_y - c1) >= upperlefth[1] && (PixIdx_y - c1) <= lowerrighth[1]){
                
            
            //projection.coeffRef(DetIdx_x) += oX * (dyPix/cos_alpha) * img(PixIdx_x-c1 + indexSlice*Lx);

            // Store the coefficient weighting the pixel value in the projection matrix A
       
                
            //if(indicator != indexSlice + Lx * (PixIdx_y - c1)){count++;}
            indicator = indexSlice + Lx*factor * (PixIdx_y - c1);
            //cout<<count<<", "<<indicator<<endl;
            
            A.insert(DetIdx_x,  match(indicator)) = oX * (dyPix/cos_alpha);
            //cout<<"yige"<<endl;
            }
            //A.insert(DetIdx_x, indexSlice + Lx * (PixIdx_y - c1)) = oX * (dyPix/cos_alpha);

            if (xDet(DetIdx_x + 1) <= yBoundaryProj(PixIdx_y + incIdxPix)){
                DetIdx_x = DetIdx_x + 1;
                leftBoundary_x = xDet(DetIdx_x);
            } else{
                PixIdx_y = PixIdx_y + incIdxPix;
                leftBoundary_x = yBoundaryProj(PixIdx_y);
            }
        }

    }

    A.makeCompressed();

    return make_pair(projection, A);
}












