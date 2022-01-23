#ifndef SETTINGS_H
#define SETTINGS_H

//  File Name : settings.h    Purpose : Settings of the 2D fan-beam geometry
#include "../Eigen/Core" //Eigen dependencies

namespace sett {
    // Object selection
    const string objSelection ="SheppLogan256_UltraGrad";//"ReversedSheppLogan256_Grad";//"ReversedRanGradSqauresTer256";//"SquaresTer256";////"03-concrete_block";//"04-steel_hole_plate";//"01-steel_step_cylinder";//"02-small_steel_cast_part";//ModifiedSheppLogan512

    // Image size (in pixels)
    const int Lx =256; //400; //64;
    const int Ly =256; //400; //64;

    //Lp norm 
    const double p = 0.5; 

    // Image size (in mm)
    const double Sx = 200.;
    const double Sy = 200.;

    // Number of pixels of the line segment detector
    const int Nu = 450;//1728;//450; //128;

    // Detector size (in mm)
    const double Su = 500;

    // Distance Source to Detector
    const double DSD = 800;

    // Distance Source to Center of the image
    const double DSO =  400;

    // Pixel detector size (in mm)
    const double duDet = Su/Nu;

    // Angle setting
    //const Eigen::ArrayXd deg = (Eigen::ArrayXd(10) << Eigen::ArrayXd::LinSpaced(5,0,144), Eigen::ArrayXd::LinSpaced(5,198,342)-360).finished();
    const Eigen::ArrayXd deg = (Eigen::ArrayXd(12) << Eigen::ArrayXd::LinSpaced(6,0,152), Eigen::ArrayXd::LinSpaced(6,195,344)-360).finished();
    //const Eigen::ArrayXd deg = (Eigen::ArrayXd(14) << Eigen::ArrayXd::LinSpaced(7,0,156), Eigen::ArrayXd::LinSpaced(7,193,346)-360).finished();
    //const Eigen::ArrayXd deg = (Eigen::ArrayXd(16) << Eigen::ArrayXd::LinSpaced(8,0,157.5), Eigen::ArrayXd::LinSpaced(8,191.25,348.75)-360).finished();
    //const Eigen::ArrayXd deg = (Eigen::ArrayXd(18) << Eigen::ArrayXd::LinSpaced(9,0,160), Eigen::ArrayXd::LinSpaced(9,190,350)-360).finished();
    //const Eigen::ArrayXd deg = (Eigen::ArrayXd(20) << Eigen::ArrayXd::LinSpaced(10,0,162), Eigen::ArrayXd::LinSpaced(10,189,351)-360).finished();
    //const Eigen::ArrayXd deg = (Eigen::ArrayXd(22) << Eigen::ArrayXd::LinSpaced(11,0,163.64), Eigen::ArrayXd::LinSpaced(11,188.18,351)-360).finished();
    //const Eigen::ArrayXd deg = (Eigen::ArrayXd(24) << Eigen::ArrayXd::LinSpaced(12,0,165), Eigen::ArrayXd::LinSpaced(12,187.5,352.5)-360).finished();
    //const Eigen::ArrayXd deg = (Eigen::ArrayXd(28) << Eigen::ArrayXd::LinSpaced(14,0,167.1), Eigen::ArrayXd::LinSpaced(14,186.5,353.5)-360).finished();
    //const Eigen::ArrayXd deg = (Eigen::ArrayXd(30) << Eigen::ArrayXd::LinSpaced(15,0,168), Eigen::ArrayXd::LinSpaced(15,186,354)-360).finished();
    //const Eigen::ArrayXd deg = (Eigen::ArrayXd(32) << Eigen::ArrayXd::LinSpaced(16,0,168.75), Eigen::ArrayXd::LinSpaced(16,185.625,354.375)-360).finished();
    //const Eigen::ArrayXd deg = (Eigen::ArrayXd(36) << Eigen::ArrayXd::LinSpaced(18,0,170), Eigen::ArrayXd::LinSpaced(18,185,355)-360).finished();
    //const Eigen::ArrayXd deg = (Eigen::ArrayXd(40) << Eigen::ArrayXd::LinSpaced(20,0,171), Eigen::ArrayXd::LinSpaced(20,184.5,355.5)-360).finished();
    //const Eigen::ArrayXd deg = (Eigen::ArrayXd(50) << Eigen::ArrayXd::LinSpaced(25,0,172.8), Eigen::ArrayXd::LinSpaced(25,183.6,356.4)-360).finished()
    //const Eigen::ArrayXd deg = (Eigen::ArrayXd(60) << Eigen::ArrayXd::LinSpaced(30,0,174), Eigen::ArrayXd::LinSpaced(30,183,357)-360).finished();
    //const Eigen::ArrayXd deg = (Eigen::ArrayXd(120) << Eigen::ArrayXd::LinSpaced(60,0,177), Eigen::ArrayXd::LinSpaced(60,181.5,358.5)-360).finished();
	//const Eigen::ArrayXd deg = (Eigen::ArrayXd(180) << Eigen::ArrayXd::LinSpaced(90, 0, 178), Eigen::ArrayXd::LinSpaced(90, 180.5, 359) - 360).finished();
    //const Eigen::ArrayXd deg = (Eigen::ArrayXd(360) << Eigen::ArrayXd::LinSpaced(180,0,179), Eigen::ArrayXd::LinSpaced(180,180,359)-360).finished();
    //const Eigen::ArrayXd deg = (Eigen::ArrayXd(1440) << Eigen::ArrayXd::LinSpaced(720,0,179.75), Eigen::ArrayXd::LinSpaced(720,180,359.75)-360).finished();

    int nProj = deg.size();
    
    //Multi-Resolution Settings ()
    int factor = 1;

    int upperleft[2] = {55-1,55-1};
    int lowerright[2] = {56-1,56-1};
    int Mx = lowerright[0] - upperleft[0] + 1;
    int My = lowerright[1] - upperleft[1] + 1;
    int upperlefth[2] = {(upperleft[0])*factor,(upperleft[1])*factor};
    int lowerrighth[2] = {(lowerright[0]+1)*factor-1,(lowerright[1]+1)*factor-1};
    const int Mxh = lowerrighth[0] - upperlefth[0] + 1;
    const int Myh = lowerrighth[1] - upperlefth[1] + 1;
    int degree = 2;
    Eigen::ArrayXd match(Lx*Ly*factor*factor);
    int Nr = 8,Nc = 8,Margin = 16;
        

    // Geometry calculation
    // Pixel boundaries of the detector
    const Eigen::ArrayXd xDet = (Eigen::ArrayXd::LinSpaced(Nu + 1, 0, Nu)) * duDet;

    // Source coordinates
    const double xSource = Su/2.;
    const double ySource = DSD;

    // Coordinates of the image center
    const double xImgCenter = Su/2.;
    const double yImgCenter = DSD - DSO;
    
    // IRCD Settings
    //double SigmaZ=0.01;
    //bool AddNoiseToProj = false;
    //bool Jeffery = false;
    int max_num = 10000, display = 100; // Maximum iteration number and display interval
    bool Ran_order = false; // Randomize update order of pixels
    bool intermediate = true; // save intermediate results
    bool visulize_original = true; // generate a png for original raw image
}
#endif
