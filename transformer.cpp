#include<iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "../Eigen/Core" //Eigen dependencies
#include "include/RW_ArrayXd.h"
using namespace std;

int main(){
    // Simple function to transform normal images to Eigen arrays.
    cout<<"Please put the image in the folder objects and give its name in the following:"<<endl;
    string name;
    cin >> name;
    cv::Mat img = cv::imread("../objects/" + name, 0);
    if(img.empty())
    {
        cout << "Could not read the image: " << endl;
        return 1;
    }
    int R = img.rows;
    int C = img.cols;
    //imshow("Display window", img);
    Eigen::ArrayXd IMG = Eigen::ArrayXd::Zero(R*C);
    int j = 0;
    for(int r = 0; r < R; r++){
        for(int c = 0; c < C; c++){
            j = r + R * c;
            IMG(j) = double(img.at<uint8_t>(r,c))/255;
            //cout<<IMG(j)<<endl;
            // M.at<float>(r,c) = (*img)(j)*255;
        }
    }
    cout<<IMG.maxCoeff()<<endl;
    size_t lastindex = name.find_last_of("."); 
    string rawname = name.substr(0, lastindex); 
    saveArrayXd(IMG, "../objects/" + rawname);
    cout<<"Done! Now your have a new file ../objects/" + rawname<<endl;

    return 0;
}