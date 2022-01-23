# README of SNUV_TOMO#

## Project tree and description ##
| File/folder   | Description                    |
| ------------- | ------------------------------ |
| **Eigen**     | [C++ template library for linear algebra](http://eigen.tuxfamily.org/index.php?title=Main_Page)   |
| **include**   | Folder comprising self-made functions to save (resp. load) a sparse matrix or an array into (resp. from) a binary file; Mp.h contains the main IRCD routine.|
| **objects**   | Contains different synthetic image to test reconstruction|
| settings.h   | Specifies the studied object and some important quantities for the forward projection|
| projection.cpp  | Perform the distance-driven projection for the 2D fan-beam geometry |
| mainSNUV_IRCD.cpp   | Perform reconstruction using Iterative Reweighted Coordinate Descent (IRCD) Algorithm with Smoothed-NUV priors |
| mainASDPOCS.cpp   | Perform ASDPOCS reconstruction|
| transformer.cpp   | Transform normal images into binary files|
| CMakeLists.txt   | CMake file used to generate a Makefile |

## Before running the applications ##
* Create 4 folders **build**, **projections**, **reconstructions** and **reconstructions/intermediate_results**.

* Open a terminal and check the version of CMake installed on your computer with the command *cmake -version*. If the version is below 3.5 then modify the first line of **CMakeLists.txt**.

## Example ##
The **objects** folder contains different, for example, the binary file **SheppLogan256_UltraGrad** corresponds to a Shepp-Logan phantom with gradient of dimension 256 x 256, flattened column-wise (if a pixel has coordinates (m,n) with 0 <= m,n <= 255 then it has index m + 255*n in the flattened image).

As a first example you can project this phantom and reconstruct it from the projections using **SNUV_IRCD** or **ASDPOCS** reconstruction methods.

First you have to compile the different applications listed in the CMake file (**CMakeLists.txt**):

1. Open a Terminal window

2. Access the **build** folder: *cd path/to/SNUV_TOMO/build/*

3. Run the command *cmake ..* to create **Makefile** from **CMakeLists.txt**.

4. Run the command *make* to compile all the applications at once. If you want to compile just one application, say PROJECTION, then run the command *make PROJECTION*.

After step 4, you should have 4 applications in the **build** folder (the files **PROJECTION**, **ASDPOCS**, **SNUV_IRCD**, **TRANSFORMER**).

Now you can project the phantom and use the projection to reconstruct the image:

1. Open a Terminal window.

2. Access the **build** folder: *cd path/to/ct-2d/build/*.

3. If the input image is not binary, please first run *./TRANSFORMER* to transform it into binary files.

4. Run the command *./PROJECTION* to perform the distance-driven projection for the different angles specified in the header **settings.h**. The resulting projection vector and projection matrix are saved in 2 separate binary files located in the **projections** folder.

5. The maximal iteration number of IRCD, and whether the intermediate results should be saved, can be specified in **settings.h**.

6. Run the commands *./SNUV_IRCD* and *./ASDPOCS*to reconstruct the image from the projection vector using SNUV-IRCD and TV-ASDPOCS, respectively. The reconstructions are saved in binary&png files located in the **reconstructions** folder.


