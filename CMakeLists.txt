cmake_minimum_required(VERSION 3.0)
project(SNUV_IRCD)
find_package(OpenCV)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -w -O3 -lopencv_imgcodecs -lopencv_core -lopencv_imgproc -lopencv_highgui")

include_directories("Eigen" "include")

set(SOURCE_FILES settings.h mainSNUV_IRCD.cpp mainASDPOCS.cpp transformer.cpp projection.cpp)

add_executable(SNUV_IRCD
        settings.h
        include/RW_SparseMat.h
        include/RW_ArrayXd.h
        include/MP.h
        mainSNUV_IRCD.cpp)
        
add_executable(ASDPOCS
        settings.h
        include/RW_SparseMat.h
        include/RW_ArrayXd.h
        mainASDPOCS.cpp)
        
#add_executable(NUV_EM
#        settings.h
#        include/RW_SparseMat.h
#        include/RW_ArrayXd.h
#        mainGaMP.cpp)


add_executable(PROJECTION
        settings.h
        include/RW_SparseMat.h
        include/RW_ArrayXd.h
        projection.cpp)

add_executable(TRANSFORMER
        include/RW_ArrayXd.h
        transformer.cpp)

