#ifndef CT_2D_CONVERT_ARRAYXD_CSV_H
#define CT_2D_CONVERT_ARRAYXD_CSV_H

#include "../Eigen/Core"
#include <fstream>

Eigen::ArrayXd loadCSV(std::string dir, int size) {

    std::ifstream in(dir);

    std::string line;

    int col = 0;

    Eigen::ArrayXd res(size);

    if (in.is_open()) {

        std::getline(in, line);
        char *ptr = (char *) line.c_str();
        int len = line.length();

        char *start = ptr;
        for (int i = 0; i < len; i++) {

            if (ptr[i] == ',') {
                res(col) = atof(start);
                col++;
                start = ptr + i + 1;
            }
        }
        res(col) = atof(start);
        in.close();
    }
    return res;
}

void saveToCSV(Eigen::ArrayXd array, std::string name){
    Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", ", ");
    std::ofstream file (name);
    file << array.format(CSVFormat);
    file.close();
}

#endif //CT_2D_CONVERT_ARRAYXD_CSV_H
