#pragma once
#ifndef LIBCORE_LIBRARY_H
#define LIBCORE_LIBRARY_H

#include <vector>

using namespace std;

class mat {
private:
    // Matrix object, it take the form row1, row2, ...
    vector<double > data;

    // Number of rows
    int N;

    // Number of columns
    int M;

public:

    // Constructor
    mat();

    // Overloaded Constructor
    mat(int &N, int &M);

    // Destructor
    ~mat();

    // Get Entry
    double getEntry(int i, int j);

    // Get entry reference
    double& getReference(int i, int j);

    // Set Entry
    void setEntry(int i, int j, double value);

    // Get number of rows
    int nRows();

    // Get number of columns
    int nCols();

    // Get Row
    void getRow(int i, vector<double> &res);

    // Get Column
    void getCol(int j, vector<double> &res);

    // Matrix multiplication
    void matmul(mat &xmat, mat &res);


};


#endif //LIBCORE_LIBRARY_H
