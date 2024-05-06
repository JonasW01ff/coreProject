#pragma once
#ifndef LIBCORE_LIBRARY_H
#define LIBCORE_LIBRARY_H

#include <vector>

using namespace std;

class mat {
private:
    // Matrix object
    vector<vector<double> > data;

public:

    // Constructor
    mat();

    // Overloaded Constructor
    mat(int N, int M);

    // Destructor
    ~mat();

    // Get index
    double getEntry(int i, int j);


};


#endif //LIBCORE_LIBRARY_H
