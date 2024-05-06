
/**
 * @file matrix.cpp
 * @brief Implementation file for the matrix class.
 */

#include "matrix.h"
#include <iostream>

/**
 * @brief Default constructor for the matrix class.
 */
mat::mat() {

};

/**
 * @brief Overloaded constructor for the matrix class.
 * @param N The number of rows in the matrix.
 * @param M The number of columns in the matrix.
 */
mat::mat(int N, int M) {

    // Initialize matrix with size NxM
    data.resize(N);
    for (int i=0; i<N; i++){
        data[i].resize(M);
        for (int j=0; j<M; j++){
            // Set all entries to 0
            data[i][j] = 0;
        }
    }
};

/**
 * @brief Destructor for the matrix class.
 */
mat::~mat() {

};

/**
 * @brief Get the value of a specific entry in the matrix.
 * @param i The row index of the entry.
 * @param j The column index of the entry.
 * @return The value of the entry at position (i, j).
 * @throws std::runtime_error if the entry is outside the matrix bounds.
 */
double mat::getEntry(int i, int j) {

    // Check that entry is inside matrix bounds.
    if (i >= data.size()){
        throw std::runtime_error("The matrix has too few rows.");
    } else if (j >= data[i].size()){
        throw std::runtime_error("The matrix has too few columns.");
    }

    // Return entry.
    return data[i][j];
}