
/**
 * @file matrix.cpp
 * @brief Implementation file for the matrix class.
 */

#include "matrix.h"
#include <Accelerate/Accelerate.h>
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
mat::mat(int &N, int &M) {

    // Initialize matrix with size NxM
    data = vector<double>(N*M, 0.);

    this->N = N;
    this->M = M;
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
    if (i*M+j > data.size()){
        throw std::runtime_error("The matrix has too few entries.");
    }

    // Return entry.
    return data[i*M+j];
}

double& mat::getReference(int i, int j) {

    // Check that entry is inside matrix bounds.
    if (i*M+j > data.size()){
        throw std::runtime_error("The matrix has too few entries.");
    }

    // Return entry.
    return data[i*M+j];
}


void mat::setEntry(int i, int j, double value) {

    // Check that entry is inside matrix bounds.
    if (i*M+j > data.size()){
        throw std::runtime_error("The matrix has too few entries.");
    }

    // Set entry value
    data[i*M+j] = value;

}

int mat::nRows() {
    return N;
}

int mat::nCols(){
    return M;
}

void mat::getRow(int i,vector<double> &res) {

    // Check that the dimensions match
    if (res.size() != M){
        throw runtime_error("Vector needs to be the same size as the number of columns");
    }

    // Allocate row
    for (int k=0; k<res.size(); k++){
        res[k] = this->data[i*N+k];
    }

}

void mat::getCol(int j, vector<double> &res) {

    //Check that dimension match
    if (res.size() != N){
        throw runtime_error("Vector needs to be the same size as the number of rows");
    }

    // Allocate Column
    for (int k=0; k<res.size(); k++){
        res[k] = this->data[k*N+j];
    }

}


void mat::matmul(mat &xmat, mat &res) {

    // Check that dimensions match
    if(M != xmat.nRows()){
        throw std::runtime_error("Matrix multiplication can only take the format MxN * N*K");
    }

    // perform matrix multiplication
    for (int i=0; i<N; i++){
        for (int j=0; j<xmat.nCols(); j++){

            // Perform vector sum product
            double &x = res.getReference(i, j);
            for (int k=0; k<M; k++){
                x += this->getEntry(i,k)*xmat.getEntry(k,j);
            }

        }
    }

}