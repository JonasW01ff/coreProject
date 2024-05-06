//
// Created by Jonas Wolff on 28/02/2024.
//

#include <stdexcept>

#include "distribution_normal.h"
#include "math.h"

using namespace std;

// Constructor
distribution_normal::distribution_normal(double pmu, double psigma){

    // Check for degenerate distribution
    if (psigma <= 0){
        throw runtime_error("The normal distribution can not have variance <= 0.");
    }

    // Save parameters
    mu = pmu;
    sigma = psigma;

}

// Destructor
distribution_normal::~distribution_normal() {

}

// PDF, take reference to x for better performance
double distribution_normal::pdf(double &x) {

    // Normalise x
    double x_std;
    x_std = (x - mu)/sigma;

    // Normalization constant
    double c_norm;
    c_norm = 1/(sigma*sqrt(2*M_PI));

    // Return density
    return c_norm*exp(-0.5*pow(x_std,2));

}

// CDF, take reference to x for better performance
double distribution_normal::cdf(double &x) {

    // Normalise x
    double x_std;
    x_std = (x - mu)/sigma;

    // Return probability
    return 0.5*(1+erf(x_std/sqrt(2)));

}

