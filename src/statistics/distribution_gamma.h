#include <vector>

#pragma once
#ifndef LIBRARYCORE_DISTRIBUTION_GAMMA_H
#define LIBRARYCORE_DISTRIBUTION_GAMMA_H

using namespace std;

// class object for the gamma distribution function
class distribution_gamma {
private:

    // shape parameter
    double alpha;

    // rate parameter
    double beta;

public:
    
        // Constructor
        distribution_gamma(double palpha, double pbeta);

        // Constructor Overload
        distribution_gamma(vector<double> bins, vector<int> freqs);

    
        // Destructor
        ~distribution_gamma();
    
        // PDF
        double pdf(double &x);
    
        // CDF
        double cdf(double &x);

        // Get alpha
        double getAlpha(){
            return alpha;
        };

        // Get beta
        double getBeta(){
            return beta;
        };
    
    
    };

#endif //LIBRARYCORE_DISTRIBUTION_GAMMA_H