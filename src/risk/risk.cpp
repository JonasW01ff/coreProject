//
// Created by Jonas Wolff on 26/02/2024.
//

#include "risk.h"
#include "../statistics/distribution_binomial.h"

#include <algorithm>
#include <vector>
#include <math.h>

using namespace std;


// Constructor
risk::risk(){

}

// Overloaded Constructor
risk::risk(vector<double> pLosses){
    loadLosses(pLosses);
}


// Destructor
risk::~risk(){

}


// loadLosses function
void risk::loadLosses(vector<double> pLosses){

    // Check for empty vector
    if (pLosses.size() == 0){
        throw runtime_error("at least 1 loss has to be given.");
    }


    std::sort(begin(pLosses), end(pLosses), std::greater<double>{});
    losses = pLosses;
}

// Calculate VaR
double risk::calcEmpiricalVaR(double alpha) {

    if (losses.size() == 0){
        throw runtime_error("looses need to be loaded before calculating empirical VaR");
    }

    double q_float = losses.size();
    q_float = q_float*(1-alpha);

    // +1 is because we rounded down and -1 because we count from zero.
    int q_idx = floor(q_float);
    return losses[q_idx];
}

// Calculate exact confidence intervals for Empirical VaR
vector<double> risk::calcCiEmpiricalVaR(double alpha, double probability){

    // Check that the losses variable has been initialized
    if (losses.size() == 0){
        throw runtime_error("No losses has been given to VaR calculations.");
    }

    // Find probability that each observation is lower than the VaR
    double prop;
    double n = losses.size();
    int loweridx, upperidx;
    distribution_binomial quantile_dist(n, 1.-alpha);
    for (int i=0; i<n; i++){

        // Probability that VaR is less than i+1 biggest loss
        prop = quantile_dist.cdf(i);

        // Find index pais giving less than (1-p)/2 probability to VaR falling outside on both
        // ends of the interval
        if (i == 0 && prop > (1.0-probability)/2.0){
            throw runtime_error("CI probability or alpha too high to find upper bound.");
        }
        if( prop <= (1.0-probability)/2.0){
            upperidx = i;
        } else if (1.0-prop <= (1.0-probability)/2.0) {
            loweridx = i;
            // Choose the pair that give lowest probability to VaR being in the interval
            break;
        }
    }

    // Creat confidence interval
    vector<double> CI;
    CI.push_back(losses[loweridx]);
    CI.push_back(losses[upperidx]);

    return CI;

}