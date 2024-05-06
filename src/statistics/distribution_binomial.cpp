//
// Created by Jonas Wolff on 28/02/2024.
//

#include <stdexcept>
#include "math.h"

#include "distribution_binomial.h"
#include "distribution_normal.h"
#include "../functions/functions_math.h"

using namespace std;

// Constructor
/**
 * @brief Constructs a binomial distribution object with the given parameters.
 * 
 * @param pn The number of trials in the binomial distribution.
 * @param pp The probability of success in each trial.
 * 
 * @throws std::runtime_error if pn is less than or equal to 0.
 */
distribution_binomial::distribution_binomial(int pn, double pp){

    // Check sanity of parameters
    if (pn <= 0){
        throw runtime_error("The binomial distribution can not have negative trials");
    }


    // Save parameters
    n = pn;
    p = pp;

};

// Destructor
distribution_binomial::~distribution_binomial(){

};


/**
 * Calculates the probability density function (PDF) of the binomial distribution for a given value of k.
 *
 * @param k The number of successes.
 * @return The probability density at k.
 * @throws std::runtime_error if k is negative or greater than n.
 */
double distribution_binomial::pdf(int &k) {

    // Check parameter is not negative or greater than n
    if (k < 0){
        throw runtime_error("Binomial cannot have negative number of successes.");
    } else if (k > n){
        throw runtime_error("Binomial cannot have more succeses than trials.");
    }

    if (n>10){
        distribution_normal binomial_approx(n*p, sqrt(n*p*(1-p)));
        double pk = k;
        return binomial_approx.pdf(pk);
    }

    // binomial coefficient
    // Only the tails of the coefficient survives cancellation
    int val = functions_math::binomial_coeficient(n,k);

    // return mass probability
    return val*pow(p, k)*pow(p, n-k);

};


/**
 * Calculates the cumulative distribution function (CDF) of the binomial distribution.
 * 
 * @param k The number of successes.
 * @return The probability of having k or fewer successes in the binomial distribution.
 * @throws std::runtime_error if k is negative or greater than the number of trials (n).
 */
double distribution_binomial::cdf(int &k){

    // Check parameter is not negative or greater than n
    if (k < 0){
        throw runtime_error("Binomial cannot have negative number of successes.");
    } else if (k > n){
        throw runtime_error("Binomial cannot have more succeses than trials.");
    }

    // Use the normal approximation if n is large if have tested that for n<100 the incomplete beta function is still accurate
    if (n>max(250*(1-p)/p, 250*p/(1-p)) and n > 100){
        distribution_normal binomial_approx(n*p, sqrt(n*p*(1-p)));
        double pk = k+.5;
        return binomial_approx.cdf(pk);

    } else {
    
        // Use the Beta Distribution to calculate the CDF
        double pa = n-k;
        double pb = k+1.;
        double px = 1.-p;
        double eps = 1e-18;
        double res = functions_math::func_beta_incomplete(pa,pb, px,eps);

        return res;
    }

};