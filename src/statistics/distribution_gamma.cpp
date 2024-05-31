#include "distribution_gamma.h"
#include <stdexcept>
#include <vector>
#include <math.h>
#include "../functions/functions_math.h"
#include <functional>
#include "../statistics/cdflib.hpp"


using namespace std;
using namespace functions_math;

// Constructor
distribution_gamma::distribution_gamma(double palpha, double pbeta){
    alpha = palpha;
    beta = pbeta;
}

// Constructor Overload
// Log likelihood estimation of binned data
distribution_gamma::distribution_gamma(vector<double> bins, vector<int> freqs){

    // Test if bins and freqs are of the same length
    if (bins.size() == freqs.size()){
        bins.insert(bins.begin(), 0.0); // Add 0.0 to the front of freqs
    } else if (bins.size() != freqs.size()+1){
        throw runtime_error("The number of bins and frequencies must be equal or differ by 1.");   
    }

    // Test if bins are in ascending order
    for (int i=0; i<bins.size()-1; i++){
        if (bins[i] >= bins[i+1]){
            throw runtime_error("Bins must be in ascending order.");
        }
    }

    // Solve for alpha using bisection
    double denominator = 0;
    for (int i=0; i<freqs.size(); i++){
        denominator += freqs[i];
    }
    double numerator = 0;
    for (int i=0; i<freqs.size(); i++){
        numerator += freqs[i]*log(bins[i+1]-bins[i]);
    }
    double func_target = numerator/denominator;

    // alpha can only be positive and now we have a function with support on the whole real line
    function<double(double)> func = [=](double alpha)  {return func_target - log(beta) + func_digamma(alpha) - log(alpha-1);};


    // Use the lambda function to solve for alpha using bisection or other numerical methods
    double palpha;
    double eps = 1e-10;

    func_solve_root_bisection(func, 0.3, 0.7, palpha,  eps);
    // set alpha
    alpha = log(palpha);

    func_target = 0;
    for (int i=0; i<bins.size()-1; i++){
        func_target += freqs[i]*(alpha-log(bins[i+1] - bins[i]));   
    }
    func = [=](double pbeta)  {
        double sum = 0;
        for (int i=0; i<bins.size()-1; i++){
            double p1 = alpha +1.;
            double p2 = bins[i]/pbeta;
            double p3 = bins[i+1]/pbeta;
            double p4 = alpha;
            sum += freqs[i]*pbeta*alpha*(func_gamma_incomplete(p1, p3) - func_gamma_incomplete(p1, p2)/(func_gamma_incomplete(p4, p3) - func_gamma_incomplete(p4, p2)));  
        }
        return sum - func_target;
    };
    double pbeta;
    // Find lower beta and upper beta where func has opposite signs
    func_solve_root_bisection(func, 1e-6, 1e+18, pbeta,  eps);
    beta = pbeta;
}

// Destructor
distribution_gamma::~distribution_gamma(){
}

// PDF
double distribution_gamma::pdf(double &x){
    return pow(beta, alpha)/tgamma(alpha)*pow(x, alpha-1)*exp(-beta*x);
}

// CDF
double distribution_gamma::cdf(double &x){
    double px = beta*x;
    double palpha = alpha;
    if (px <= 0){
        return 0;
    }
    if (palpha <= 0){
        return 0;
    }
    if (px == 1){
        return 1;
    }
    double res1 =func_gamma_incomplete(px, palpha);
    double res2 = tgammal(palpha);
    return res1/res2;
}
