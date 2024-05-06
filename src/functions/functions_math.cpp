//
// Created by Jonas Wolff on 28/02/2024.
//

#include <numeric>
#include "algorithm"
#include "math.h"
#include "../algorithms//toms179.h"

#include "functions_math.h"
#include "../statistics/distribution_normal.h"
#include "../statistics/distribution_gamma.h"
#include "../statistics/cdflib.hpp"

using namespace std;

double functions_math::func_digamma(double x){
    // Digamma function see J. M. Bernardo, AS 80, 580 (1985)

    double digamma = 0.;
    int ifault = 1;
    double r;
    // Constants
    // s is a small number, c is a threshold for the stirling formula
    // s3, s4, s5 are stirlngs nth coifficients, d1 is digamma(1.)
    double s = 1e-5, c = 8.5, s3 = 8.333333333e-2, s4 = 8.333333333e-3, s5 = 8.968253968e-3, d1 = -0.5772156649;
    if (x<=0){
        throw runtime_error("Digamma function is not defined for x <= 0.");
    }

    // Calulate digamma function only if x <= a
    if (x >= s){ goto S10;}
    digamma =  d1 - 1.0/x;
    return digamma;
    // Increse x for assymptotic accuracy
S10:
    if (x >= c){ goto S20;}
    digamma = digamma - 1.0/x;
    x = x + 1.0;
    goto S10;

    // Use stirling formula for large x
S20:
    r = 1.0/x;
    digamma = digamma + log(x) - 0.5*r;
    r *= r;
    digamma = digamma - r*(s3 - r*(s4 - r*s5));
    return digamma;
}

double func_binomial_coeficient_logmethod(int n , int k){

    // Set start of product
    double val = log(1);

    // Only the tails of the factorials survive
    for (double i = 1; i <= k or i <= n-k; i++){
        val -= log(i);
    }
    for (double i = n; i >= k or i >= n-k; i--){
        val += log(i);
    }

    // transform to natural logarithm
    val = exp(val);

    // Return val
    return val;
}

void func_NM_2D_assign_order(vector<double> &func_simplex, vector<int> &simplex_order){
    // Order the simplex
    for (int i_=2;i_>0; i_--){
        if (func_simplex[simplex_order[1]] > func_simplex[simplex_order[2]]){
            swap(simplex_order[1], simplex_order[2]);
        }
        if (func_simplex[simplex_order[0]] > func_simplex[simplex_order[1]]){
            swap(simplex_order[0], simplex_order[1]);
        }
    }
}

// Minimize function with Nelder-Mead in 2D
void functions_math::func_minimize_NM_2D(std::function<double(double, double)> func, double x0, double y0, double &x1, double &y1, double &eps){

    //Define parameters
    double alpha = 1.;
    double beta = 0.5;
    double gamma = 1.;

    // Define the simplex
    vector<vector<double>> simplex = {{x0, y0}, {x0+1., y0}, {x0, y0+1.}};
    vector<int> simplex_order = {2, 1, 0};

    // Calculate the function values
    vector<double> func_simplex;
    for (int j=0; j<3; j++){
        func_simplex.push_back(func(simplex[j][0], simplex[j][1]));
    }

    func_NM_2D_assign_order(func_simplex, simplex_order);

 
    for (int i=1; i<100000; i++){
            
        // Calculate the centroid
        vector<double> centroid = {0., 0.};
        for (int j=0; j<3; j++){
            if (j != simplex_order[2]){
                centroid[0] += simplex[j][0];
                centroid[1] += simplex[j][1];
            }
        }
        centroid[0] /= 2.;
        centroid[1] /= 2.;

        // Reflect
        vector<double> reflect = {0., 0.};
        reflect[0] = (1.+alpha)*centroid[0] - alpha*simplex[simplex_order[2]][0];
        reflect[1] = (1.+alpha)*centroid[1] - alpha*simplex[simplex_order[2]][1];
        double f_reflect = func(reflect[0], reflect[1]);

        // Check if reflect is the best
        if (f_reflect < func_simplex[simplex_order[0]]){
            // Expand
            vector<double> expand = {0., 0.};
            expand[0] = (1.+gamma)*reflect[0] - gamma*centroid[0];
            expand[1] = (1.+gamma)*reflect[1] - gamma*centroid[1];
            double f_expand = func(expand[0], expand[1]);

            if (f_expand < f_reflect){
                simplex[simplex_order[2]] = expand;
                func_simplex[simplex_order[2]] = f_expand;
                func_NM_2D_assign_order(func_simplex, simplex_order);
            } else {
                simplex[simplex_order[2]] = reflect;
                func_simplex[simplex_order[2]] = f_reflect;
                func_NM_2D_assign_order(func_simplex, simplex_order);
            }
        } else if (f_reflect < func_simplex[simplex_order[1]]){
            simplex[simplex_order[2]] = reflect;
            func_simplex[simplex_order[2]] = f_reflect;
            func_NM_2D_assign_order(func_simplex, simplex_order);
        } else {
            // Contract
            vector<double> contract = {0., 0.};
            contract[0] = beta*simplex[simplex_order[2]][0] + (beta-1.)*centroid[0];
            contract[1] = beta*simplex[simplex_order[2]][1] + (beta-1.)*centroid[1];
            double f_contract = func(contract[0], contract[1]);

            if (f_reflect <= func_simplex[simplex_order[2]]){
                simplex[simplex_order[2]] = reflect;
                func_simplex[simplex_order[2]] = f_reflect;
                func_NM_2D_assign_order(func_simplex, simplex_order);
            } 
            if (f_contract < func_simplex[simplex_order[2]]){
                simplex[simplex_order[2]] = contract;
                func_simplex[simplex_order[2]] = f_contract;
                func_NM_2D_assign_order(func_simplex, simplex_order);
            } else {
                // Shrinks
                for (int j=0; j<3; j++){
                    if (j != simplex_order[0]){
                        simplex[j][0] = simplex[simplex_order[0]][0] + 0.5*(simplex[j][0] - simplex[simplex_order[0]][0]);
                        simplex[j][1] = simplex[simplex_order[0]][1] + 0.5*(simplex[j][1] - simplex[simplex_order[0]][1]);
                    }
                }
                func_simplex[simplex_order[1]] = func(simplex[simplex_order[1]][0], simplex[simplex_order[1]][1]);
                func_simplex[simplex_order[2]] = func(simplex[simplex_order[2]][0], simplex[simplex_order[2]][1]);
                func_NM_2D_assign_order(func_simplex, simplex_order);

            }
        }

        // test for convergence
        if (abs(func_simplex[simplex_order[2]] - func_simplex[simplex_order[0]])/(abs(func_simplex[simplex_order[2]])+abs(func_simplex[simplex_order[0]])) < eps){
            break;
        }

    }

    // test for convergence
    if (abs(func_simplex[simplex_order[2]] - func_simplex[simplex_order[0]])/(abs(func_simplex[simplex_order[2]])+abs(func_simplex[simplex_order[0]]))  > eps){
        throw runtime_error("Nelder-Mead did not converge within 100000 iterations.");
    }
    
    // Set the minimum
    x1 = simplex[0][0];
    y1 = simplex[0][1];
}

void functions_math::func_fit_histogram_gamma(vector<double> bins, vector<double> freqs, double &alpha, double &beta){
    // Fit a histogram to a gamma distribution
    double eps = 1e-16;
    
    // Check that bins and freqs are the same size
    if (bins.size() != freqs.size()){
        throw runtime_error("Bins and freqs are not the same size.");
    }

    // Check that bins are sorted
    for (int i=0; i<bins.size()-1; i++){
        if (bins[i] > bins[i+1]){
            throw runtime_error("Bins are not sorted.");
        }
    }

    // Check that bins are positive
    for (int i=0; i<bins.size(); i++){
        if (bins[i] < 0){
            throw runtime_error("Bins are not positive.");
        }
    }

    // Check that freqs are positive
    for (int i=0; i<freqs.size(); i++){
        if (freqs[i] < 0){
            throw runtime_error("Freqs are not positive.");
        }
    }


    // Define the function to fit
    auto func = [&](double x, double a, double b){
        return b*pow(x,a-1.)*exp(-x*b);
    };

    // Define the function to minimize
    auto func_min = [&](double a, double b){

        if (b <= 0){
            return 1e-18;
        } 
        
        distribution_gamma gamma_dist(abs(a), abs(b));  

        double sum = 1.;
        for (int i=0; i<bins.size()-1; i++){
            if (freqs[i] == 0){
                continue;
            }
            sum *= pow((gamma_dist.cdf(bins[i+1]) - gamma_dist.cdf(bins[i]))/freqs[i],freqs[i]);
        }
        if (freqs[freqs.size()-1] != 0){
            sum *= pow((1.-gamma_dist.cdf(bins[bins.size()-1]))/freqs[freqs.size()-1],freqs[freqs.size()-1]);
        }
        return -sum;
    };

    // Define the starting point
    double a0 = 5.;
    double b0 = 7.;

    // Minimize the function
    functions_math::func_minimize_NM_2D(func_min, a0, b0, alpha, beta, eps);
}

// Finds root of function via secant method
void functions_math::func_solve_root_secant(std::function<double(double)> func, double x1, double x0, double &x2, double &eps, double iterationscale){
    
        double f0 = func(x0);
        double f1 = func(x1);
    
        for (int i = 1; i < 10000; i++){
    
            x2 = x1 - iterationscale*f1*(x1-x0)/(f1-f0);
            f0 = f1;
            f1 = func(x2);
            x0 = x1;
            x1 = x2;
    
            if (abs(f1) < eps){
                return;
            }
    
        }
    
        throw runtime_error("Secant method did not convert within 10000 iterations.");
}

// Finds root of function via bisection
void functions_math::func_solve_root_bisection(std::function<double(double)> func,double xL, double xU, double &xN, double &eps){

    double fmid;
    double fL = func(xL);
    double fU = func(xU);
    // Order interval end points so that func(xU) is positive
    if (fL >  fU){
        swap(xL, xU);
        swap(fL, fU);
    }
    // check if function is guaranteed a root
    if (fU < 0 or fL > 0){
        throw runtime_error("Bisection needs interval end points to have both negative and positive function values.");
    }

    
    for (int i = 1; i < 10000; i++){

        xN = (xU+xL)/2.;
        fmid = func(xN);

        if (fmid >= 0){
            xU = xN;
        } else {
            xL = xN;
        }

        if (abs(fmid) < eps) {
            return;
        }

    }

    throw runtime_error("Bisection did not convert within 10000 iterations.");
}

double func_riemann_zeta(double s){

    double sum = 0.;
    int i;
    for (i = 1; i < 10000; i++){
        sum += 1./pow(i,s);

        if (abs(1./pow(i,s)) < 1e-18*sum){
            break;
        
        }
    }

    if (i >= 10000){
        throw runtime_error("Riemann Zeta function did not converge within 10000 iterations.");
    }

    return sum;
}

// this is B_n^+
double func_bernoulli_numbers(int n){
    
    if (n == 0){
        return 0.;
    }

    if (n == 1){
        return .5;
    }

    if (n % 2 == 1){
        return 0.;
    }
    
    // See David Harvey's paper for the algorithm
    double res = pow(-1.,n/2+1.)*2.*func_riemann_zeta(n);
    for (int i = 1; i <= n-1; i++){
        res *= i*M_1_PI/2.;
    }
   
    return res;

}

double functions_math::binomial_coeficient(int n, int k){
    // binomial coefficient

    // Set start of product
    double val = 1;

    // For large numbers do log transform
    if ( k > 10){
        // Use Normal Distribution instead
        return func_binomial_coeficient_logmethod(n, k);
    }

    // Only the tails of the factorials survive
    for (double i = 1; i <= k and i <= n-k; i++){
        val /= i;
    }
    for (double i = n; i > k and i > n-k; i--){
        val *= i;
    }

    return val;

};





// Calculates the minimax approximation for lgamma(a) - (a - .5)*ln(a) + a - .5*ln(2*M_PI)
double func_beta_DELTA(double &a){
    

    if (a >= 20.){
        return  1./(12.*a) - 1./(360.*pow(a,3.)) + 1./(1260.*pow(a,5.)) - 1./(1680.*pow(a,7.));
        return lgamma(a) - (a - 0.5)*log(a) + a - 0.5*log(2.*M_PI);
    } else if (a >= 15.) {
        double c0 =  .833333333333333e-1;
        double c1 = -.277777777770481e-2;
        double c2 =  .793650663183693e-3;
        double c3 = -.595156336428591e-3;
        double c4 =  .820756370353826e-3;
        return c0/pow(a,1.)+c1/pow(a,3.)+c2/pow(a,5.)+c3/pow(a,7.)+c4/pow(a,9.);

    } else if (a >= 8.) {
        double c0 =  .833333333333333e-1;
        double c1 = -.277777777760991e-2;
        double c2 =  .793650666825390e-3;
        double c3 = -.595202931351870e-3;
        double c4 =  .837308034031215e-3;
        double c5 = -.125322962780713e-2;
        return c0/pow(a,1.)+c1/pow(a,3.)+c2/pow(a,5.)+c3/pow(a,7.)+c4/pow(a,9.)+
                c5/pow(a,11.);

    } else {
        throw runtime_error("minimax approximation was invoked with parameter < 8. This is not optimal behavior.");
    }
}

// Calculates ALGDIV for ln(tgamma(b)/tgamma(a+b))
double func_beta_ALGDIV(double &a, double &b){

    // Check for optimality
    if (a < 0 or b <8){
        throw runtime_error("ALGDIV(a,b) is not optimal for a < 0 or b <8.");
    }

    double c[6] =  {.833333333333333e-1,
                    -.277777777760991e-2,
                    .793650666825390e-3,
                    -.595202931351870e-3,
                    .837308034031215e-3,
                    .125322962780713e-2};
    double w = 0.;
    for (int i = 0.; i <= 5; i++){
        w += c[i] *(1./pow(b,i*2.+1.) - 1/pow(a+b,i*2.+1.));
    }
    return w - (a + b - 0.5)*log(1. + a/b) - a*(log(b) - 1.);
}

// Calculate auxiliary function for BETALN and  BRCOMP
double func_beta_BCORR(double &a, double &b){

    // Check for optimality
    if (a < 8 or b <8){
        throw runtime_error("BCORR(a,b) is not optimal for a < 8 or b <8.");
    }

    double a0 = min(a,b);
    double b0 = max(a,b);
    double c[6] =  {.833333333333333e-1,
                    -.277777777760991e-2,
                    .793650666825390e-3,
                    -.595202931351870e-3,
                    .837308034031215e-3,
                    .125322962780713e-2};
    double w = 0.;
    for (int i = 0; i <= 5; i++){
    w += c[i]*(1./pow(b0,i*2.+1.) - 1/pow(a0+b0,i*2.+1.));
    }
    double checkval2  = func_beta_DELTA(a0);
    return  checkval2 + w;
}

// Calculate -log(func_beta_G(a,b))
double func_beta_BETALN(double a, double b){

    if (a > b){
        swap(a,b);
    }

    // Check for optimality
    if (a <= 0 or b <= 0){
        throw runtime_error("BETALN(a,b) cannot be calculated efficiently for a <= 0 or b <= 0");
    }

    if (a >= 8){
        double v = b*log(1. + a/b);
        double u = -(a - 0.5)*log(a/(a + b));
        return (.5*log(2.*M_PI) - 0.5*log(b)) + func_beta_BCORR(a,b) - u - v;
    } else if (2 < a && a < 8){
        return (a-1.)/(a+b-1.)*func_beta_BETALN(a-1.,b);
    } else if (b >= 8){
        return lgamma(a) + func_beta_ALGDIV(a,b);
    } else if (2 < b && b < 8){
        return (b-1.)/(a+b-1.)*func_beta_BETALN(a,b-1.);
    } else if (a >= 1) {
        return lgamma(a) + lgamma(b) - lgamma(a+b);
    } else {
        throw runtime_error("BETALN(a, b) is not accurate for a < 1");
    }

}

// Calculates the beta function
double functions_math::func_beta(double &a, double &b) {

    /*if (fmod(a,1)<0.000001 && fmod(b,1)<0.000001) {
        double B = 1;
        for (int i = 1; i <= a+b; i++){
            double di = i;
            if (i <= min(a,b)){
                B /= di;
            } else if (i > min(a,b) and i <= max(a,b)){
                continue;
            } else {
                B *= di;
            }
        }
        return B;
    } else */
    double res = (tgammal(a)*tgammal(b))/tgammal(a+b);

    if (std::isnan(res)) {
        if (max(a,b)>8 && a > 1){
            res = exp(func_beta_BETALN(a,b));
            //res = exp((lgammal(a)+lgammal(b))-lgammal(a+b)); 
        } else {
            res = exp((lgammal(a)+lgammal(b))-lgammal(a+b)); 
        }
    }
        
    return res;
}

// Auxiliary function
double func_beta_phi(double &x){
    if (x <= 0){
        throw runtime_error("func_beta_phi(x) only takes argument x > 0.");
    }
    return x - 1. - log(x);
}

// Calculates the inverse of the beta function
long double func_beta_G(double &a, double &b) {
    
    if (min(a,b) < 8){ 
        long double G = 1;
        double valab, vala, valb;
        if (a+b >= 1){
            valab = tgammal(1.+fmod(a+b,1.));
        } else {
            valab = tgammal(fmod(a+b,1.));
        }
        if (a >= 1){
            vala = tgammal(1.+fmod(a,1.));
        } else {
            vala = tgammal(fmod(a,1.));
        }
        if (b >= 1){
            valb = tgammal(1.+fmod(b,1.));
        } else {
            valb = tgammal(fmod(b,1.));
        }

        for (int i = 1; i <= a+b-1.; i++){
            double di = i;
            if (i <= a-1.){
                G /= di+fmod(a,1.);
            } 
            if (i <= b-1.){
                G /= di+fmod(b,1.);
            } 
            if (i <= a+b-1.) {
                G *= di+fmod(a+b,1.);
            }
        }
        

        double frac = valab/vala/valb;
        return G*frac;
    } else { 

        double res = tgammal(a+b)/(tgammal(a)*tgammal(b));

        if (std::isnan(res)) {
            if (max(a,b)>8 && a > 1){
                res = exp(-func_beta_BETALN(a,b));
            } else {
                res = exp(lgammal(a+b) -lgammal(a) - lgammal(b)); 
            }
        }
            
        return res;
        

    }
}

double func_beta_BRCOMP(double &a, double &b, double &x, double &y){

    // Optimality check
    if (a >= 0 and b >= 0){}
    else if (0 < x and x < 1){}
    else if (y == (1.-x)){}
    else {
        throw runtime_error("BRCOMP(a,b) only accurate for a,b >= 0 or 0 < x <1 or y = 1-x");
    }

    if (min(a,b) < 8){
        return pow(x,a)*pow(y,b)*func_beta_G(a,b);
    }

    //double lambda = a - (a + b)*x;
    double p = a/(a+b);
    double lambda = (a+b)*(p-x);
    double p1 = 1. - lambda/a;
    double p2 = 1. + lambda/b;
    double z = (a* func_beta_phi(p1) + b* func_beta_phi(p2)) + func_beta_BCORR(a,b);
    double checkval3 = log(1.+lambda/b);
    double checkval4 = log(1.-lambda/a);
    double my_bcorr = func_beta_BCORR(a,b);
    double true_bcorr = bcorr ( &a, &b );
    //double z = ( - a*log(1.-lambda/a)  -b*log(1.+lambda/b)) + func_beta_BCORR(a,b);
    
    return sqrt(a/(a+b)*b)*M_SQRT1_2*M_2_SQRTPI*exp(-z) ;

}


// Calculate the incomplete beta function by using taylor series on the (1-t)^(b-1) term
double func_beta_BPSER(double &a, double &b, double &x, double &eps){

    // Calculate taylorsum
    double taylorsum = 0.;
    double numerator = 1.;
    double frac = 1./1.;
    double summand;
    int maxexpansion = 10000;
    for (double j=1.; j<maxexpansion; j++){

        numerator *= (j-b);
        frac *= (j-b)/j;
        summand = frac/(a+j)*pow(x,j);
        //summand = numerator*exp(-lgamma(j+1.)-log(a+j))*pow(x,j);
        taylorsum += summand;

        // Check if sum is adding anything to precision
        if (abs(summand) <= eps*abs(taylorsum) && j > 10){
            break;
        }
    }

    // Check that we got wanted precision
    if (abs(summand) > eps*abs(taylorsum)){
        throw runtime_error("Beta function calculation with taylor approximation did not converge fast enough.");
    }

    // Return beta function with taylor approximation
    return func_beta_G(a,b)*pow(x,a)/a*(1. + taylorsum*a);

};

// Calculate the series expansion of the incomplete beta function by integrating the remaining tail
double func_beta_BPSER_reminder(double &a, double &b, double &x, double &eps){
    double py = 1. - x;
    return 1. - func_beta_BPSER(b,a,py, eps);
};

// Calculates the maclaurin expansion of the normalized incomplete  gamma function
void func_beta_BGRAT(double &a, double &b, double &x, double &y, double &w, double &eps, bool &IERR){
    
    if (a<=b){
        throw runtime_error("BGRAT(a,b) only works for a > b.");
    }
    
    IERR = false;
    double w_0 = w;
    double T = a + (b-1.)/2.;
    double u = -T*log(x);
    function<double(double, double)> H = [=](double c, double u)  {return exp(-u)*pow(u,c)/ tgamma(c);};
    //double M = H(b, u)* tgamma(a+b)/ (tgamma(a)*pow(T,b));
    double M = exp(-u)*pow(u,b)*func_beta_G(a,b) *pow(T,-b);

    // start values
    double _ = H(b,u);
    double J_i = functions_math::func_gamma_incomplete(b, u)/H(b,u);
    vector<double> p;
    double p_i;
    p.push_back(1.);
    double I_x = M*p[0]*J_i;

    // Increase precision
    for (int i=1; i<10000; i++){

        // Get next p_i
        p_i = (b-1.)/tgamma(2.*i+2.);
        for (int m=1; m<=i-1; m++){
            p_i += 1./i*(m*b-i)/tgamma(2.*m+2.)*p[i-m];
        }
        p.push_back(p_i);

        // Get next J_i
        double j = i-1;
        J_i = (b+2.*j)*(b+2.*j+1.)/(4.*pow(T,2.))*J_i;
        J_i += (u+b+2.*j+1.)/(4.*pow(T,2.))*pow(log(x)/2.,2.*j);

        // Check if machine epsilon has been hit
        if (M*p_i*J_i <= eps*I_x){
            I_x += M*p_i*J_i;
            IERR = true;
            break;
        }

        // Add to result
        I_x += M*p_i*J_i;

    }

    // Check if we ended prematurely
    if (IERR == false){
        // Tell the algorithm that we did not encounter underflow
        throw runtime_error("The incomplete beta function could not be calculated do to BGRAT not encountering underflow.");
    }

    // return maclaurin polynomial
    w = I_x +w_0;
}

// Calculate effect of changes to a in the incomplete beta function
long double func_beta_BUP(double &a, double &b, double &x, double &y, double &n, double &eps){

    // Check n >= 1
    if (n < 1){
        throw runtime_error("func_beta_BUP cannot take n < 1.");
    }

    // Calculate sum
    long double d_i = 1.;
    long double h_i = 1.0;
    long double r_i;
    //long double lognormalize = log(func_beta_BRCOMP(a,b,x,y)/a);
    long double sumval = h_i;
    bool h_isdecreasing = b<=1.;
    for (long double i = 1; i<=n-1.; i++){

        d_i *= (a+b+i-1.)/(a+1.-1.+i);
        h_i = d_i*pow(x,i);
        r_i = (a+1.-1.+i)/(a+b+i-1.);
        sumval += h_i;


        // Stump sum if we have hit machine epsilon and all the rest summands are under machine epsilon
        if (h_isdecreasing && h_i <= eps*sumval){
            break;
        }

        // Check if h has started to decrease
        if (h_isdecreasing == false && x >= r_i){
            h_isdecreasing = true;
        }
    }

    return sumval*func_beta_BRCOMP(a,b,x,y)/a;
};


// Calculate the continued fraction approximation, with numerical stabilization trick
double func_beta_BFRAC(double &a, double &b, double &x, double &y, double &lambda, double &eps){

    double p_ = a/(a+b);
    // Check if function approximates well
    if (x > p_){
        throw runtime_error("BFRAC does not approximate accurately for x > a/(a+b)");
    }
    double frac_ = func_beta_BRCOMP(a,b,x,y);
    double frac_true = beta_rcomp ( &a, &b, &x, &y );
    if (frac_ == 0){
        return 0.;
    }
    double lambda_ = (a+b)*(p_-x);
    double ai = 1.;
    double bi = a/(a+1.)*(lambda_+1.);
    double Aj = 0.; // A_{i-1}
    double Ai = ai;
    double Bj = 1.; // B_{i-1}
    double Bi = bi;
    double A,B; // A_{i+1}, B_{i+1}

    for (double i = 1.; i < 10000; i++){
        ai = (a+i-1.)*(a+b+i-1.)/pow(a+2.*i-1.,2.)*i*(b-i)*pow(x,2.);
        bi = i + i*(b-i)*x/(a+2.*i-1.) + (a+i)/(a+2.*i +1.)*(lambda_ + 1. + i*(1.+y));
        A = bi*Ai + ai*Aj;
        B = bi*Bi + ai*Bj;
        Aj = Ai;
        Bj = Bi;
        Ai = A;
        Bi = B;

        // Check if we have hit wanted precision.
        if (abs(Ai/Bi-Aj/Bj) < eps*abs(Ai/Bi)){
            break;
        }
    }
    double checkval1 = Ai/Bi;
    double res = frac_*Ai/Bi;
    double res_true =beta_frac ( &a, &b, &x, &y , &lambda, &eps);
    return res;

}




double func_beta_calculateLnU(double a, double b) {
    double ab = a + b;
    double delta_ab = func_beta_DELTA(ab);
    double delta_a = func_beta_DELTA(a);
    double delta_b = func_beta_DELTA(b);
    double lnU = delta_ab - delta_a - delta_b;
    return lnU;
}


double func_beta_BASYM_b(double n, double r, vector<double> an){
    if (n==0){
        return 1.;
    } else if (n==1){
        return r*an[1];
    } else {
        double b = r*an[1];
        for (int i = 1; i <= n-1; i++){
            b += r*an[n] + (1.0/n)*((n-i)*r-i)*an[n-i]*b;
        }
        return b;
    }
    
    double b = 1.;
    for (int i = 1; i <= n-1; i++){
        b *= (n-i+1.)/i;
    }
    return b;

}

// Assymptotic expansion of the incomplete beta function
double func_beta_BASYM(double &a, double &b, double &x, double &y, double &lambda, double &eps){
    double p = a/(a+b);
    double q = 1 - p;

    if (x > p){
        throw runtime_error("BASYM does not approximate accurately for x > a/(a+b)");
    }

    double brcomp = func_beta_BRCOMP(a, b, p, q);
    
    //double U = brcomp * M_SQRT2 * sqrt( M_PI * ( (a + b) / a / b));
    double U = exp(func_beta_calculateLnU(a, b));
    //double U = pow(p,a)*pow(q,b)*func_beta_G(a,b)*sqrt(2*M_PI*(a+b)/a/b);

    double beta_gamma;
    if (a <= b) {
        beta_gamma = std::sqrt(q / a);
    } else {
        beta_gamma = std::sqrt(p / b);
    }

    auto func_phi = [&](double t) {
        return -(a * log(t / p) + b * log((1 - t) / q));
    };
    double z = sqrt(func_phi(x));


    vector<double> an;
    vector<double> bn;
    vector<double> en;
    vector<double> cn;
    vector<double> Ln;
    
    double cursum = 0;
    for (int n = 0; n <= 10000; n++) {
        double an_d;
        if (a <= b) {
            an_d = (2.0 / (2. + n)) * q * (1. + pow(-1., n) * pow(a / b, n + 1.));
        } else {
            an_d = (2.0 / (2. + n)) * p * (pow(-1., n) + pow(b / a, n + 1.));
        }
        an.push_back(an_d);
        double r = -n/2.;
        if (n==0){
            bn.push_back(1.);
            en.push_back(1.);
            // cn is not used in the first iteration
            cn.push_back(-1.);
            cn.push_back(1.);
            Ln.push_back((sqrt(M_PI) / 4.) * exp(pow(z, 2.)) * erfc(z));
        }
        else if (n==1){
            bn.push_back(r * an[1]);
            Ln.push_back(M_SQRT1_2/2.); 
            //Ln.push_back(pow(2.,-3./2.)); 
        }
        else { 
            bn.push_back(func_beta_BASYM_b(n, -n/2., an));
        }
        
        if (n>=1){
            cn.push_back( (1.0 / (1.+n))  * bn[n]);
            int i = -1;
            en.push_back( -std::accumulate(en.begin(), en.end(), 0.0, [&](double sum, double ei) {
                i++;
                return sum + ei * cn[n - i + 1];
            }));

            Ln.push_back(M_SQRT1_2/2. * pow(M_SQRT2 * z, n - 1.) + (n - 1.) * Ln[n-2]  );
        }

        if (abs(en[n]*Ln[n]*pow(beta_gamma,n)) <= abs(eps*cursum)){
            cursum += en[n] * Ln[n] * pow(beta_gamma, n);
            break;
        } else {
            cursum += en[n] * Ln[n] * pow(beta_gamma, n);
        }
        
    }

    
    double res = M_2_SQRTPI *U* exp(-pow(z, 2.)) * cursum;

    return res;
}



double functions_math::func_beta_incomplete(double a, double b, double x, double &eps){

    double y = 1.-x;

    // Test for the real values of the incomplete beta function
    double w = 0, wl = 0;
    int ierr = 0;
    beta_inc(&a,&b,&x, &y, &w, &wl, &ierr);
    return w;
    double p = a/(a+b);
    double q = 1. - p;


    // run methods for min(a,b) <= 1
    if (min(a,b) <= 1) {

        bool invuse = false;
        double return_res;

        /*
        if (x>0.5){
            swap(a,b);
            swap(x, y);
            p = a/(a+b);
            invuse = true;
        }*/

        // Check if taylor approximation is enough
        double res;
        if (max(a, b) <= 1 && a >= min(0.2, b)) {
            res = func_beta_BPSER(a, b, x, eps);

        } else if (max(a, b) <= 1 && a < min(0.2, b) && pow(x, a) <= 0.9 ) {
            res = func_beta_BPSER(a, b, x, eps);

        } else if (max(a, b) > 1 && b <= 1 ) {
            res = func_beta_BPSER(a, b, x, eps);

        } else if (max(a, b) > 1 && b > 1 && x < 0.1 && pow(b * x, a) <= 0.7 ) {
            res = func_beta_BPSER(a, b, x, eps);

        } else if (max(a, b) <= 1 && a < min(0.2, b) && pow(x, a) > 0.9 && x >= 0.3) {
            res = func_beta_BPSER_reminder(a, b, x, eps);

        } else if (max(a, b) > 1 && b > 1 && x >= 0.3) {
            res = func_beta_BPSER_reminder(a, b, x, eps);
        } else if (max(a,b) > 1 && b>15 && 0.1 <= x && x < 0.3){
            double w = 0.;
            double eps15 = 15.*eps;
            bool IERR = false;
            func_beta_BGRAT(b,a,y,x,w,eps15,IERR);
            res = 1. - w;
        } else if (max(a,b) > 1 && b>15 && 0.1 < x && pow(b*x,a) > 0.7){
            double w = 0.;
            double eps15 = 15.*eps;
            bool IERR = false;
            func_beta_BGRAT(b,a,y,x,w,eps15,IERR);
            res = 1. - w;
        } else if (max(a,b)>1 && b>1 && 0.1<=x && x < 0.3 && b <= 15){
            double n = 20.;
            double bpn = b + n;
            double w = func_beta_BUP(b, a, y, x, n, eps);
            double eps15 = 15.*eps;
            bool IERR = false;
            func_beta_BGRAT(bpn,a,y,x,w,eps15,IERR);
            res = 1. - w;
        } else if (max(a,b)> 1 && b > 1 && x < 0.1 && pow(b*x, a)>0.7 && b <= 15){
            double n = 20.;
            double bpn = b + n;
            double w = func_beta_BUP(b, a, y, x, n, eps);
            double eps15 = 15.*eps;
            bool IERR = false;
            func_beta_BGRAT(bpn,a,y,x,w,eps15,IERR);
            res = 1. - w;
        } else if (max(a,b) <= 1 && a < min(0.2, b) && pow(x,a) > 0.9 && x < 0.3){
            double n = 20.;
            double bpn = b + n;
            double w = func_beta_BUP(b, a, y, x, n, eps);
            double eps15 = 15.*eps;
            bool IERR = false;
            func_beta_BGRAT(bpn,a,y,x,w,eps15,IERR);
            res = 1. - w;
        } else {
            throw runtime_error("Incomplete beta function has not appropriate approximation");
        }
        
        if (invuse){
            return 1. - res;
        } else {
            return res;
        }
    } else if (min(a,b) > 1){

        bool invuse = false;
        double return_res;

        if (x>p){
            swap(a,b);
            swap(x, y);
            p = a/(a+b);
            invuse = true;
        }

        double n = floor(b);
        double b_ = b - n;
        if (n==b){
            n = b -1.;
            b_= 1.;
        }

        double w;
        bool IERR;
        double eps15 = 15.*eps;
        double eps100 = 100.*eps;
        double lambda = (a+b)*(p-x);

        if ( b < 40 && b*x < 0.7){
            return_res = func_beta_BPSER(a, b, x, eps);
        } else if (b <40 && b*x >0.7 && x <= 0.7) {
            return_res = func_beta_BUP(b_,a,y,x,n,eps) + func_beta_BPSER(a, b_,x, eps);
        } else if (b < 40 && x > 0.7 && a > 15){
                IERR = false;
                w = func_beta_BUP(b_,a,y,x,n,eps);
                func_beta_BGRAT(a, b_, x, y, w, eps15, IERR);
                return_res = w;
        } else if (b < 40&& x > 0.7 && a <= 15){
            double m = 20.;
            double apm = a+m;
            w = func_beta_BUP(b_,a,y,x,n,eps) + func_beta_BUP(a,b_,x,y,m,eps);
            func_beta_BGRAT(apm,b_,x,y,w,eps15,IERR);
            return_res = w;
        } else if (b >= 40 && a <= b && a <= 100){
            return_res = func_beta_BFRAC(a,b,x,y,lambda,eps15);
        } else if (b >= 40 && 100 < a && a <= b && x < .97*p){
            return_res = func_beta_BFRAC(a,b,x,y,lambda,eps15);
        } else if (b >= 40 && a > b && b <= 100){
            return_res = func_beta_BFRAC(a,b,x,y,lambda,eps15);
        } else if (b >= 40 && 100 < b && b <= a && y > 1.03*q){
            return_res = func_beta_BFRAC(a,b,x,y,lambda,eps15);
        } else if (b >= 40 && 100 < a && a <= b && x >= .97*p){
            return_res = func_beta_BASYM(a,b,x,y,lambda,eps100);
        } else if (b >= 40 && 100 < b && b <= a && y <= 1.03*q){
            return_res = func_beta_BASYM(a,b,x,y,lambda,eps100);
        } else {
            throw runtime_error("Incomplete beta function has not appropriate approximation");
        }

        if (invuse){
            return 1. - return_res;
        } else {
            return return_res;
        }
    }

    // throw error if no method is good enough
    throw runtime_error("No method found for accurately calculating the incomplete beta function.");

};


double func_gamma_incomplete_continued_fraction_G(double &p, double &x){

    double eps = 2.22e-16;
    double a_i,b_i,f,C,D, d_m, DELTA;
    a_i = 1.;
    if (x <= p) {
        b_i = p - 1. + 1.;
    } else {
        b_i =x+2.-1.-p;
    }
    d_m = 10e-300;
    f = a_i/b_i;
    C = a_i/d_m;
    D = 1./b_i;
    int maxiterations = 10000;

    // Calculate continued fraction
    for (int i =2; i<maxiterations;  i++){
        if (i%2==0 and x <= p){
            a_i = -(p-1.+i*1.0/2.);
            b_i = p-1+i;
        } else if (i%2!=0 and x <= p){
            a_i = (i-1.0)/2.*x;
            b_i = p-1.+i;
        } else {
            a_i = -(i-1.)*(i-p-1.);
            b_i = x+2.*i-1.-p;
        }
        D = D*a_i+b_i;
        if (D == 0.){
            D = d_m;
        }
        C = b_i + a_i/C;
        if (C == 0.){
            C  = d_m;
        }
        D = 1./D;
        DELTA = C*D;
        f = f*DELTA;
        if (abs(DELTA -1.)< eps){
            break;
        }
    }

    // Check if continued fraction got wanted precision
    if (abs(DELTA-1.)> eps){
        throw runtime_error("Continued fraction did not converge when calculating the incomplete gamma function.");
    };

    // Return result
    return f;
}

// Calculates the incomplete gamma function (normalized)
/*double functions_math::func_gamma_incomplete(double &p, double &x){
    // https://hal.science/hal-01329669v2/document

    // Check that x is positive
    if (x < 0.){
        throw runtime_error("The incomplete gamma function does not take negative x.");
    }

    // G is lower gamma if x <= p and G is the upper gamma if x > p
    double G = func_gamma_incomplete_continued_fraction_G(p,x)*exp(p*log(x)-x);

    // return incomplete gamma function
    if (x<=p) {
        return tgamma(p) - G;
    } else {
        return G;
    }


};*/

// Calculates the incomplete gamma function (normalized) upper
// Algorithm AS 239 Appl. Statist. (1988) Vol. 37, No. 3 pp. 466-473 (8 pages)
double functions_math::func_gamma_incomplete(double &p, double &x){
  distribution_normal stdN(0,1);
  int *ifault = 0;
  double a;
  double an;
  double arg;
  double b;
  double c;
  double elimit = - 88.0;
  double oflo = 1.0E+37;
  double plimit = 1000.0;
  double pn1;
  double pn2;
  double pn3;
  double pn4;
  double pn5;
  double pn6;
  double rn;
  double tol = 1.0E-14;
  bool upper;
  double value;
  double xbig = 1.0E+08;

  value = 0.0;
//
//  Check the input.
//
  if ( x < 0.0 )
  {
    throw runtime_error("X < 0");
    *ifault = 1;
    return value;
  }

  if ( p <= 0.0 )
  {
    throw runtime_error("P <= 0");
    *ifault = 1;
    return value;
  }

  ifault = 0;

  if ( x == 0.0 )
  {
    value = 0.0;
    return 1.-value;
  }
//
//  If P is large, use a normal approximation.
//
  if ( plimit < p )
  {
    pn1 = 3.0 * sqrt ( p ) * ( pow ( x / p, 1.0 / 3.0 ) 
    + 1.0 / ( 9.0 * p ) - 1.0 );

    upper = false;
    if (upper==true){
        value = 1- stdN.cdf( pn1 );
    } else {
        value = stdN.cdf( pn1 ); 
    }
    return 1. - value;
  }
//
//  If X is large set value = 1.
//
  if ( xbig < x )
  {
    value = 1.0;
    return 1. - value;
  }
//
//  Use Pearson's series expansion.
//
  if ( x <= 1.0 || x < p )
  {
    arg = p * log ( x ) - x - lgamma ( p + 1.0 );
    c = 1.0;
    value = 1.0;
    a = p;

    for ( ; ; )
    {
      a = a + 1.0;
      c = c * x / a;
      value = value + c;

      if ( c <= tol )
      {
        break;
      }
    }

    arg = arg + log ( value );

    if ( elimit <= arg )
    {
      value = exp ( arg );
    }
    else
    {
      value = 0.0;
    }
  }
//
//  Use a continued fraction expansion.
//
  else 
  {
    arg = p * log ( x ) - x - lgamma ( p );
    a = 1.0 - p;
    b = a + x + 1.0;
    c = 0.0;
    pn1 = 1.0;
    pn2 = x;
    pn3 = x + 1.0;
    pn4 = x * b;
    value = pn3 / pn4;

    for ( ; ; )
    {
      a = a + 1.0;
      b = b + 2.0;
      c = c + 1.0;
      an = a * c;
      pn5 = b * pn3 - an * pn1;
      pn6 = b * pn4 - an * pn2;

      if ( pn6 != 0.0 )
      {
        rn = pn5 / pn6;

        if ( (fabs ( value - rn ) <=  tol) && (fabs ( value - rn ) <= tol * rn )) 
        {
          break;
        }
        value = rn;
      }

      pn1 = pn3;
      pn2 = pn4;
      pn3 = pn5;
      pn4 = pn6;
//
//  Re-scale terms in continued fraction if terms are large.
//
      if ( oflo <= abs ( pn5 ) )
      {
        pn1 = pn1 / oflo;
        pn2 = pn2 / oflo;
        pn3 = pn3 / oflo;
        pn4 = pn4 / oflo;
      }
    }

    arg = arg + log ( value );

    if ( elimit <= arg )
    {
      value = 1.0 - exp ( arg );
    }
    else
    {
      value = 1.0;
    }
  }

  return 1. - value;

}


