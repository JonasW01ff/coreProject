
#include "functions.h"
#include <cmath>
#include "math.h"
#include "functions_math.h"
#include <functional>
#include <numeric>
#include "algorithm"

using namespace std;

double alphahage(double &R, double &n){
    return (1-pow(1+R,-n))/R;
}

double AnnuityRate(double pD0, double pC, double pn){
    // Find the annuity's interest rate through bisection
    double R;
    function<double(double)> func = [=](double R)  {return pD0/pC-(1-pow(1+R,-pn))/R;};
    double eps = 1e-10;
    functions_math::func_solve_root_bisection(func,-0.9, 1., R, eps);
    return R;
}

double AnnuityLength(double pD0, double pC, double pR){
    return -log(1-pD0*pR/pC)/log(1+pR);
}