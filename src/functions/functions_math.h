//
// Created by Jonas Wolff on 28/02/2024.
//
#pragma once
#ifndef LIBRARYCORE_FUNCTIONS_MATH_H
#define LIBRARYCORE_FUNCTIONS_MATH_H

#endif //LIBRARYCORE_FUNCTIONS_MATH_H

#include <functional>

using namespace std;

namespace functions_math{

    // Calculates the binomial coefficient
    double binomial_coeficient(int n, int k);

    // Calculate the incomplete beta function
    double func_beta_incomplete(double a, double b, double x, double &eps);

    // Calculates the beta function
    double func_beta(double &a, double &b);

    // Calculates the incomplete gamma function
    double func_gamma_incomplete(double &a, double &x);

    // Find root of function by bisection
    void func_solve_root_bisection(std::function<double(double)> func,double xL, double xU, double &xN, double &eps);

    // Find root of function by secant method
    void func_solve_root_secant(std::function<double(double)> func,double x1, double x0, double &xres, double &eps,  double iterationscale = 1.0);

    double func_digamma(double x);

    void func_minimize_NM_2D(function<double(double, double)> func, double x0, double y0, double &x1, double &y1, double &eps);

    void func_fit_histogram_gamma(vector<double> bins, vector<double> freqs, double &alpha, double &beta);

}
