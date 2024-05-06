//
// Created by Jonas Wolff on 28/02/2024.
//
#pragma once
#ifndef LIBRARYCORE_DISTRIBUTION_NORMAL_H
#define LIBRARYCORE_DISTRIBUTION_NORMAL_H

// class object for the normal distribution function
class distribution_normal {
private:
    // mean
    double mu;
    // standard deviation
    double sigma;

public:

    // Constructor
    distribution_normal(double pmu, double psigma);

    // Destructor
    ~distribution_normal();

    // PDF
    double pdf(double &x);

    // CDF
    double cdf(double &x);

};


#endif //LIBRARYCORE_DISTRIBUTION_NORMAL_H
