//
// Created by Jonas Wolff on 28/02/2024.
//
#pragma once
#ifndef LIBRARYCORE_DISTRIBUTION_BINOMIAL_H
#define LIBRARYCORE_DISTRIBUTION_BINOMIAL_H


// class object for the binomial distribution function
class distribution_binomial {
private:

    // trials
    int n;

    // Success probability
    double p;

public:

    // Constructor
    distribution_binomial(int pn, double pp);

    // Destructor
    ~distribution_binomial();

    // PDF
    double pdf(int &k);

    // CDF
    double cdf(int &k);

};
#endif //LIBRARYCORE_DISTRIBUTION_BINOMIAL_H