//
// Created by Jonas Wolff on 26/02/2024.
//
#pragma once
#ifndef LIBRARYCORE_RISK_H
#define LIBRARYCORE_RISK_H

#include <vector>

using namespace std;

class risk {
private:
    vector<double> losses;

public:

    // Constructor
    risk();

    // Overloaded Constructor
    risk(vector<double> losses);

    // Destructor
    ~risk();

    // Load losses
    void loadLosses(vector<double> losses);

    // Calculate Empirical VaR
    double calcEmpiricalVaR(double alpha);

    // Calculate exact confidence intervals for Empirical VaR
    vector<double> calcCiEmpiricalVaR(double alpha, double probability);

};

#endif //LIBRARYCORE_RISK_H
