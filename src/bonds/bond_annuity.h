//
// Created by Jonas Wolff on 15/03/2024.
//
#pragma once
#ifndef LIBRARYCORE_BOND_ANNUITY_H
#define LIBRARYCORE_BOND_ANNUITY_H


#include <stdexcept>
#include <functional>
#include "math.h"


using namespace std;

/****************************************************************************/
//                Class representing an annuity                             //
// One needs to all characteristics but one for the class to determine      //
// the type of annuity.                                                     //
/****************************************************************************/
class bond_annuity {
private:

    // The following characteristics are respectively the interest rate R,
    // the debt in time 0 value, the fixed payment each year C, and the
    // number of years to maturity.
    double R;
    double D0;
    double C;
    double n;

public:

    /*******************************************************/
    //                    CONSTRUCTORS                     //
    // One can give full parametrization of the annuity or //
    // One can also leave one parameter to be determined.  //
    // When leaving one parameter to be determined it is   //
    // set to NULL as a void pointer                       //
    /*******************************************************/
    bond_annuity(double &pR, double &pD0, double &pC, double &pn);

    /*************************************************************/
    //              No interest rate is given                    //
    /*************************************************************/
    bond_annuity(void *ptrR, double &pD0, double &pC, double &pn);

    /*************************************************************/
    //               No principal amount is given                //
    /*************************************************************/
    bond_annuity(double &pR, void *ptrD0, double &pC, double &pn);

    /*************************************************************/
    //               No payment size is given                    //
    /*************************************************************/
    bond_annuity(double &pR, double &pD0, void *ptrC, double &pn);

    /*************************************************************/
    //               No loan maturity is given                   //
    /*************************************************************/
    bond_annuity(double &pR, double &pD0, double &pC, void *ptrn);

    // Destructor
    ~bond_annuity();

    /*************************************************************/
    //       Get back a vector of the annuity's parameters       //
    /*************************************************************/
    vector<double> getparams();





};


#endif //LIBRARYCORE_BOND_ANNUITY_H
