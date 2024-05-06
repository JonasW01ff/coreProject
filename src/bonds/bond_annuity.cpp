
#include "bond_annuity.h"
#include <stdexcept>
#include <functional>
#include "math.h"

#include "../functions/functions_math.h"

bond_annuity::bond_annuity(double &pR, double &pD0, double &pC, double &pn){
    
        // Save all parameters
        R=pR;
        D0=pD0;
        C=pC;
        n=pn;
}

bond_annuity::bond_annuity(void *ptrR, double &pD0, double &pC, double &pn){

        // Find the annuity's interest rate through bisection
        double pR;
        function<double(double)> func = [=](double R)  {return pD0/pC-(1-pow(1+R,-pn))/R;};
        double eps = 1e-10;
        functions_math::func_solve_root_bisection(func,-0.9, 1., pR, eps);

        // Save all parameters
        R=pR;
        D0=pD0;
        C=pC;
        n=pn;
    }

bond_annuity::~bond_annuity(){
    // Destructor
}


bond_annuity::bond_annuity(double &pR, void *ptrD0, double &pC, double &pn){

        // Find the principal amount by use of alpha hage
        double pD0;
        pD0=pC*(1-pow(1+pR,-pn))/pR;


        // Save all parameters
        R=pR;
        D0=pD0;
        C=pC;
        n=pn;
    }


bond_annuity::bond_annuity(double &pR, double &pD0, void *ptrC, double &pn){

        // Find the payment size by use of alpha hage
        double pC;
        pC = pD0*pR/(1-pow(1+pR,-pn));

        // Save all parameters
        R=pR;
        D0=pD0;
        C=pC;
        n=pn;
    }

bond_annuity::bond_annuity(double &pR, double &pD0, double &pC, void *ptrn){

        // Find the loan maturity by use of alpha hage.
        double pn;
        pn=-log(1-pD0*pR/pC)/log(1+pR);

        // Check if loan is paid of in whole years.
        if (pn - round(pn) > 0.01){
            throw runtime_error("The annuity does not run for a whole number of years.");
        }

        // Save all parameters.
        R=pR;
        D0=pD0;
        C=pC;
        n=pn;
    }

vector<double> bond_annuity::getparams(){
            vector<double> res;
            res.push_back(R);
            res.push_back(D0);
            res.push_back(C);
            res.push_back(n);
            return res;
        }
