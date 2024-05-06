/* File: Ratings.c
/* Compile with -fvisibility=hidden.
/* To use it in excel place the .dylib file in /Users/jonaswolff/Library/Containers/com.microsoft.Excel/Data/tmp
/* The below example show how to use the functions in c++
/*     // Open the library.
/*    char *lib_name = "./bin/libCORE.dylib";
/*    void *lib_handle = dlopen(lib_name, RTLD_NOW);
/*
/*    // Get the symbol addresses.
/*    typedef void (*hello_t)(char*);
/*    hello_t addRating = (hello_t) dlsym(lib_handle, "addRating");
/*    if (!lib_handle) {
/*        cerr << "Cannot open library: " << dlerror() << '\n';
/*        return 1;
/*    }
/*    char * arg1 = "AA";
/*    addRating(arg1);
/*
/*    // Get the symbol addresses.
/*    typedef void (*hello_k)(void);
/*    hello_k clearRatings = (hello_k) dlsym(lib_handle, "clearRatings");
/*    clearRatings();
/*    typedef int (*hello_p)(void);
/*    hello_p ratings = (hello_p) dlsym(lib_handle, "ratings");
/*    int x = ratings();
/*    cout << x << endl;
/*
/*    
/*
/*    dlclose(lib_handle);
 **********************************/

#include "functions_export.h"
#include "functions/functions_math.h"
#include "functions/functions.h"
#include "statistics/distribution_binomial.h"
#include "statistics/distribution_gamma.h"
#include "bonds/bond_annuity.h"
#include "risk/risk.h"
#include <stdio.h>
#include <string.h>
 
#define EXPORT __attribute__((visibility("default")))
#define MAX_NUMBERS 99
 
#define function_completed 100
#define function_failed -1

static int _number_list[MAX_NUMBERS];
static int _numbers = 0;
 
// Initializer.
__attribute__((constructor))
static void initializer(void) {                             // 2
    printf("[%s] initializer()\n", __FILE__);
}
 
// Finalizer.
__attribute__((destructor))
static void finalizer(void) {                               // 3
    printf("[%s] finalizer()\n", __FILE__);
}
 
// Used by meanRating, middleRating, frequentRating.
static char *_char_rating(int rating) {
    char result[10] = "";
    int int_rating = rating;
    for (int i = 0; i < int_rating; i++) {
        strncat(result, "*", sizeof(result) - strlen(result) - 1);
    }
    return strdup(result);
}
 
// Used by addRating.
void _add(int number) {                                     // 4
    if (_numbers < MAX_NUMBERS) {
        _number_list[_numbers++] = number;
    }
}
 
// Used by meanRating.
int _mean(void) {
    int result = 0;
    if (_numbers) {
        int sum = 0;
        int i;
        for (i = 0; i < _numbers; i++) {
            sum += _number_list[i];
        }
        result = sum / _numbers;
    }
    return result;
}
 
EXPORT
void addRating(char *rating) {                            // 5
    if (rating != NULL) {
        int numeric_rating = 0;
        int pos = 0;
        while (*rating++ != '\0' && pos++ < 5) {
            numeric_rating++;
        }
        _add(numeric_rating);
    }
}
 
EXPORT
char *meanRating(void) {
    return _char_rating(_mean());
}
 
EXPORT
int ratings(void) {
    return _numbers;
}
 
EXPORT
void clearRatings(void) {
    _numbers = 0;
}

EXPORT
double xIncomepleteBeta(double a, double b, double x){
    double eps = 1e-16;
    try {
        return functions_math::func_beta_incomplete(a, b, x, eps);
    } catch (const std::runtime_error& e) {
        return -1.;
    }
}

EXPORT
double xBinomialCDF(int k, int n, double p){
    try {
        distribution_binomial myDist(n,p);
        return myDist.cdf(k);
    } catch (const std::runtime_error& e) {
        return -1.;
    }

}

EXPORT 
int xGammaMLE(double *bins, short bins_size, int *freqs, int freqs_size, double *alpha, double *beta){
    try {
        vector<double> bins_vec;
        vector<int> freqs_vec;
        for (int i=0; i<bins_size; i++){
            bins_vec.push_back(bins[i]);
        }
        for (int i=0; i<freqs_size; i++){
            freqs_vec.push_back(freqs[i]);
        }
        distribution_gamma myDist(bins_vec, freqs_vec);
        *alpha = myDist.getAlpha();
        *beta = myDist.getBeta();
        return function_completed;
    } catch (const std::runtime_error& e) {
        return function_failed;
    }
}

EXPORT
double xAlphaHage(double R, double n){
    return alphahage(R, n);
}

EXPORT
double xAnnuityRate(double pD0, double pC, double pn){
    return AnnuityRate(pD0, pC, pn);
}

EXPORT
double xAnnuityLength(double pD0, double pC, double pR){
    return AnnuityLength(pD0, pC, pR);
}

EXPORT
int xVaRCI(double *losses, int n, double alpha, double prob, double *res){

    try {
        vector<double> losses_vec;
        for (int i=0; i<n; i++){
            losses_vec.push_back(losses[i]);
        }

        risk myrisk(losses_vec);
        vector<double> ci = myrisk.calcCiEmpiricalVaR(alpha, prob);

        static const size_t result_size = 2;
        res[0] = ci[0];
        res[1] = ci[1];
        return function_completed;
    } catch (const std::runtime_error& e) {
        return function_failed;
    }

};
