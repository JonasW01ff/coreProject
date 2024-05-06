//
// Created by Jonas Wolff on 15/03/2024.
//
#pragma once
#ifndef LIBRARYCORE_TOMS179_H
#define LIBRARYCORE_TOMS179_H

double alogam ( double x, int *ifault );
void beta_cdf_values ( int *n_data, double *a, double *b, double *x,
                       double *fx );
void gamma_log_values ( int *n_data, double *x, double *fx );
double mdbeta ( double x, double p, double q, int *ier );


#endif //LIBRARYCORE_TOMS179_H
