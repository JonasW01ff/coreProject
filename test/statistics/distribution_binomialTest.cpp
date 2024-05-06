//
// Created by Jonas Wolff on 02/05/2024.
//

#include "gtest/gtest.h"
#include "statistics/distribution_binomial.h"
#include "math.h"

using namespace std;

TEST(statistics_test, distribution_binomial_test_1){
    int n,k;
    double p;
    double abs_error = 0.00001;
    distribution_binomial exampledist = distribution_binomial(n, p);
    n = 163, k = 41, p= 0.536076355658045;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(exampledist.cdf(k), 1.35951568641937E-13, abs_error);
    n = 440, k = 121, p= 0.0890400037997268;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 1,  abs_error);
    n = 293, k = 15, p= 0.588787251401285;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 1.01089301226686E-86,  abs_error);
    n = 735, k = 98, p= 0.56491038226123;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 3.1874697853076E-131,  abs_error);
    n = 150, k = 95, p= 0.430536741785791;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 0.999999805664099,  abs_error);
    n = 308, k = 40, p= 0.1249846856892;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 0.642168653022687,  abs_error);
    n = 98, k = 96, p= 0.760118358585749;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 0.99999999993231,  abs_error);
    n = 987, k = 188, p= 0.162400333671768;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 0.991591048341306,  abs_error);
    n = 211, k = 137, p= 0.49255997079321;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 0.999998371396081,  abs_error);
    n = 963, k = 281, p= 0.0468627511457321;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 1,  abs_error);
    n = 925, k = 878, p= 0.353762367004671;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 1,  abs_error);
    n = 149, k = 92, p= 0.296606549150753;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 1,  abs_error);
    n = 882, k = 54, p= 0.0896405252663913;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 0.00120259251562635,  abs_error);
    n = 25, k = 22, p= 0.23323349691393;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 0.999999999999479,  abs_error);
    n = 414, k = 214, p= 0.335286493607483;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 0.999999999999989,  abs_error);
    n = 668, k = 275, p= 0.969377120187492;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 0,  abs_error);
    n = 712, k = 195, p= 0.857857382609225;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 1.1533104333113E-271,  abs_error);
    n = 7, k = 3, p= 0.294065396226026;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 0.881527695861613,  abs_error);
    n = 3, k = 1, p= 0.5750643716944;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 0.388249366868291,  abs_error);
    n = 9, k = 3, p= 0.65538840494824;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 0.049772228413173,  abs_error);
    n = 9, k = 2, p= 0.454776609116707;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 0.142859828609298,  abs_error);
    n = 7, k = 7, p= 0.854169274979956;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 1,  abs_error);
    n = 10, k = 1, p= 0.588291179400762;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 0.00213934636701296,  abs_error);
    n = 5, k = 5, p= 0.214678492807644;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 1,  abs_error);
    n = 8, k = 7, p= 0.694511329588108;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 0.945870417382839,  abs_error);
    n = 11, k = 0, p= 0.476917665290616;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 0.000802177706700666,  abs_error);
    n = 6000, k = 4853, p= 0.514046198770022;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 1,  abs_error);
    n = 6507, k = 899, p= 0.613560236409762;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 0,  abs_error);
    n = 6171, k = 6120, p= 0.855614005756625;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 1,  abs_error);
    n = 1388, k = 574, p= 0.292105029723995;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 1,  abs_error);
    n = 5603, k = 3273, p= 0.442106858345011;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 1,  abs_error);
    n = 6817, k = 52, p= 0.135426844864049;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 0,  abs_error);
    n = 2354, k = 2177, p= 0.145462384685106;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 1,  abs_error);
    n = 6483, k = 6130, p= 0.920540865808698;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 0.999999999999998,  abs_error);
    n = 5385, k = 146, p= 0.0681008938297409;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(   exampledist.cdf(k), 1.55928162828391E-41,  abs_error);
    n = 7873, k = 1872, p= 0.86835995993961;
    exampledist = distribution_binomial(n, p);
    ASSERT_NEAR(exampledist.cdf(k), 0,  abs_error);


}


