cmake_minimum_required(VERSION 3.26)
project(coreDLL)

set(CMAKE_CXX_STANDARD 23)

add_library(coreDLL SHARED functions_export.cpp
        functions/functions.cpp
        functions/functions_math.cpp
        algebra/linear_algebra/matrix.cpp
        algorithms/cacm222.cpp
        algorithms/toms179.cpp
        bonds/bond_annuity.cpp
        risk/risk.cpp
        statistics/cdflib.cpp
        statistics/distribution_binomial.cpp
        statistics/distribution_gamma.cpp
        statistics/distribution_normal.cpp)
