//
// Created by pierfied on 11/21/18.
//

#include <gtest/gtest.h>

#include "test_multidim_array.cpp"
#include "test_cic_grid.cpp"
#include "test_leapfrog_integrator.cpp"

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}