//
// Created by pierfied on 11/21/18.
//

#include <stupid.h>
#include <gtest/gtest.h>

#include "test_multidim_array.cpp"

int main(int argc, char **argv) {
    hello();
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}