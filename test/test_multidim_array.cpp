//
// Created by pierfied on 11/21/18.
//

#include <gtest/gtest.h>
#include <multidim_array.h>

TEST(multidim_array, dims_2d) {
    int n = 10;

    array_2d<double> arr(n, n);

    EXPECT_EQ(n, arr.nx);
    EXPECT_EQ(n, arr.ny);
}

TEST(multidim_array, write_read_2d) {
    int n = 10;

    array_2d<double> arr(n, n);

    int count = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            arr(i, j) = count++;
        }
    }

    count = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            EXPECT_EQ(count++, arr(i, j));
        }
    }
}

TEST(multidim_array, zero_2d) {
    int n = 10;

    array_2d<double> arr(n, n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            arr(i, j) = 1;
        }
    }

    arr.reset_zero();

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            EXPECT_EQ(0, arr(i, j));
        }
    }
}

TEST(multidim_array, dims_3d) {
    int n = 10;

    array_3d<double> arr(n, n, n);

    EXPECT_EQ(n, arr.nx);
    EXPECT_EQ(n, arr.ny);
    EXPECT_EQ(n, arr.nz);
}

TEST(multidim_array, write_read_3d) {
    int n = 10;

    array_3d<double> arr(n, n, n);

    int count = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                arr(i, j, k) = count++;
            }
        }
    }

    count = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                EXPECT_EQ(count++, arr(i, j, k));
            }
        }
    }
}

TEST(multidim_array, zero_3d) {
    int n = 10;

    array_3d<double> arr(n, n, n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                arr(i, j, k) = 1;
            }
        }
    }

    arr.reset_zero();

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                EXPECT_EQ(0, arr(i, j, k));
            }
        }
    }
}