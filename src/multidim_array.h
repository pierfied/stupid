//
// Created by pierfied on 11/19/18.
//

#ifndef STUPID_MULTIDIM_ARRAY_H
#define STUPID_MULTIDIM_ARRAY_H

#include <stdlib.h>

template<typename T>
class array_2d {
    T *data;

public:
    const int nx;
    const int ny;

    array_2d(int nx, int ny) : nx(nx), ny(ny) {
        data = new T[nx * ny]();
    }

    inline T &operator()(int i, int j) {
        return data[i * ny + j];
    }

    inline T &index(int i, int j) {
        return operator()(i, j);
    }

    void reset_zero() {
        for (int i = 0; i < nx * ny; ++i) {
            data[i] = T();
        }
    }

    ~array_2d() {
        delete[] data;
    }
};


template<typename T>
class array_3d {
    T *data;

public:
    const int nx;
    const int ny;
    const int nz;

    array_3d(int nx, int ny, int nz) : nx(nx), ny(ny), nz(nz) {
        data = new T[nx * ny * nz]();
    }

    inline T &operator()(int i, int j, int k) {
        return data[(i * ny + j) * nz + k];
    }

    inline T &index(int i, int j, int k) {
        return operator()(i, j, k);
    }

    void reset_zero() {
        for (int i = 0; i < nx * ny * nz; ++i) {
            data[i] = T();
        }
    }

    ~array_3d() {
        delete[] data;
    }
};


#endif //STUPID_MULTIDIM_ARRAY_H
