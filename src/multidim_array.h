//
// Created by pierfied on 11/19/18.
//

#ifndef STUPID_MULTIDIM_ARRAY_H
#define STUPID_MULTIDIM_ARRAY_H

#include <stdlib.h>

template<typename T>
class array_2d {
public:
    T *data;

    const long nx;
    const long ny;

    array_2d(long nx, long ny) : nx(nx), ny(ny) {
        data = new T[nx * ny]();
    }

    inline T &operator()(long i, long j) {
        return data[i * ny + j];
    }

    inline T &index(long i, long j) {
        return operator()(i, j);
    }

    void reset_zero() {
        for (long i = 0; i < nx * ny; ++i) {
            data[i] = T();
        }
    }

    ~array_2d() {
        delete[] data;
    }
};


template<typename T>
class array_3d {
public:
    T *data;

    const long nx;
    const long ny;
    const long nz;

    array_3d(long nx, long ny, long nz) : nx(nx), ny(ny), nz(nz) {
        data = new T[nx * ny * nz]();
    }

    inline T &operator()(long i, long j, long k) {
        return data[(i * ny + j) * nz + k];
    }

    inline T &index(long i, long j, long k) {
        return operator()(i, j, k);
    }

    void reset_zero() {
        for (long i = 0; i < nx * ny * nz; ++i) {
            data[i] = T();
        }
    }

    ~array_3d() {
        delete[] data;
    }
};


#endif //STUPID_MULTIDIM_ARRAY_H
