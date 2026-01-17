#ifndef ALGEBRALINEAR_H
#define ALGEBRALINEAR_H

#include "Vector.h"
#include <cmath>

class AlgebraLinear {
public:
    AlgebraLinear() = default;

    Vector somar(const Vector& a, const Vector& b) const {
        return Vector(a.i + b.i, a.j + b.j, a.k + b.k);
    }

    Vector subtrair(const Vector& a, const Vector& b) const {
        return Vector(a.i - b.i, a.j - b.j, a.k - b.k);
    }

    Vector escalar(float s, const Vector& v) const {
        return Vector(s * v.i, s * v.j, s * v.k);
    }

    float produtoEscalar(const Vector& a, const Vector& b) const {
        return a.i * b.i + a.j * b.j + a.k * b.k;
    }

    Vector produtoVetorial(const Vector& a, const Vector& b) const {
        return Vector(
            a.j * b.k - a.k * b.j,
            a.k * b.i - a.i * b.k,
            a.i * b.j - a.j * b.i
        );
    }

    float norma(const Vector& v) const {
        return std::sqrt(produtoEscalar(v, v));
    }

    Vector normalizar(const Vector& v) const {
        float n = norma(v);
        if (n <= 1e-8f) return Vector(0.0f, 0.0f, 0.0f);
        return escalar(1.0f / n, v);
    }
};

#endif // ALGEBRALINEAR_H
