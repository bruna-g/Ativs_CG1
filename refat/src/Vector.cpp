#include "../include/Vector.h"
#include "../include/Point.h"
#include "../include/Cilindro.h"
#include "../include/Matriz3x3.h"
#include <cmath>

Vector::Vector(float i_v, float j_v, float k_v) {
    this->i = i_v;
    this->j = j_v;
    this->k = k_v;
    this->q = 0;
}

Vector cross(const Vector& a, const Vector& b) {
    return Vector(a.j * b.k - a.k * b.j,
                  a.k * b.i - a.i * b.k,
                  a.i * b.j - a.j * b.i);
}

Vector matrizVetorProduto(const Matriz3x3& M, const Vector& v) {
    float i_res = M.m[0][0] * v.i + M.m[0][1] * v.j + M.m[0][2] * v.k;
    float j_res = M.m[1][0] * v.i + M.m[1][1] * v.j + M.m[1][2] * v.k;
    float k_res = M.m[2][0] * v.i + M.m[2][1] * v.j + M.m[2][2] * v.k;
    return Vector(i_res, j_res, k_res);
}

Vector matrizPorMatriz(const Matriz3x3& A, const Matriz3x3& B) {
    float i_res = A.m[0][0] * B.m[0][0] + A.m[0][1] * B.m[1][0] + A.m[0][2] * B.m[2][0];
    float j_res = A.m[1][0] * B.m[0][0] + A.m[1][1] * B.m[1][0] + A.m[1][2] * B.m[2][0];
    float k_res = A.m[2][0] * B.m[0][0] + A.m[2][1] * B.m[1][0] + A.m[2][2] * B.m[2][0];
    return Vector(i_res, j_res, k_res);
}

