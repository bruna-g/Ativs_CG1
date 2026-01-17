#include "../include/Matriz3x3.h"
#include "../include/Vector.h"

Matriz3x3::Matriz3x3(float m00, float m01, float m02,
                     float m10, float m11, float m12,
                     float m20, float m21, float m22) {
    m[0][0] = m00; m[0][1] = m01; m[0][2] = m02;
    m[1][0] = m10; m[1][1] = m11; m[1][2] = m12;
    m[2][0] = m20; m[2][1] = m21; m[2][2] = m22;
}

Matriz3x3 outerProduto(const Vector& d) {
    Matriz3x3 M(0.0f, 0.0f, 0.0f,
                0.0f, 0.0f, 0.0f,
                0.0f, 0.0f, 0.0f);
    M.m[0][0] = d.i * d.i;
    M.m[0][1] = d.i * d.j;
    M.m[0][2] = d.i * d.k;

    M.m[1][0] = d.j * d.i;
    M.m[1][1] = d.j * d.j;
    M.m[1][2] = d.j * d.k;

    M.m[2][0] = d.k * d.i;
    M.m[2][1] = d.k * d.j;
    M.m[2][2] = d.k * d.k;
    return M;
}

Matriz3x3 matrizSubtrai(const Matriz3x3& A, const Matriz3x3& B) {
    Matriz3x3 R(0.0f, 0.0f, 0.0f,
                0.0f, 0.0f, 0.0f,
                0.0f, 0.0f, 0.0f);
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c)
            R.m[r][c] = A.m[r][c] - B.m[r][c];
    return R;
}

Matriz3x3 escalarMatriz(float s, const Matriz3x3& M) {
    Matriz3x3 R(0.0f, 0.0f, 0.0f,
                0.0f, 0.0f, 0.0f,
                0.0f, 0.0f, 0.0f);
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c)
            R.m[r][c] = s * M.m[r][c];
    return R;
}

Matriz3x3 transpostaVetor(const Vector& v) {
    Matriz3x3 M(0.0f, 0.0f, 0.0f,
                0.0f, 0.0f, 0.0f,
                0.0f, 0.0f, 0.0f);
    M.m[0][0] = v.i;
    M.m[1][0] = v.j;
    M.m[2][0] = v.k;
    return M;
}

float transpostoPorVetor(const Matriz3x3& M, const Vector& v) {
    float res = M.m[0][0] * v.i + M.m[1][0] * v.j + M.m[2][0] * v.k;
    return res;
}
