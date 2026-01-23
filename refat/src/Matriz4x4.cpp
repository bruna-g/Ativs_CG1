#include "../include/Matriz4x4.h"
#include "../include/Utils.h"
#include <cmath>

#define M_PI 3.14159265358979323846

Matriz4x4::Matriz4x4(float valores[4][4]) {
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            m[i][j] = valores[i][j];
        }
    }
}

Matriz4x4 Matriz4x4::rotacaoX(float angulo) {
    float rad = angulo * M_PI / 180.0f;
    float c = cos(rad);
    float s = sin(rad);
    float rot[4][4] = {
        {1, 0, 0, 0},
        {0, c, -s, 0},
        {0, s, c, 0},
        {0, 0, 0, 1}
    };
    return Matriz4x4(rot);
}

Matriz4x4 Matriz4x4::rotacaoY(float angulo) {
    float rad = angulo * M_PI / 180.0f;
    float c = cos(rad);
    float s = sin(rad);
    float rot[4][4] = {
        {c, 0, s, 0},
        {0, 1, 0, 0},
        {-s, 0, c, 0},
        {0, 0, 0, 1}
    };
    return Matriz4x4(rot);
}

Matriz4x4 Matriz4x4::rotacaoZ(float angulo) {
    float rad = angulo * M_PI / 180.0f;
    float c = cos(rad);
    float s = sin(rad);
    float rot[4][4] = {
        {c, -s, 0, 0},
        {s, c, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    };
    return Matriz4x4(rot);
}

Matriz4x4 Matriz4x4::rotacaoArbitraria(Vector eixo, float angulo) {
    Vector n = normalizar(eixo);
    float rad = angulo * M_PI / 180.0f;
    float c = cos(rad);
    float s = sin(rad);
    float rot[4][4] = {
        {c + n.i*n.i*(1-c), n.i*n.j*(1-c) - n.k*s, n.i*n.k*(1-c) + n.j*s, 0},
        {n.j*n.i*(1-c) + n.k*s, c + n.j*n.j*(1-c), n.j*n.k*(1-c) - n.i*s, 0},
        {n.k*n.i*(1-c) - n.j*s, n.k*n.j*(1-c) + n.i*s, c + n.k*n.k*(1-c), 0},
        {0, 0, 0, 1}
    };
    return Matriz4x4(rot);
}

Point Matriz4x4::multiplicarPonto(const Point& p) const {
    float x = m[0][0]*p.x + m[0][1]*p.y + m[0][2]*p.z + m[0][3];
    float y = m[1][0]*p.x + m[1][1]*p.y + m[1][2]*p.z + m[1][3];
    float z = m[2][0]*p.x + m[2][1]*p.y + m[2][2]*p.z + m[2][3];
    float w = m[3][0]*p.x + m[3][1]*p.y + m[3][2]*p.z + m[3][3];
    if (w != 1.0f) {
        x /= w;
        y /= w;
        z /= w;
    }
    return Point(x, y, z);
}