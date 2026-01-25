#include "../include/Matriz4x4.h"
#include "../include/Utils.h"
#include <cmath>

#define M_PI 3.14159265358979323846

Matriz4x4::Matriz4x4() {
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            m[i][j] = (i == j) ? 1.0f : 0.0f;
        }
    }
}

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
    Vector eixo_norm = normalizar(eixo);
    Vector temp(1.f, 0.f, 0.f);
    if (std::fabs(eixo_norm.i > 0.9f)) {
        temp = Vector(0.f, 1.f, 0.f);
    }

    Vector local_x = normalizar(produto_vetorial(eixo_norm, temp));
    Vector local_y = normalizar(produto_vetorial(eixo_norm, local_x));

    float M_data[4][4] = {
        {local_x.i, local_x.j, local_x.k, 0},
        {local_y.i, local_y.j, local_y.k, 0},
        {eixo_norm.i, eixo_norm.j, eixo_norm.k, 0},
        {0, 0, 0, 1}
    };
    Matriz4x4 M(M_data);

    Matriz4x4 M_inversa = inversa(M);
    Matriz4x4 R = rotacaoZ(angulo);

    Matriz4x4 M_inv_R = M_inversa.multiplicarMat(R);
    Matriz4x4 resultado = M_inv_R.multiplicarMat(M);

    return resultado;
}

Matriz4x4 Matriz4x4::multiplicarMat(Matriz4x4& outra) const {
    Matriz4x4 resultado;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            resultado.m[i][j] = 0.0f;
            for (int k = 0; k < 4; ++k) {
                resultado.m[i][j] += m[i][k] * outra.m[k][j];
            }
        }
    }
    return resultado;
}

Point Matriz4x4::multiplicarPonto(const Point& p) const {
    float x = m[0][0] * p.x + m[0][1] * p.y + m[0][2] * p.z + m[0][3];
    float y = m[1][0] * p.x + m[1][1] * p.y + m[1][2] * p.z + m[1][3];
    float z = m[2][0] * p.x + m[2][1] * p.y + m[2][2] * p.z + m[2][3];
    float w = m[3][0] * p.x + m[3][1] * p.y + m[3][2] * p.z + m[3][3];
    if (w != 1.0f) {
        x /= w;
        y /= w;
        z /= w;
    }
    return Point(x, y, z);
}

Vector Matriz4x4::multiplicarVetor(const Vector& v) const {
    float i = m[0][0] * v.i + m[0][1] * v.j + m[0][2] * v.k;
    float j = m[1][0] * v.i + m[1][1] * v.j + m[1][2] * v.k;
    float k = m[2][0] * v.i + m[2][1] * v.j + m[2][2] * v.k;
    return Vector(i, j, k);
}


Matriz4x4 Matriz4x4::identidade() {
    float id[4][4] = {
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    };
    return Matriz4x4(id);
}

Matriz4x4 Matriz4x4::inversa(const Matriz4x4& matriz) {
    Matriz4x4 inv;

    float s0 = matriz.m[0][0] * matriz.m[1][1] - matriz.m[1][0] * matriz.m[0][1];
    float s1 = matriz.m[0][0] * matriz.m[1][2] - matriz.m[1][0] * matriz.m[0][2];
    float s2 = matriz.m[0][0] * matriz.m[1][3] - matriz.m[1][0] * matriz.m[0][3];
    float s3 = matriz.m[0][1] * matriz.m[1][2] - matriz.m[1][1] * matriz.m[0][2];
    float s4 = matriz.m[0][1] * matriz.m[1][3] - matriz.m[1][1] * matriz.m[0][3];
    float s5 = matriz.m[0][2] * matriz.m[1][3] - matriz.m[1][2] * matriz.m[0][3];

    float c5 = matriz.m[2][2] * matriz.m[3][3] - matriz.m[3][2] * matriz.m[2][3];
    float c4 = matriz.m[2][1] * matriz.m[3][3] - matriz.m[3][1] * matriz.m[2][3];
    float c3 = matriz.m[2][1] * matriz.m[3][2] - matriz.m[3][1] * matriz.m[2][2];
    float c2 = matriz.m[2][0] * matriz.m[3][3] - matriz.m[3][0] * matriz.m[2][3];
    float c1 = matriz.m[2][0] * matriz.m[3][2] - matriz.m[3][0] * matriz.m[2][2];
    float c0 = matriz.m[2][0] * matriz.m[3][1] - matriz.m[3][0] * matriz.m[2][1];

    float det =
        s0 * c5 - s1 * c4 + s2 * c3
        + s3 * c2 - s4 * c1 + s5 * c0;

    if (fabs(det) < 1e-6f) {
        return Matriz4x4::identidade();
    }

    float invDet = 1.0f / det;

    inv.m[0][0] = (matriz.m[1][1] * c5 - matriz.m[1][2] * c4 + matriz.m[1][3] * c3) * invDet;
    inv.m[0][1] = (-matriz.m[0][1] * c5 + matriz.m[0][2] * c4 - matriz.m[0][3] * c3) * invDet;
    inv.m[0][2] = (matriz.m[3][1] * s5 - matriz.m[3][2] * s4 + matriz.m[3][3] * s3) * invDet;
    inv.m[0][3] = (-matriz.m[2][1] * s5 + matriz.m[2][2] * s4 - matriz.m[2][3] * s3) * invDet;

    inv.m[1][0] = (-matriz.m[1][0] * c5 + matriz.m[1][2] * c2 - matriz.m[1][3] * c1) * invDet;
    inv.m[1][1] = (matriz.m[0][0] * c5 - matriz.m[0][2] * c2 + matriz.m[0][3] * c1) * invDet;
    inv.m[1][2] = (-matriz.m[3][0] * s5 + matriz.m[3][2] * s2 - matriz.m[3][3] * s1) * invDet;
    inv.m[1][3] = (matriz.m[2][0] * s5 - matriz.m[2][2] * s2 + matriz.m[2][3] * s1) * invDet;

    inv.m[2][0] = (matriz.m[1][0] * c4 - matriz.m[1][1] * c2 + matriz.m[1][3] * c0) * invDet;
    inv.m[2][1] = (-matriz.m[0][0] * c4 + matriz.m[0][1] * c2 - matriz.m[0][3] * c0) * invDet;
    inv.m[2][2] = (matriz.m[3][0] * s4 - matriz.m[3][1] * s2 + matriz.m[3][3] * s0) * invDet;
    inv.m[2][3] = (-matriz.m[2][0] * s4 + matriz.m[2][1] * s2 - matriz.m[2][3] * s0) * invDet;

    inv.m[3][0] = (-matriz.m[1][0] * c3 + matriz.m[1][1] * c1 - matriz.m[1][2] * c0) * invDet;
    inv.m[3][1] = (matriz.m[0][0] * c3 - matriz.m[0][1] * c1 + matriz.m[0][2] * c0) * invDet;
    inv.m[3][2] = (-matriz.m[3][0] * s3 + matriz.m[3][1] * s1 - matriz.m[3][2] * s0) * invDet;
    inv.m[3][3] = (matriz.m[2][0] * s3 - matriz.m[2][1] * s1 + matriz.m[2][2] * s0) * invDet;

    return inv;
}

Matriz4x4 Matriz4x4::transposta(const Matriz4x4& matriz) {
    Matriz4x4 transposta;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            transposta.m[i][j] = matriz.m[j][i];
        }
    }
    return transposta;
}