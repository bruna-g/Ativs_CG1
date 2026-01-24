#ifndef MATRIZ4X4_H
#define MATRIZ4X4_H

#include "Vector.h"
#include "Point.h"

class Matriz4x4 {
public:
    float m[4][4];

    Matriz4x4();
    Matriz4x4(float valores[4][4]);

    static Matriz4x4 rotacaoX(float angulo);
    static Matriz4x4 rotacaoY(float angulo);
    static Matriz4x4 rotacaoZ(float angulo);
    static Matriz4x4 rotacaoArbitraria(Vector eixo, float angulo);
    static Matriz4x4 identidade();
    static Matriz4x4 inversa(const Matriz4x4& matriz);
    static Matriz4x4 transposta(const Matriz4x4& matriz);

    Point multiplicarPonto(const Point& p) const;
    Vector multiplicarVetor(const Vector& v) const;
    Matriz4x4 multiplicarMat(Matriz4x4& outra) const;

};

#endif