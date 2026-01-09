#ifndef CUBO_H
#define CUBO_H

#include "Point.h"
#include "Vector.h"

class Cubo {
public:
    Point centro_base;
    float aresta;
    
    Cubo(const Point& centro_c, const float aresta_c);
};

struct IntersecaoCubo {
    bool intercepta;
    float t;
    Vector normal;
    Point ponto;
};

// Função relacionada ao Cubo
IntersecaoCubo calcula_intersecao_cubo_completa(const Cubo& cubo, const Vector& dr, Point& Po);

#endif // CUBO_H
