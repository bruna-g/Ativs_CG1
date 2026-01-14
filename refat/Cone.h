#ifndef CONE_H
#define CONE_H

#include "Point.h"
#include "Vector.h"

class Cone {
public:
    Point cb;
    Point v;
    float raio;

    // Derivados (computados no construtor)
    Vector dc;     // eixo (unitário) da base para o vértice
    float altura;  // |v - cb|
    float teta;    // atan(raio/altura)

    Cone(const Point& cb_c, const Point& v_c, const float raio_c);

    float CalcularIntersecao(const Point& origem, const Vector& dir) const;
    float CalcularIntersecao(const Point& origem, const Point& pontoJanela) const;
};

// Função relacionada ao Cone
float calcula_t_cone(Point& Pj);

#endif // CONE_H
