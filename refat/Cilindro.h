#ifndef CILINDRO_H
#define CILINDRO_H

#include "Point.h"
#include "Vector.h"

class Cilindro {
public:
    Point cb;
    float raio;
    float altura;
    Vector dc;

    Cilindro(const Point& cb_c, const float raio_c, const float altura_c, const Vector& dc_c);

    float CalcularIntersecao(const Point& origem, const Vector& dir) const;
    Vector CalcularNormal(const Point& P) const;
};

// Funções relacionadas ao Cilindro
Vector calcula_n_cilindro(Point P, Cilindro cilindro);
float calcula_t_Cilindro(Cilindro cilindro, Vector dr, Point Po);

#endif // CILINDRO_H
