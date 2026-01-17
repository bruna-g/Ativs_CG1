#ifndef ESFERA_H
#define ESFERA_H

#include "Point.h"
#include "Vector.h"
#include "Material.h"

class Color;
struct Cena;

class Esfera {
public:
    Point centro;
    float raio;
    Material material;

    Esfera(const Point& centro_e, const float raio_e);
    Esfera(const Point& centro_e, const float raio_e, const Material& material_p);

    float CalcularIntersecao(const Point& origem, const Vector& dir) const;
    Vector CalcularNormal(const Point& P) const;

    Color CalcularCor(const Cena& cena, float t, const Vector& dir) const;
    Color CalcularCor(const Cena& cena, float t, const Vector& dir,
        const Color& K_a, const Color& K_d, const Color& K_e) const;
};

// Funções relacionadas à Esfera
Vector calcula_n_esfera(Point P, Esfera esfera);
float calcula_t_esfera(Esfera esfera, Vector dr, Point Po);

#endif // ESFERA_H
