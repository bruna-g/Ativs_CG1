#ifndef CONE_H
#define CONE_H

#include "Point.h"
#include "Vector.h"
#include "Material.h"
#include "Objeto.h"

class Color;
struct Cena;

class Cone : public Objeto {
public:
    Point cb;
    Point v;
    float raio;

    Material material;

    // Derivados (computados no construtor)
    Vector dc;     // eixo (unitário) da base para o vértice
    float altura;  // |v - cb|
    float teta;    // atan(raio/altura)

    Cone(const Point& cb_c, const Point& v_c, const float raio_c);
    Cone(const Point& cb_c, const Point& v_c, const float raio_c, const Material& material_p);

    float CalcularIntersecao(const Point& origem, const Vector& dir) const;
    float CalcularIntersecao(const Point& origem, const Point& pontoJanela) const;

    Color CalcularCor(const Cena& cena, float t, const Vector& dir) const;
    Color CalcularCor(const Cena& cena, float t, const Vector& dir, const Color& K) const;

    bool verificarIntersecao(Vetor p0, Vetor dr) override;
    Vetor calcularNormal(Vetor posicao) override;
};

// Função relacionada ao Cone
float calcula_t_cone(Point& Pj);

#endif // CONE_H
