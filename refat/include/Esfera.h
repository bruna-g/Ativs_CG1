#ifndef ESFERA_H
#define ESFERA_H

#include "Objeto.h"
#include "Point.h"
#include "Vector.h"
#include "Material.h"

class Color;
struct Cena;

// Esfera com suporte à interface OO (Objeto) e também aos métodos antigos
// usados por alguns arquivos (CalcularCor/CalcularNormal/CalcularIntersecao).
class Esfera : public Objeto {
public:
    Point centro;
    float raio;
    Material material;

    Esfera(const Point& centro_e, const float raio_e);
    Esfera(const Point& centro_e, const float raio_e, const Material& material_p);

    void aplicarEscalaUniforme(float s);

    // Interseção (origem + t*dir) com a esfera.
    // Observação: Vetor é alias de Vector (Vetor.h), então uma única assinatura basta.
    float CalcularIntersecao(const Point& origem, const Vector& dir) const;

    // Normal em um ponto da superfície.
    Vector CalcularNormal(const Point& P) const;

    // Shading (pipeline antigo)
    Color CalcularCor(const Cena& cena, float t, const Vector& dir) const;
    Color CalcularCor(const Cena& cena, float t, const Vector& dir,
        const Color& K_a, const Color& K_d, const Color& K_e) const;

    // Interface Objeto
    bool verificarIntersecao(Vetor p0, Vetor dr) override;
    Vetor calcularNormal(Vetor posicao) override;
    void aplicarEscalaNoPivoObjeto(const Vetor& escala, const Point& pivo) override;

    Point aplicarTranslacao(const Point& ponto, const Vetor& translacao) const {
        return Objeto::aplicarTranslacao(ponto, translacao);
    }

    Point aplicarEscala(const Point& ponto, const Vetor& escala) const {
        return Objeto::aplicarEscala(ponto, escala);
    }
};

// Funções relacionadas à Esfera
Vector calcula_n_esfera(Point P, Esfera esfera);
float calcula_t_esfera(Esfera esfera, Vector dr, Point Po);

#endif // ESFERA_H
