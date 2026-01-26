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

    void recalcularDerivados();
    void aplicarEscalaUniforme(float s) override;

    float CalcularIntersecao(const Point& origem, const Vector& dir) const;
    float CalcularIntersecao(const Point& origem, const Point& pontoJanela) const;

    Color CalcularCor(const Cena& cena, float t, const Vector& dir) const;
    Color CalcularCor(const Cena& cena, float t, const Vector& dir, const Color& K) const;

    bool verificarIntersecao(Vetor p0, Vetor dr) override;
    Vetor calcularNormal(Vetor posicao) override;
    void aplicarEscalaNoPivoObjeto(const Vetor& escala, const Point& pivo) override;

    Point aplicarTranslacao(const Point& ponto, const Vetor& translacao) const {
        return Objeto::aplicarTranslacao(ponto, translacao);
    }

    Point aplicarEscala(const Point& ponto, const Vetor& escala) const {
        return Objeto::aplicarEscala(ponto, escala);
    }

    Vector CalcularNormal(const Point& Pi, const Vector& dir) const;

    void rotacionarX(float angulo);
    void rotacionarY(float angulo);
    void rotacionarZ(float angulo);
    void rotacaoArbitraria(Vector eixo, float angulo);
};

// Função relacionada ao Cone
float calcula_t_cone(Point& Pj);

#endif // CONE_H
