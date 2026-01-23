#ifndef CILINDRO_H
#define CILINDRO_H

#include "Point.h"
#include "Vector.h"
#include "Material.h"
#include "Objeto.h"

class Color;
struct Cena;

class Cilindro : public Objeto {
public:
    Point cb;
    float raio;
    float altura;
    Vector dc;
    Material material;

    Cilindro(const Point& cb_c, const float raio_c, const float altura_c, const Vector& dc_c);
    Cilindro(const Point& cb_c, const float raio_c, const float altura_c, const Vector& dc_c, const Material& material_p);

    void aplicarEscalaUniforme(float s) override;

    float CalcularIntersecao(const Point& origem, const Vector& dir) const;
    Vector CalcularNormal(const Point& P) const;

    Color CalcularCor(const Cena& cena, float t, const Vector& dir) const;
    Color CalcularCor(const Cena& cena, float t, const Vector& dir,
        const Color& K_a, const Color& K_d, const Color& K_e) const;

    bool verificarIntersecao(Vetor p0, Vetor dr) override;
    Vetor calcularNormal(Vetor posicao) override;
    void aplicarEscalaNoPivoObjeto(const Vetor& escala) override;

    Point aplicarTranslacao(const Point& ponto, const Vetor& translacao) const {
        return Objeto::aplicarTranslacao(ponto, translacao);
    }

    Point aplicarEscala(const Point& ponto, const Vetor& escala) const {
        return Objeto::aplicarEscala(ponto, escala);
    }

    void rotacionarX(float angulo, Cilindro cilindro);
    void rotacionarY(float angulo);
    void rotacionarZ(float angulo);

};

// Funções relacionadas ao Cilindro
Vector calcula_n_cilindro(Point P, Cilindro cilindro);
float calcula_t_Cilindro(Cilindro cilindro, Vector dr, Point Po);
void rotacionarArbitrario(const Vector& eixo, float angulo);

Point calcularCentro(Cilindro cilindro);

#endif // CILINDRO_H
