#ifndef PLANO_H
#define PLANO_H

#include "Point.h"
#include "Vector.h"
#include "Material.h"
#include "Objeto.h"

class Color;
struct Cena;

class Plano : public Objeto {
public:
    Point p_pi;
    Vector n;
    bool tem_textura;
    Material material;

    Plano(const Point& p_pi_p, const Vector& n_v, bool tem_textura_p = false);
    Plano(const Point& p_pi_p, const Vector& n_v, const Material& material_p);

    // Interseção do raio (origem + t*dir) com o plano.
    // Retorna t (pode ser negativo se o plano estiver "atrás" do raio).
    float CalcularIntersecao(const Point& origem, const Vector& dir) const;

    // Calcula a cor no ponto de interseção do raio com este plano.
    Color CalcularCor(const Cena& cena, const Vector& dir) const;

    // Compatibilidade: versão antiga que recebe Ks (agora ignorados)
    Color CalcularCor(const Cena& cena, const Vector& dir, const Color& K_e, const Color& K_d) const;

    bool verificarIntersecao(Vetor p0, Vetor dr) override;
    Vetor calcularNormal(Vetor posicao) override;
    void aplicarEscalaNoPivoObjeto(const Vetor& escala, const Point& pivo) override;
    void aplicarEscalaUniforme(float s) override;

    Point aplicarTranslacao(const Point& ponto, const Vetor& translacao) const {
        return Objeto::aplicarTranslacao(ponto, translacao);
    }

    Point aplicarEscala(const Point& ponto, const Vetor& escala) const {
        return Objeto::aplicarEscala(ponto, escala);
    }
};

// Função relacionada ao Plano
float calcula_t_plano(Vector w, Vector n, Vector dr);

#endif // PLANO_H
