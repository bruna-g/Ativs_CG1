#ifndef OBJETO_H
#define OBJETO_H

#include "Point.h"
#include "Vetor.h"

#include <cmath>

class Objeto {
private:
    Vetor Ke;
    Vetor Kd;
    Vetor Ka;
    double shininess;

    double distancia;
    Vetor pontoIntersecao;

    bool temIntersecao;

public:
    virtual ~Objeto() = default;

    Objeto()
        : Ke(0.0f, 0.0f, 0.0f),
        Kd(0.0f, 0.0f, 0.0f),
        Ka(0.0f, 0.0f, 0.0f),
        shininess(0.0),
        distancia(0.0),
        pontoIntersecao(0.0f, 0.0f, 0.0f),
        temIntersecao(false) {
    }

    virtual bool verificarIntersecao(Vetor p0, Vetor dr) = 0;
    virtual Vetor calcularNormal(Vetor posicao) = 0;

    // Padronizado: aplica escala (possivelmente não-uniforme) no PRÓPRIO objeto
    // em torno de um pivô (pivo), atualizando seus parâmetros internos.
    virtual void aplicarEscalaNoPivoObjeto(const Vetor& escala, const Point& pivo) = 0;

    virtual void aplicarEscalaUniforme(float s) = 0;

    Point aplicarTranslacao(const Point& ponto, const Vetor& translacao) const {
        // Matriz homogênea 4x4 de translação:
        // [1 0 0 tx]
        // [0 1 0 ty]
        // [0 0 1 tz]
        // [0 0 0  1]
        const float tx = translacao.i;
        const float ty = translacao.j;
        const float tz = translacao.k;

        const float x = ponto.x + tx * ponto.p;
        const float y = ponto.y + ty * ponto.p;
        const float z = ponto.z + tz * ponto.p;
        const float w = ponto.p;

        if (std::fabs(w) > 1e-6f && std::fabs(w - 1.0f) > 1e-6f) {
            return Point(x / w, y / w, z / w);
        }
        return Point(x, y, z);
    }

    Point aplicarEscala(const Point& ponto, const Vetor& escala) const {
        // Matriz homogênea 4x4 de translação:
        // [Sx 0 0 0]
        // [0 Sy 0 0]
        // [0 0 Sz 0]
        // [0 0 0  1]
        const float sx = escala.i;
        const float sy = escala.j;
        const float sz = escala.k;

        const float x = ponto.x * sx;
        const float y = ponto.y * sy;
        const float z = ponto.z * sz;
        const float w = 1 * ponto.p;

        if (std::fabs(w) > 1e-6f && std::fabs(w - 1.0f) > 1e-6f) {
            return Point(x / w, y / w, z / w);
        }
        return Point(x, y, z);
    }

    // Escala NÃO-uniforme em torno de um pivô (pivo):
    // P' = T(pivo) * S(escala) * T(-pivo) * P
    Point aplicarEscalaNoPivo(const Point& ponto, const Vetor& escala, const Point& pivo) const {
        const Vetor paraOrigem(-pivo.x, -pivo.y, -pivo.z);
        const Vetor deVolta(pivo.x, pivo.y, pivo.z);

        Point p = aplicarTranslacao(ponto, paraOrigem);
        p = aplicarEscala(p, escala);
        p = aplicarTranslacao(p, deVolta);
        return p;
    }

    void setKe(Vetor Ke_p) { Ke = Ke_p; }
    Vetor getKe() const { return Ke; }

    void setKd(Vetor Kd_p) { Kd = Kd_p; }
    Vetor getKd() const { return Kd; }

    void setKa(Vetor Ka_p) { Ka = Ka_p; }
    Vetor getKa() const { return Ka; }

    void setShininess(double shininess_p) { shininess = shininess_p; }
    double getShininess() const { return shininess; }

    void setDistancia(double distancia_p) { distancia = distancia_p; }
    double getDistancia() const { return distancia; }

    void setPontoIntersecao(Vetor pontoIntersecao_p) { pontoIntersecao = pontoIntersecao_p; }
    Vetor getPontoIntersecao() const { return pontoIntersecao; }

    void setTemIntersecao(bool temIntersecao_p) { temIntersecao = temIntersecao_p; }
    bool getTemIntersecao() const { return temIntersecao; }
};

#endif // OBJETO_H
