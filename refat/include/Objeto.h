#ifndef OBJETO_H
#define OBJETO_H

#include "Vetor.h"

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
