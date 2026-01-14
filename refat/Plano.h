#ifndef PLANO_H
#define PLANO_H

#include "Point.h"
#include "Vector.h"

class Plano {
public:
    Point p_pi;
    Vector n;
    bool tem_textura;

    Plano(const Point& p_pi_p, const Vector& n_v, bool tem_textura_p = false);

    // Interseção do raio (origem + t*dir) com o plano.
    // Retorna t (pode ser negativo se o plano estiver "atrás" do raio).
    float CalcularIntersecao(const Point& origem, const Vector& dir) const;
};

// Função relacionada ao Plano
float calcula_t_plano(Vector w, Vector n, Vector dr);

#endif // PLANO_H
