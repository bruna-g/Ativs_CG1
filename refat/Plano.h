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
};

// Função relacionada ao Plano
float calcula_t_plano(Vector w, Vector n, Vector dr);

#endif // PLANO_H
