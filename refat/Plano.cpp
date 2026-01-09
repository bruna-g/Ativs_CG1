#include "Plano.h"
#include "Utils.h"

Plano::Plano(const Point& p_pi_p, const Vector& n_v, bool tem_textura_p) 
    : p_pi(p_pi_p), n(n_v), tem_textura(tem_textura_p) {}

float calcula_t_plano(Vector w, Vector n, Vector dr) {
    float result = calcula_prod_esc(calcula_esc_por_vetor(-1, w), n) / calcula_prod_esc(dr, n);
    return result;
}
