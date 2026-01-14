#include "Plano.h"
#include "Utils.h"

Plano::Plano(const Point& p_pi_p, const Vector& n_v, bool tem_textura_p)
    : p_pi(p_pi_p), n(n_v), tem_textura(tem_textura_p) {
}

float Plano::CalcularIntersecao(const Point& origem, const Vector& dir) const {
    Vector w = subtrai_pontos(origem, p_pi);
    return calcula_t_plano(w, n, dir);
}

float calcula_t_plano(Vector w, Vector n, Vector dr) {
    float result = calcula_prod_esc(calcula_esc_por_vetor(-1, w), n) / calcula_prod_esc(dr, n);
    return result;
}
