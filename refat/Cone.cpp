#include "Cone.h"
#include "Utils.h"
#include <cmath>
#include <algorithm>

Cone::Cone(const Point& cb_c, const Point& v_c, const float raio_c)
    : cb(cb_c), v(v_c), raio(raio_c), dc(0.0f, 1.0f, 0.0f), altura(0.0f), teta(0.0f) {
    Vector cbv = subtrai_pontos(v, cb);
    altura = calcula_norma(cbv);
    if (altura > 0.0f) {
        dc = Vector(cbv.i / altura, cbv.j / altura, cbv.k / altura);
    }
    teta = atan(raio / (altura > 0.0f ? altura : 1.0f));
}

float Cone::CalcularIntersecao(const Point& origem, const Vector& dir) const {
    // Cone finito com base em cb e vértice em v.
    // Eixo dc aponta de cb -> v; altura = |v-cb|; semi-ângulo teta.
    float cosT = cos(teta);
    float cos2 = cosT * cosT;

    // co = origem - ápice
    Vector co = subtrai_pontos(origem, v);

    float dv = calcula_prod_esc(dir, dc);
    float cov = calcula_prod_esc(co, dc);

    float a = dv * dv - cos2 * calcula_prod_esc(dir, dir);
    float b = dv * cov - cos2 * calcula_prod_esc(dir, co); // coeficiente do termo 2*b*t
    float c = cov * cov - cos2 * calcula_prod_esc(co, co);

    float tEscolhido = INFINITY;
    bool temIntersecao = false;

    // Interseção com o disco da base (plano passando em cb, normal dc)
    float dn = calcula_prod_esc(dc, dir);
    if (fabs(dn) > 1e-8f) {
        Vector w_base = subtrai_pontos(origem, cb);
        float t_base = -(calcula_prod_esc(w_base, dc) / dn);

        if (t_base > 1e-4f) {
            Point p_base = calcula_eq_ray(origem, t_base, dir);
            Vector dist_centro = subtrai_pontos(p_base, cb);
            float dist2 = calcula_prod_esc(dist_centro, dist_centro);
            if (dist2 <= raio * raio) {
                tEscolhido = std::min(tEscolhido, t_base);
                temIntersecao = true;
            }
        }
    }

    // Interseção com a lateral (cone infinito + recorte pela altura)
    float delta = b * b - a * c;
    if (delta < 0.0f) {
        return temIntersecao ? tEscolhido : -1.0f;
    }

    // Caso degenerado: a ~ 0 (equação linear)
    if (fabs(a) < 1e-8f && fabs(b) > 1e-8f) {
        float t_lin = -c / (2.0f * b);
        if (t_lin > 1e-4f) {
            Point Pi = calcula_eq_ray(origem, t_lin, dir);
            Vector v_pi = subtrai_pontos(v, Pi);
            float h = calcula_prod_esc(v_pi, dc);
            if (h >= 0.0f && h <= altura) {
                tEscolhido = std::min(tEscolhido, t_lin);
                temIntersecao = true;
            }
        }
        return temIntersecao ? tEscolhido : -1.0f;
    }

    float sqrtD = sqrt(delta);
    float t1 = (-b - sqrtD) / (a);
    float t2 = (-b + sqrtD) / (a);

    auto validarT = [&](float t) -> bool {
        if (t <= 1e-4f) return false;
        Point Pi = calcula_eq_ray(origem, t, dir);
        Vector v_pi = subtrai_pontos(v, Pi);
        float h = calcula_prod_esc(v_pi, dc);
        return (h >= 0.0f && h <= altura);
        };

    bool t1_ok = validarT(t1);
    bool t2_ok = validarT(t2);

    if (t1_ok) {
        tEscolhido = std::min(tEscolhido, t1);
        temIntersecao = true;
    }
    if (t2_ok) {
        tEscolhido = std::min(tEscolhido, t2);
        temIntersecao = true;
    }

    return temIntersecao ? tEscolhido : -1.0f;
}

float Cone::CalcularIntersecao(const Point& origem, const Point& pontoJanela) const {
    Vector dir = calcula_dr(origem, pontoJanela);
    return CalcularIntersecao(origem, dir);
}

// Wrapper legado (mantido para compatibilidade com código antigo)
extern Point Po;
extern Cone cone;
float calcula_t_cone(Point& Pj) {
    Vector dir = calcula_dr(Po, Pj);
    return cone.CalcularIntersecao(Po, dir);
}
