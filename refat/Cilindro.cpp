#include "Cilindro.h"
#include "Utils.h"
#include <cmath>
#include <algorithm>

Cilindro::Cilindro(const Point& cb_c, const float raio_c, const float altura_c, const Vector& dc_c) 
    : cb(cb_c), raio(raio_c), altura(altura_c), dc(dc_c) {}

Vector calcula_n_cilindro(Point P, Cilindro cilindro) {
    Vector P_B = subtrai_pontos(P, cilindro.cb);
    float P_B_u = calcula_prod_esc(P_B, cilindro.dc);
    Vector P_B_u_u = calcula_esc_por_vetor(P_B_u, cilindro.dc);
    Vector n = subtrai_vetores(P_B, P_B_u_u);
    float norma = calcula_norma(n);
    return Vector(n.i / norma, n.j / norma, n.k / norma);
}

float calcula_t_Cilindro(Cilindro cilindro, Vector dr, Point Po) {
    Vector Po_B = subtrai_pontos(Po, cilindro.cb);
    Vector Po_B_u_u = calcula_esc_por_vetor(calcula_prod_esc(Po_B, cilindro.dc), cilindro.dc);
    Vector v = subtrai_vetores(Po_B, Po_B_u_u);

    Vector d_u_u = calcula_esc_por_vetor(calcula_prod_esc(dr, cilindro.dc), cilindro.dc);
    Vector w = subtrai_vetores(dr, d_u_u);

    float a_delta = calcula_prod_esc(w, w);
    float b_delta = calcula_prod_esc(v, w);
    float c_delta = calcula_prod_esc(v, v) - cilindro.raio * cilindro.raio;

    float delta = b_delta * b_delta - a_delta * c_delta;
    if (delta >= 0) {
        float t1 = (-b_delta + sqrt(delta)) / (a_delta);
        float t2 = (-b_delta - sqrt(delta)) / (a_delta);

        Point P1 = calcula_eq_ray(Po, t1, dr);
        Point P2 = calcula_eq_ray(Po, t2, dr);

        float P1_B_u = calcula_prod_esc(subtrai_pontos(P1, cilindro.cb), cilindro.dc);
        float P2_B_u = calcula_prod_esc(subtrai_pontos(P2, cilindro.cb), cilindro.dc);

        float t = -1.0f;
        if ((P1_B_u >= 0 && P1_B_u <= cilindro.altura && t1 > 0) &&
            (P2_B_u >= 0 && P2_B_u <= cilindro.altura && t2 > 0))
            t = std::min(t1, t2);
        else if (P1_B_u >= 0 && P1_B_u <= cilindro.altura && t1 > 0)
            t = t1;
        else if (P2_B_u >= 0 && P2_B_u <= cilindro.altura && t2 > 0)
            t = t2;
        return t;
    }
    else {
        return INFINITY;
    }
}
