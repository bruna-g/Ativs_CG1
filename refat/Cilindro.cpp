#include "Cilindro.h"
#include "Utils.h"
#include <cmath>
#include <algorithm>

Cilindro::Cilindro(const Point& cb_c, const float raio_c, const float altura_c, const Vector& dc_c)
    : cb(cb_c), raio(raio_c), altura(altura_c), dc(dc_c) {
}

Vector Cilindro::CalcularNormal(const Point& P) const {
    Point cbCopy = cb;
    Point pCopy = P;
    Vector P_B = subtrai_pontos(pCopy, cbCopy);
    float P_B_u = calcula_prod_esc(P_B, dc);
    Vector P_B_u_u = calcula_esc_por_vetor(P_B_u, dc);
    Vector n = subtrai_vetores(P_B, P_B_u_u);
    float norma = calcula_norma(n);
    return Vector(n.i / norma, n.j / norma, n.k / norma);
}

float Cilindro::CalcularIntersecao(const Point& origem, const Vector& dir) const {
    Point cbCopy = cb;
    Point origemCopy = origem;
    Vector Po_B = subtrai_pontos(origemCopy, cbCopy);
    Vector Po_B_u_u = calcula_esc_por_vetor(calcula_prod_esc(Po_B, dc), dc);
    Vector v = subtrai_vetores(Po_B, Po_B_u_u);

    Vector d_u_u = calcula_esc_por_vetor(calcula_prod_esc(dir, dc), dc);
    Vector w = subtrai_vetores(dir, d_u_u);

    float a_delta = calcula_prod_esc(w, w);
    float b_delta = calcula_prod_esc(v, w);
    float c_delta = calcula_prod_esc(v, v) - raio * raio;

    float delta = b_delta * b_delta - a_delta * c_delta;
    if (delta < 0) return INFINITY;

    float t1 = (-b_delta + sqrt(delta)) / (a_delta);
    float t2 = (-b_delta - sqrt(delta)) / (a_delta);

    Point P1 = calcula_eq_ray(origem, t1, dir);
    Point P2 = calcula_eq_ray(origem, t2, dir);

    Point cbTmp1 = cb;
    Point cbTmp2 = cb;
    float P1_B_u = calcula_prod_esc(subtrai_pontos(P1, cbTmp1), dc);
    float P2_B_u = calcula_prod_esc(subtrai_pontos(P2, cbTmp2), dc);

    float t = -1.0f;
    if ((P1_B_u >= 0 && P1_B_u <= altura && t1 > 0) &&
        (P2_B_u >= 0 && P2_B_u <= altura && t2 > 0))
        t = std::min(t1, t2);
    else if (P1_B_u >= 0 && P1_B_u <= altura && t1 > 0)
        t = t1;
    else if (P2_B_u >= 0 && P2_B_u <= altura && t2 > 0)
        t = t2;

    return t;
}

Vector calcula_n_cilindro(Point P, Cilindro cilindro) {
    return cilindro.CalcularNormal(P);
}

float calcula_t_Cilindro(Cilindro cilindro, Vector dr, Point Po) {
    return cilindro.CalcularIntersecao(Po, dr);
}
