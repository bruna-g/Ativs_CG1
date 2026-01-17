#include "../include/Esfera.h"
#include "../include/Cena.h"
#include "../include/Color.h"
#include "../include/Utils.h"
#include <cmath>
#include <algorithm>

Esfera::Esfera(const Point& centro_e, const float raio_e)
    : centro(centro_e), raio(raio_e), material() {
}

Esfera::Esfera(const Point& centro_e, const float raio_e, const Material& material_p)
    : centro(centro_e), raio(raio_e), material(material_p) {
}

Vector Esfera::CalcularNormal(const Point& P) const {
    Vector n = subtrai_pontos(P, centro);
    float norma = calcula_norma(n);
    return Vector(n.i / norma, n.j / norma, n.k / norma);
}

float Esfera::CalcularIntersecao(const Point& origem, const Vector& dir) const {
    Point origemCopy = origem;
    Point centroCopy = centro;
    Vector w = subtrai_pontos(origemCopy, centroCopy);

    float a_delta = calcula_prod_esc(dir, dir);
    float b_delta = 2.0f * calcula_prod_esc(dir, w);
    float c_delta = calcula_prod_esc(w, w) - raio * raio;

    float delta = b_delta * b_delta - 4.0f * a_delta * c_delta;
    if (delta < 0) return -1.0f;

    float t1 = (-b_delta + sqrt(delta)) / (2.0f * a_delta);
    float t2 = (-b_delta - sqrt(delta)) / (2.0f * a_delta);

    float t = -1.0f;
    if (t1 > 0.0f && t2 > 0.0f)
        t = std::min(t1, t2);
    else if (t1 > 0.0f)
        t = t1;
    else if (t2 > 0.0f)
        t = t2;

    return t;
}

Color Esfera::CalcularCor(const Cena& cena, float t, const Vector& dir,
    const Color& K_a, const Color& K_d, const Color& K_e) const {
    Point P = calcula_eq_ray(cena.observador, t, dir);

    bool naSombra = false;
    float dist_Pi_Luz = calcula_norma(subtrai_pontos(cena.luz.pos, P));
    if (t < dist_Pi_Luz) naSombra = true;

    Color Ia(cena.luzAmbiente.r * K_a.r, cena.luzAmbiente.g * K_a.g, cena.luzAmbiente.b * K_a.b);
    if (naSombra) {
        Color I(lidarExcecao(Ia.r), lidarExcecao(Ia.g), lidarExcecao(Ia.b));
        int R = static_cast<int>(roundf(I.r * 255.0f));
        int G = static_cast<int>(roundf(I.g * 255.0f));
        int B = static_cast<int>(roundf(I.b * 255.0f));
        return Color(R, G, B);
    }

    Vector n = CalcularNormal(P);
    Vector l = calcula_l(cena.luz.pos, P);
    Vector v(-dir.i, -dir.j, -dir.k);
    Vector r1 = calcula_esc_por_vetor(2.0f, calcula_esc_por_vetor(calcula_prod_esc(n, l), n));
    Vector r(r1.i - l.i, r1.j - l.j, r1.k - l.k);

    float fd = lidarExcecao(calcula_prod_esc(n, l));
    float cosAlpha = lidarExcecao(calcula_prod_esc(r, v));
    float fe = pow(cosAlpha, cena.expoenteEspecular);

    Color Id(cena.luz.intensidade.r * K_d.r * fd,
        cena.luz.intensidade.g * K_d.g * fd,
        cena.luz.intensidade.b * K_d.b * fd);
    Color Ie(cena.luz.intensidade.r * K_e.r * fe,
        cena.luz.intensidade.g * K_e.g * fe,
        cena.luz.intensidade.b * K_e.b * fe);
    Color I(lidarExcecao(Id.r + Ie.r + Ia.r), lidarExcecao(Id.g + Ie.g + Ia.g), lidarExcecao(Id.b + Ie.b + Ia.b));

    int R = static_cast<int>(roundf(I.r * 255.0f));
    int G = static_cast<int>(roundf(I.g * 255.0f));
    int B = static_cast<int>(roundf(I.b * 255.0f));
    return Color(R, G, B);
}

Color Esfera::CalcularCor(const Cena& cena, float t, const Vector& dir) const {
    Point P = calcula_eq_ray(cena.observador, t, dir);

    bool naSombra = false;
    float dist_Pi_Luz = calcula_norma(subtrai_pontos(cena.luz.pos, P));
    if (t < dist_Pi_Luz) naSombra = true;

    Color Ia(cena.luzAmbiente.r * material.Ka.r, cena.luzAmbiente.g * material.Ka.g, cena.luzAmbiente.b * material.Ka.b);
    if (naSombra) {
        Color I(lidarExcecao(Ia.r), lidarExcecao(Ia.g), lidarExcecao(Ia.b));
        int R = static_cast<int>(roundf(I.r * 255.0f));
        int G = static_cast<int>(roundf(I.g * 255.0f));
        int B = static_cast<int>(roundf(I.b * 255.0f));
        return Color(R, G, B);
    }

    Vector n = CalcularNormal(P);
    Vector l = calcula_l(cena.luz.pos, P);
    Vector v(-dir.i, -dir.j, -dir.k);
    Vector r1 = calcula_esc_por_vetor(2.0f, calcula_esc_por_vetor(calcula_prod_esc(n, l), n));
    Vector r(r1.i - l.i, r1.j - l.j, r1.k - l.k);

    float fd = lidarExcecao(calcula_prod_esc(n, l));
    float cosAlpha = lidarExcecao(calcula_prod_esc(r, v));
    float fe = pow(cosAlpha, material.m);

    Color Id(cena.luz.intensidade.r * material.Kd.r * fd,
        cena.luz.intensidade.g * material.Kd.g * fd,
        cena.luz.intensidade.b * material.Kd.b * fd);
    Color Ie(cena.luz.intensidade.r * material.Ke.r * fe,
        cena.luz.intensidade.g * material.Ke.g * fe,
        cena.luz.intensidade.b * material.Ke.b * fe);
    Color I(lidarExcecao(Id.r + Ie.r + Ia.r), lidarExcecao(Id.g + Ie.g + Ia.g), lidarExcecao(Id.b + Ie.b + Ia.b));

    int R = static_cast<int>(roundf(I.r * 255.0f));
    int G = static_cast<int>(roundf(I.g * 255.0f));
    int B = static_cast<int>(roundf(I.b * 255.0f));
    return Color(R, G, B);
}

Vector calcula_n_esfera(Point P, Esfera esfera) {
    return esfera.CalcularNormal(P);
}

float calcula_t_esfera(Esfera esfera, Vector dr, Point Po) {
    return esfera.CalcularIntersecao(Po, dr);
}
