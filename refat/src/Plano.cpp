#include "../include/Plano.h"
#include "../include/Utils.h"
#include "../include/Cena.h"
#include "../include/Color.h"
#include "../include/Textura.hpp"
#include "../include/Cone.h"
#include "../include/Vector.h"

#include <cmath>
#include <limits>

Plano::Plano(const Point& p_pi_p, const Vector& n_v, bool tem_textura_p)
    : p_pi(p_pi_p), n(n_v), tem_textura(tem_textura_p), material() {
    material.usarTextura = tem_textura_p;
}

Plano::Plano(const Point& p_pi_p, const Vector& n_v, const Material& material_p)
    : p_pi(p_pi_p), n(n_v), tem_textura(material_p.usarTextura), material(material_p) {
}

float Plano::CalcularIntersecao(const Point& origem, const Vector& dir) const {
    Vector w = subtrai_pontos(origem, p_pi);
    return calcula_t_plano(w, n, dir);
}

Color Plano::CalcularCor(const Cena& cena, const Vector& dir) const {
    bool naSombraEsf = false;
    bool naSombraCil = false;
    bool naSombraCone = false;

    float ti = CalcularIntersecao(cena.observador, dir);
    Point Pi = calcula_eq_ray(cena.observador, ti, dir);

    Vector l = calcula_l(cena.luz.pos, Pi);
    Point Pi_mod(Pi.x + l.i * 1e-4f, Pi.y + l.j * 1e-4f, Pi.z + l.k * 1e-4f);

    Vetor kaV = getKa();
    Vetor kdV = getKd();
    Vetor keV = getKe();
    Color Ka(kaV.i, kaV.j, kaV.k);
    Color KdBase(kdV.i, kdV.j, kdV.k);
    Color Ke(keV.i, keV.j, keV.k);

    Color Kd_textured = KdBase;
    if (material.usarTextura && material.textura != nullptr) {
        Vector nLocal = n;
        Vector arbitrary = (fabs(nLocal.i) < 0.9f) ? Vector(1.0f, 0.0f, 0.0f) : Vector(0.0f, 1.0f, 0.0f);
        Vector u_axis = cross(arbitrary, nLocal);
        float nu = calcula_norma(u_axis);
        if (nu == 0.0f) nu = 1.0f;
        u_axis = calcula_esc_por_vetor(1.0f / nu, u_axis);
        Vector v_axis = cross(nLocal, u_axis);

        Vector vecPi = subtrai_pontos(Pi, p_pi);
        float u = calcula_prod_esc(vecPi, u_axis) * material.texturaScale;
        float v = calcula_prod_esc(vecPi, v_axis) * material.texturaScale;
        u = u - floor(u);
        if (u < 0) u += 1.0f;
        v = v - floor(v);
        if (v < 0) v += 1.0f;

        size_t tx = static_cast<size_t>(u * material.textura->get_largura_pixels()) % material.textura->get_largura_pixels();
        size_t ty = static_cast<size_t>(v * material.textura->get_altura_pixels()) % material.textura->get_altura_pixels();

        rgb px = material.textura->get_cor_pixel(ty, tx);
        Color texCol(px[0] / 255.0f, px[1] / 255.0f, px[2] / 255.0f);
        Kd_textured = Color(texCol.r, texCol.g, texCol.b);
    }

    float dist_Pi_Luz = calcula_norma(subtrai_pontos(cena.luz.pos, Pi_mod));

    // Sombra da esfera
    Point centroEsf = cena.centroEsfera;
    Vector w_sombra = subtrai_pontos(Pi_mod, centroEsf);
    float a_delta = calcula_prod_esc(l, l);
    float b_delta = 2.0f * calcula_prod_esc(l, w_sombra);
    float c_delta = calcula_prod_esc(w_sombra, w_sombra) - cena.raioEsfera * cena.raioEsfera;
    float delta = b_delta * b_delta - 4.0f * a_delta * c_delta;
    float s1 = INFINITY;
    float s2 = INFINITY;
    if (delta > 0.f) {
        s1 = (-b_delta - sqrt(delta)) / (2.0f * a_delta);
        s2 = (-b_delta + sqrt(delta)) / (2.0f * a_delta);
    }
    if ((s1 > 1e-4f && s1 < dist_Pi_Luz) || (s2 > 1e-4f && s2 < dist_Pi_Luz)) naSombraEsf = true;

    // Sombra do cone
    if (!naSombraEsf && cena.cone != nullptr) {
        float t_cone = cena.cone->CalcularIntersecao(Pi_mod, l);
        if (t_cone > 1e-4f && t_cone < dist_Pi_Luz) naSombraCone = true;
    }

    // Se tem textura, aplica também no termo ambiente.
    Color Ka_textured = material.usarTextura ? Kd_textured : Ka;
    Color Ia(cena.luzAmbiente.r * Ka_textured.r, cena.luzAmbiente.g * Ka_textured.g, cena.luzAmbiente.b * Ka_textured.b);
    if (naSombraEsf || naSombraCil || naSombraCone) {
        Color I(lidarExcecao(Ia.r), lidarExcecao(Ia.g), lidarExcecao(Ia.b));
        int R = static_cast<int>(roundf(I.r * 255.0f));
        int G = static_cast<int>(roundf(I.g * 255.0f));
        int B = static_cast<int>(roundf(I.b * 255.0f));
        return Color(R, G, B);
    }

    Vector v = Vector(-dir.i, -dir.j, -dir.k);
    Vector r1 = calcula_esc_por_vetor(2.0f, calcula_esc_por_vetor(calcula_prod_esc(n, l), n));
    Vector r(r1.i - l.i, r1.j - l.j, r1.k - l.k);

    float fd = lidarExcecao(calcula_prod_esc(n, l));
    float cosAlpha = lidarExcecao(calcula_prod_esc(r, v));
    float fe = pow(cosAlpha, static_cast<float>(getShininess()));

    Color Id(cena.luz.intensidade.r * Kd_textured.r * fd,
        cena.luz.intensidade.g * Kd_textured.g * fd,
        cena.luz.intensidade.b * Kd_textured.b * fd);
    Color Ie(cena.luz.intensidade.r * Ke.r * fe,
        cena.luz.intensidade.g * Ke.g * fe,
        cena.luz.intensidade.b * Ke.b * fe);
    Color I(lidarExcecao(Id.r + Ie.r + Ia.r), lidarExcecao(Id.g + Ie.g + Ia.g), lidarExcecao(Id.b + Ie.b + Ia.b));

    int R = static_cast<int>(roundf(I.r * 255.0f));
    int G = static_cast<int>(roundf(I.g * 255.0f));
    int B = static_cast<int>(roundf(I.b * 255.0f));
    return Color(R, G, B);
}

Color Plano::CalcularCor(const Cena& cena, const Vector& dir, const Color& /*K_e*/, const Color& /*K_d*/) const {
    return CalcularCor(cena, dir);
}

bool Plano::verificarIntersecao(Vetor p0, Vetor dr) {
    // Converte Vetor -> Point (origem) e usa a interseção OO existente.
    Point origem(p0.i, p0.j, p0.k);
    Vector dir(dr.i, dr.j, dr.k);

    float t = CalcularIntersecao(origem, dir);
    if (!(t > 1e-4f) || !std::isfinite(t)) {
        setTemIntersecao(false);
        setDistancia(std::numeric_limits<double>::infinity());
        return false;
    }

    Point Pi = calcula_eq_ray(origem, t, dir);
    setTemIntersecao(true);
    setDistancia(static_cast<double>(t));
    setPontoIntersecao(Vetor(Pi.x, Pi.y, Pi.z));
    return true;
}

Vetor Plano::calcularNormal(Vetor /*posicao*/) {
    // Normal constante do plano.
    float nn = std::sqrt(n.i * n.i + n.j * n.j + n.k * n.k);
    if (nn <= 1e-8f) return Vetor(0.0f, 0.0f, 0.0f);
    return Vetor(n.i / nn, n.j / nn, n.k / nn);
}

float calcula_t_plano(Vector w, Vector n, Vector dr) {
    float result = calcula_prod_esc(calcula_esc_por_vetor(-1, w), n) / calcula_prod_esc(dr, n);
    return result;
}
