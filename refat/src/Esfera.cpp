#include "../include/Esfera.h"

#include "../include/Cena.h"
#include "../include/Color.h"
#include "../include/Utils.h"

#include <algorithm>
#include <cmath>
#include <limits>

Esfera::Esfera(const Point& centro_e, const float raio_e)
    : centro(centro_e), raio(raio_e), material() {
}

Esfera::Esfera(const Point& centro_e, const float raio_e, const Material& material_p)
    : centro(centro_e), raio(raio_e), material(material_p) {
}

void Esfera::aplicarEscalaUniforme(float s) {
    if (s <= 0.0f) return;
    raio *= s;
}

void Esfera::aplicarEscalaNoPivoObjeto(const Vetor& escala) {
    centro = Objeto::aplicarEscalaNoPivo(centro, escala, Point(0.0f, 0.0f, 0.0f));

    const float sx = std::fabs(escala.i);
    const float sy = std::fabs(escala.j);
    const float sz = std::fabs(escala.k);

    const float eps = 1e-6f;
    const bool uniforme = (std::fabs(sx - sy) < eps) && (std::fabs(sy - sz) < eps);
    const float s = uniforme ? sx : std::max(sx, std::max(sy, sz));

    raio *= s;
}

Vector Esfera::CalcularNormal(const Point& P) const {
    Vector n = subtrai_pontos(P, centro);
    float norma = calcula_norma(n);
    if (norma <= 1e-8f) return Vector(0.0f, 0.0f, 0.0f);
    return Vector(n.i / norma, n.j / norma, n.k / norma);
}

float Esfera::CalcularIntersecao(const Point& origem, const Vector& dir) const {
    // Mantém a assinatura antiga (Vector) e retorna +infinity quando não houver interseção.
    // Implementação mais robusta: usa half_b para reduzir perda numérica.
    const Vector w = subtrai_pontos(origem, centro);

    const double a = static_cast<double>(calcula_prod_esc(dir, dir));
    if (std::abs(a) < 1e-12) return std::numeric_limits<float>::infinity();

    const double half_b = static_cast<double>(calcula_prod_esc(w, dir));
    const double c = static_cast<double>(calcula_prod_esc(w, w)) - static_cast<double>(raio) * static_cast<double>(raio);

    const double delta = half_b * half_b - a * c;
    if (delta < 0.0) return std::numeric_limits<float>::infinity();

    const double sqrtD = std::sqrt(delta);
    const double t1 = (-half_b - sqrtD) / a;
    const double t2 = (-half_b + sqrtD) / a;

    const double eps = 1e-4;
    double t = std::numeric_limits<double>::infinity();
    if (t1 > eps) t = t1;
    if (t2 > eps) t = std::min(t, t2);
    if (!std::isfinite(t)) return std::numeric_limits<float>::infinity();

    return static_cast<float>(t);
}

Color Esfera::CalcularCor(const Cena& cena, float t, const Vector& dir,
    const Color& K_a, const Color& K_d, const Color& K_e) const {
    Point P = calcula_eq_ray(cena.observador, t, dir);

    Vector l = calcula_l(cena.luz.pos, P);
    Point P_mod(P.x + l.i * 1e-4f, P.y + l.j * 1e-4f, P.z + l.k * 1e-4f);
    float dist_Pi_Luz = calcula_norma(subtrai_pontos(cena.luz.pos, P_mod));
    bool naSombra = cena.estaNaSombra(P_mod, l, dist_Pi_Luz, const_cast<Esfera*>(this));

    Color Ia(cena.luzAmbiente.r * K_a.r, cena.luzAmbiente.g * K_a.g, cena.luzAmbiente.b * K_a.b);
    if (naSombra) {
        Color I(lidarExcecao(Ia.r), lidarExcecao(Ia.g), lidarExcecao(Ia.b));
        int R = static_cast<int>(roundf(I.r * 255.0f));
        int G = static_cast<int>(roundf(I.g * 255.0f));
        int B = static_cast<int>(roundf(I.b * 255.0f));
        return Color(R, G, B);
    }

    Vector n = CalcularNormal(P);
    Vector v(-dir.i, -dir.j, -dir.k);
    Vector r1 = calcula_esc_por_vetor(2.0f, calcula_esc_por_vetor(calcula_prod_esc(n, l), n));
    Vector r(r1.i - l.i, r1.j - l.j, r1.k - l.k);

    float fd = lidarExcecao(calcula_prod_esc(n, l));
    fd = std::clamp(fd, 0.0f, 1.0f);

    float cosAlpha = lidarExcecao(calcula_prod_esc(r, v));
    cosAlpha = std::clamp(cosAlpha, 0.0f, 1.0f);

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

    Vector l = calcula_l(cena.luz.pos, P);
    Point P_mod(P.x + l.i * 1e-4f, P.y + l.j * 1e-4f, P.z + l.k * 1e-4f);
    float dist_Pi_Luz = calcula_norma(subtrai_pontos(cena.luz.pos, P_mod));
    bool naSombra = cena.estaNaSombra(P_mod, l, dist_Pi_Luz, const_cast<Esfera*>(this));

    Vetor kaV = getKa();
    Vetor kdV = getKd();
    Vetor keV = getKe();
    Color Ka(kaV.i, kaV.j, kaV.k);
    Color Kd(kdV.i, kdV.j, kdV.k);
    Color Ke(keV.i, keV.j, keV.k);

    Color Ia(cena.luzAmbiente.r * Ka.r, cena.luzAmbiente.g * Ka.g, cena.luzAmbiente.b * Ka.b);
    if (naSombra) {
        Color I(lidarExcecao(Ia.r), lidarExcecao(Ia.g), lidarExcecao(Ia.b));
        int R = static_cast<int>(roundf(I.r * 255.0f));
        int G = static_cast<int>(roundf(I.g * 255.0f));
        int B = static_cast<int>(roundf(I.b * 255.0f));
        return Color(R, G, B);
    }

    Vector n = CalcularNormal(P);
    Vector v(-dir.i, -dir.j, -dir.k);
    Vector r1 = calcula_esc_por_vetor(2.0f, calcula_esc_por_vetor(calcula_prod_esc(n, l), n));
    Vector r(r1.i - l.i, r1.j - l.j, r1.k - l.k);

    float fd = lidarExcecao(calcula_prod_esc(n, l));
    fd = std::clamp(fd, 0.0f, 1.0f);

    float cosAlpha = lidarExcecao(calcula_prod_esc(r, v));
    cosAlpha = std::clamp(cosAlpha, 0.0f, 1.0f);

    float fe = pow(cosAlpha, static_cast<float>(getShininess()));

    Color Id(cena.luz.intensidade.r * Kd.r * fd,
        cena.luz.intensidade.g * Kd.g * fd,
        cena.luz.intensidade.b * Kd.b * fd);
    Color Ie(cena.luz.intensidade.r * Ke.r * fe,
        cena.luz.intensidade.g * Ke.g * fe,
        cena.luz.intensidade.b * Ke.b * fe);
    Color I(lidarExcecao(Id.r + Ie.r + Ia.r), lidarExcecao(Id.g + Ie.g + Ia.g), lidarExcecao(Id.b + Ie.b + Ia.b));

    int R = static_cast<int>(roundf(I.r * 255.0f));
    int G = static_cast<int>(roundf(I.g * 255.0f));
    int B = static_cast<int>(roundf(I.b * 255.0f));
    return Color(R, G, B);
}

bool Esfera::verificarIntersecao(Vetor p0, Vetor dr) {
    Point origem(p0.i, p0.j, p0.k);
    const Vector dir(dr.i, dr.j, dr.k);
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

Vetor Esfera::calcularNormal(Vetor posicao) {
    Point P(posicao.i, posicao.j, posicao.k);
    Vector n = subtrai_pontos(P, centro);
    float nn = calcula_norma(n);
    if (nn <= 1e-8f) return Vetor(0.0f, 0.0f, 0.0f);
    return Vetor(n.i / nn, n.j / nn, n.k / nn);
}

Vector calcula_n_esfera(Point P, Esfera esfera) {
    return esfera.CalcularNormal(P);
}

float calcula_t_esfera(Esfera esfera, Vector dr, Point Po) {
    return esfera.CalcularIntersecao(Po, dr);
}
