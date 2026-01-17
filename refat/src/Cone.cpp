#include "../include/Cone.h"
#include "../include/Cena.h"
#include "../include/Color.h"
#include "../include/Cilindro.h"
#include "../include/Matriz3x3.h"
#include "../include/Vector.h"
#include "../include/Utils.h"
#include <cmath>
#include <algorithm>

Cone::Cone(const Point& cb_c, const Point& v_c, const float raio_c)
    : cb(cb_c), v(v_c), raio(raio_c), material(), dc(0.0f, 1.0f, 0.0f), altura(0.0f), teta(0.0f) {
    Vector cbv = subtrai_pontos(v, cb);
    altura = calcula_norma(cbv);
    if (altura > 0.0f) {
        dc = Vector(cbv.i / altura, cbv.j / altura, cbv.k / altura);
    }
    teta = atan(raio / (altura > 0.0f ? altura : 1.0f));
}

Cone::Cone(const Point& cb_c, const Point& v_c, const float raio_c, const Material& material_p)
    : Cone(cb_c, v_c, raio_c) {
    material = material_p;
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

Color Cone::CalcularCor(const Cena& cena, float t, const Vector& dir, const Color& K) const {
    bool naSombraEsf = false;
    bool naSombraCil = false;

    Point Pi = calcula_eq_ray(cena.observador, t, dir);

    // Base do cone: plano em cb com normal dc
    Vector dist_centro = subtrai_pontos(Pi, cb);
    float altura_Pi = calcula_prod_esc(dist_centro, dc);
    bool na_base = fabs(altura_Pi) < 1e-3f;

    Vector n = na_base
        ? Vector(
            (calcula_prod_esc(dc, dir) > 0 ? -dc.i : dc.i),
            (calcula_prod_esc(dc, dir) > 0 ? -dc.j : dc.j),
            (calcula_prod_esc(dc, dir) > 0 ? -dc.k : dc.k))
        : [&]() {
        Vector V = subtrai_pontos(v, Pi);
        float vNorma = calcula_norma(V);
        if (vNorma == 0.0f) vNorma = 1.0f;
        Vector s_conjugado(V.i / vNorma, V.j / vNorma, V.k / vNorma);
        Matriz3x3 M_id(1.0f, 0.0f, 0.0f,
            0.0f, 1.0f, 0.0f,
            0.0f, 0.0f, 1.0f);
        Matriz3x3 M_e = matrizSubtrai(M_id, outerProduto(s_conjugado));
        Vector N = matrizVetorProduto(M_e, dc);
        float N_norma = calcula_norma(N);
        if (N_norma == 0.0f) N_norma = 1.0f;
        return Vector(N.i / N_norma, N.j / N_norma, N.k / N_norma);
        }();

    Vector l = calcula_l(cena.luz.pos, Pi);
    Vector vdir(-dir.i, -dir.j, -dir.k);
    Vector r1 = calcula_esc_por_vetor(2.0f, calcula_esc_por_vetor(calcula_prod_esc(n, l), n));
    Vector r(r1.i - l.i, r1.j - l.j, r1.k - l.k);
    Color Ia(cena.luzAmbiente.r * K.r, cena.luzAmbiente.g * K.g, cena.luzAmbiente.b * K.b);

    Point Pi_mod(Pi.x + l.i * 1e-4f, Pi.y + l.j * 1e-4f, Pi.z + l.k * 1e-4f);
    float dist_Pi_Luz = calcula_norma(subtrai_pontos(cena.luz.pos, Pi_mod));

    // Sombra da esfera
    Vector w_sombra = subtrai_pontos(Pi_mod, cena.centroEsfera);
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

    // Sombra do cilindro
    if (!naSombraEsf && cena.cilindro != nullptr) {
        float t_cil_shadow = cena.cilindro->CalcularIntersecao(Pi_mod, l);
        if (t_cil_shadow > 1e-4f && t_cil_shadow < dist_Pi_Luz) naSombraCil = true;
    }

    if (naSombraCil || naSombraEsf) {
        Color I(lidarExcecao(Ia.r), lidarExcecao(Ia.g), lidarExcecao(Ia.b));
        int R = static_cast<int>(roundf(I.r * 255.0f));
        int G = static_cast<int>(roundf(I.g * 255.0f));
        int B = static_cast<int>(roundf(I.b * 255.0f));
        return Color(R, G, B);
    }

    float fd = lidarExcecao(calcula_prod_esc(n, l));
    float cosAlpha = lidarExcecao(calcula_prod_esc(r, vdir));
    float fe = pow(cosAlpha, cena.expoenteEspecular);

    Color Id(cena.luz.intensidade.r * K.r * fd,
        cena.luz.intensidade.g * K.g * fd,
        cena.luz.intensidade.b * K.b * fd);
    Color Ie(cena.luz.intensidade.r * K.r * fe,
        cena.luz.intensidade.g * K.g * fe,
        cena.luz.intensidade.b * K.b * fe);
    Color I(lidarExcecao(Id.r + Ie.r + Ia.r), lidarExcecao(Id.g + Ie.g + Ia.g), lidarExcecao(Id.b + Ie.b + Ia.b));

    int R = static_cast<int>(roundf(I.r * 255.0f));
    int G = static_cast<int>(roundf(I.g * 255.0f));
    int B = static_cast<int>(roundf(I.b * 255.0f));
    return Color(R, G, B);
}

Color Cone::CalcularCor(const Cena& cena, float t, const Vector& dir) const {
    bool naSombraEsf = false;
    bool naSombraCil = false;

    Point Pi = calcula_eq_ray(cena.observador, t, dir);

    // Base do cone: plano em cb com normal dc
    Vector dist_centro = subtrai_pontos(Pi, cb);
    float altura_Pi = calcula_prod_esc(dist_centro, dc);
    bool na_base = fabs(altura_Pi) < 1e-3f;

    Vector n = na_base
        ? Vector(
            (calcula_prod_esc(dc, dir) > 0 ? -dc.i : dc.i),
            (calcula_prod_esc(dc, dir) > 0 ? -dc.j : dc.j),
            (calcula_prod_esc(dc, dir) > 0 ? -dc.k : dc.k))
        : [&]() {
        Vector V = subtrai_pontos(v, Pi);
        float vNorma = calcula_norma(V);
        if (vNorma == 0.0f) vNorma = 1.0f;
        Vector s_conjugado(V.i / vNorma, V.j / vNorma, V.k / vNorma);
        Matriz3x3 M_id(1.0f, 0.0f, 0.0f,
            0.0f, 1.0f, 0.0f,
            0.0f, 0.0f, 1.0f);
        Matriz3x3 M_e = matrizSubtrai(M_id, outerProduto(s_conjugado));
        Vector N = matrizVetorProduto(M_e, dc);
        float N_norma = calcula_norma(N);
        if (N_norma == 0.0f) N_norma = 1.0f;
        return Vector(N.i / N_norma, N.j / N_norma, N.k / N_norma);
        }();

    Vector l = calcula_l(cena.luz.pos, Pi);
    Vector vdir(-dir.i, -dir.j, -dir.k);
    Vector r1 = calcula_esc_por_vetor(2.0f, calcula_esc_por_vetor(calcula_prod_esc(n, l), n));
    Vector r(r1.i - l.i, r1.j - l.j, r1.k - l.k);
    Color Ia(cena.luzAmbiente.r * material.Ka.r, cena.luzAmbiente.g * material.Ka.g, cena.luzAmbiente.b * material.Ka.b);

    Point Pi_mod(Pi.x + l.i * 1e-4f, Pi.y + l.j * 1e-4f, Pi.z + l.k * 1e-4f);
    float dist_Pi_Luz = calcula_norma(subtrai_pontos(cena.luz.pos, Pi_mod));

    // Sombra da esfera
    Vector w_sombra = subtrai_pontos(Pi_mod, cena.centroEsfera);
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

    // Sombra do cilindro
    if (!naSombraEsf && cena.cilindro != nullptr) {
        float t_cil_shadow = cena.cilindro->CalcularIntersecao(Pi_mod, l);
        if (t_cil_shadow > 1e-4f && t_cil_shadow < dist_Pi_Luz) naSombraCil = true;
    }

    if (naSombraCil || naSombraEsf) {
        Color I(lidarExcecao(Ia.r), lidarExcecao(Ia.g), lidarExcecao(Ia.b));
        int R = static_cast<int>(roundf(I.r * 255.0f));
        int G = static_cast<int>(roundf(I.g * 255.0f));
        int B = static_cast<int>(roundf(I.b * 255.0f));
        return Color(R, G, B);
    }

    float fd = lidarExcecao(calcula_prod_esc(n, l));
    float cosAlpha = lidarExcecao(calcula_prod_esc(r, vdir));
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

// Wrapper legado (mantido para compatibilidade com código antigo)
extern Point Po;
extern Cone cone;
float calcula_t_cone(Point& Pj) {
    Vector dir = calcula_dr(Po, Pj);
    return cone.CalcularIntersecao(Po, dir);
}
