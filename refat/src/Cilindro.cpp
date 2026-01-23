#include "../include/Cilindro.h"
#include "../include/Cena.h"
#include "../include/Color.h"
#include "../include/Utils.h"
#include "../include/Matriz4x4.h"
#include <cmath>
#include <algorithm>
#include <limits>

Cilindro::Cilindro(const Point& cb_c, const float raio_c, const float altura_c, const Vector& dc_c)
    : cb(cb_c), raio(raio_c), altura(altura_c), dc(dc_c), material() {
}

Cilindro::Cilindro(const Point& cb_c, const float raio_c, const float altura_c, const Vector& dc_c, const Material& material_p)
    : cb(cb_c), raio(raio_c), altura(altura_c), dc(dc_c), material(material_p) {
}

void Cilindro::aplicarEscalaUniforme(float s) {
    if (s <= 0.0f) return;
    raio *= s;
    altura *= s;
}

void Cilindro::aplicarEscalaNoPivoObjeto(const Vetor& escala, const Point& pivo) {
    // Topo atual (cb + altura*dc)
    Point topo(cb.x + dc.i * altura, cb.y + dc.j * altura, cb.z + dc.k * altura);

    // 1) Escala base e topo em torno do pivô
    cb = Objeto::aplicarEscalaNoPivo(cb, escala, pivo);
    topo = Objeto::aplicarEscalaNoPivo(topo, escala, pivo);

    // 2) Atualiza eixo (dc) e altura
    Vector eixo = subtrai_pontos(topo, cb);
    float novaAltura = calcula_norma(eixo);
    if (novaAltura > 1e-6f) {
        altura = novaAltura;
        dc = Vector(eixo.i / novaAltura, eixo.j / novaAltura, eixo.k / novaAltura);
    }

    // 3) Atualiza raio por fator radial aproximado no plano perpendicular ao eixo.
    Vector dcLocal = dc;
    Vector ref = (std::fabs(dcLocal.i) < 0.9f) ? Vector(1.0f, 0.0f, 0.0f) : Vector(0.0f, 1.0f, 0.0f);
    Vector u = produto_vetorial(dcLocal, ref);
    float nu = calcula_norma(u);
    if (nu <= 1e-6f) nu = 1.0f;
    u = calcula_esc_por_vetor(1.0f / nu, u);
    Vector vperp = produto_vetorial(dcLocal, u);
    float nv = calcula_norma(vperp);
    if (nv <= 1e-6f) nv = 1.0f;
    vperp = calcula_esc_por_vetor(1.0f / nv, vperp);

    const float sx = escala.i;
    const float sy = escala.j;
    const float sz = escala.k;
    auto lenEsc = [&](const Vector& a) {
        Vector svec(sx * a.i, sy * a.j, sz * a.k);
        return calcula_norma(svec);
        };

    const float fu = lenEsc(u);
    const float fv = lenEsc(vperp);
    const float fatorRadial = 0.5f * (fu + fv);
    if (fatorRadial > 1e-6f) {
        raio *= fatorRadial;
    }
}

Point calcularCentro(Cilindro cilindro) {
    Vetor dcH = calcula_esc_por_vetor(cilindro.altura/2.0f, cilindro.dc);
    Point soma = soma_ponto_vetor(cilindro.cb, dcH);
    return soma;
}

void Cilindro::rotacionarX(float angulo, Cilindro cilindro) {
    Point centro = calcularCentro(cilindro);
    cilindro.cb = cilindro.aplicarTranslacao(cilindro.cb, Vector(-centro.x, -centro.y, -centro.z));
    
    Matriz4x4 R = Matriz4x4::rotacaoX(angulo);
    cilindro.cb = R.multiplicarPonto(cilindro.cb);
    
    cilindro.cb = cilindro.aplicarTranslacao(cilindro.cb, Vector(centro.x, centro.y, centro.z));
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

Color Cilindro::CalcularCor(const Cena& cena, float t, const Vector& dir,
    const Color& K_a, const Color& K_d, const Color& K_e) const {
    Point P = calcula_eq_ray(cena.observador, t, dir);

    Vector l = calcula_l(cena.luz.pos, P);
    Point P_mod(P.x + l.i * 1e-4f, P.y + l.j * 1e-4f, P.z + l.k * 1e-4f);
    float dist_Pi_Luz = calcula_norma(subtrai_pontos(cena.luz.pos, P_mod));
    bool naSombra = cena.estaNaSombra(P_mod, l, dist_Pi_Luz, const_cast<Cilindro*>(this));

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

Color Cilindro::CalcularCor(const Cena& cena, float t, const Vector& dir) const {
    Point P = calcula_eq_ray(cena.observador, t, dir);

    Vetor kaV = getKa();
    Vetor kdV = getKd();
    Vetor keV = getKe();
    Color Ka(kaV.i, kaV.j, kaV.k);
    Color Kd(kdV.i, kdV.j, kdV.k);
    Color Ke(keV.i, keV.j, keV.k);

    Vector l = calcula_l(cena.luz.pos, P);
    Point P_mod(P.x + l.i * 1e-4f, P.y + l.j * 1e-4f, P.z + l.k * 1e-4f);
    float dist_Pi_Luz = calcula_norma(subtrai_pontos(cena.luz.pos, P_mod));
    bool naSombra = cena.estaNaSombra(P_mod, l, dist_Pi_Luz, const_cast<Cilindro*>(this));

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
    float cosAlpha = lidarExcecao(calcula_prod_esc(r, v));
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

Vector calcula_n_cilindro(Point P, Cilindro cilindro) {
    return cilindro.CalcularNormal(P);
}

float calcula_t_Cilindro(Cilindro cilindro, Vector dr, Point Po) {
    return cilindro.CalcularIntersecao(Po, dr);
}

bool Cilindro::verificarIntersecao(Vetor p0, Vetor dr) {
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

Vetor Cilindro::calcularNormal(Vetor posicao) {
    Point P(posicao.i, posicao.j, posicao.k);
    Vector nLocal = CalcularNormal(P);
    return Vetor(nLocal.i, nLocal.j, nLocal.k);
}
