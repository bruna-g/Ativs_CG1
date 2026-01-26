#include "../include/Cone.h"
#include "../include/Cena.h"
#include "../include/Color.h"
#include "../include/Cilindro.h"
#include "../include/Matriz3x3.h"
#include "../include/Vector.h"
#include "../include/Utils.h"
#include "../include/Matriz4x4.h"
#include <cmath>
#include <algorithm>
#include <limits>

Cone::Cone(const Point& cb_c, const Point& v_c, const float raio_c)
    : cb(cb_c), v(v_c), raio(raio_c), material(), dc(0.0f, 1.0f, 0.0f), altura(0.0f), teta(0.0f) {
    recalcularDerivados();
}

Cone::Cone(const Point& cb_c, const Point& v_c, const float raio_c, const Material& material_p)
    : Cone(cb_c, v_c, raio_c) {
    material = material_p;
}

void Cone::rotacionarX(float angulo) {
    Point pivo = this->cb;
    // Transladar para a origem
    this->cb = this->aplicarTranslacao(this->cb, Vector(-pivo.x, -pivo.y, -pivo.z));
    this->v = this->aplicarTranslacao(this->v, Vector(-pivo.x, -pivo.y, -pivo.z));
    Matriz4x4 R = Matriz4x4::rotacaoX(angulo);

    // Aplicar rotação na base e no eixo
    this->cb = R.multiplicarPonto(this->cb);
    this->v = R.multiplicarPonto(this->v);
    this->dc = R.multiplicarVetor(this->dc);
    this->dc = normalizar(this->dc);
    // Transladar de volta para a posição original
    this->cb = this->aplicarTranslacao(this->cb, Vector(pivo.x, pivo.y, pivo.z));
    this->v = this->aplicarTranslacao(this->v, Vector(pivo.x, pivo.y, pivo.z));
}

void Cone::rotacionarY(float angulo) {
    Point pivo = this->cb;
    // Transladar para a origem
    this->cb = this->aplicarTranslacao(this->cb, Vector(-pivo.x, -pivo.y, -pivo.z));
    this->v = this->aplicarTranslacao(this->v, Vector(-pivo.x, -pivo.y, -pivo.z));
    Matriz4x4 R = Matriz4x4::rotacaoY(angulo);

    // Aplicar rotação na base e no eixo
    this->cb = R.multiplicarPonto(this->cb);
    this->v = R.multiplicarPonto(this->v);
    this->dc = R.multiplicarVetor(this->dc);
    this->dc = normalizar(this->dc);
    // Transladar de volta para a posição original
    this->cb = this->aplicarTranslacao(this->cb, Vector(pivo.x, pivo.y, pivo.z));
    this->v = this->aplicarTranslacao(this->v, Vector(pivo.x, pivo.y, pivo.z));
}

void Cone::rotacionarZ(float angulo) {
    Point pivo = this->cb;
    // Transladar para a origem
    this->cb = this->aplicarTranslacao(this->cb, Vector(-pivo.x, -pivo.y, -pivo.z));
    this->v = this->aplicarTranslacao(this->v, Vector(-pivo.x, -pivo.y, -pivo.z));
    Matriz4x4 R = Matriz4x4::rotacaoZ(angulo);

    // Aplicar rotação na base e no eixo
    this->cb = R.multiplicarPonto(this->cb);
    this->v = R.multiplicarPonto(this->v);
    this->dc = R.multiplicarVetor(this->dc);
    this->dc = normalizar(this->dc);
    // Transladar de volta para a posição original
    this->cb = this->aplicarTranslacao(this->cb, Vector(pivo.x, pivo.y, pivo.z));
    this->v = this->aplicarTranslacao(this->v, Vector(pivo.x, pivo.y, pivo.z));
}

void Cone::rotacaoArbitraria(Vector eixo, float angulo) {
    Point pivo = this->cb;
    // Transladar para a origem
    this->cb = this->aplicarTranslacao(this->cb, Vector(-pivo.x, -pivo.y, -pivo.z));
    this->v = this->aplicarTranslacao(this->v, Vector(-pivo.x, -pivo.y, -pivo.z));
    Matriz4x4 R = Matriz4x4::rotacaoArbitraria(eixo, angulo);

    // Aplicar rotação na base e no eixo
    this->cb = R.multiplicarPonto(this->cb);
    this->v = R.multiplicarPonto(this->v);
    this->dc = R.multiplicarVetor(this->dc);
    this->dc = normalizar(this->dc);
    // Transladar de volta para a posição original
    this->cb = this->aplicarTranslacao(this->cb, Vector(pivo.x, pivo.y, pivo.z));
    this->v = this->aplicarTranslacao(this->v, Vector(pivo.x, pivo.y, pivo.z));
}

void Cone::recalcularDerivados() {
    Vector cbv = subtrai_pontos(v, cb);
    altura = calcula_norma(cbv);
    if (altura > 0.0f) {
        dc = Vector(cbv.i / altura, cbv.j / altura, cbv.k / altura);
    }
    else {
        dc = Vector(0.0f, 1.0f, 0.0f);
        altura = 0.0f;
    }
    teta = atan(raio / (altura > 0.0f ? altura : 1.0f));
}

void Cone::aplicarEscalaUniforme(float s) {
    if (s <= 0.0f) return;

    // Mantém a base (cb) fixa e escala o cone em torno de cb.
    // Isso muda o tamanho (altura e raio) sem "arrastar" o objeto pelo mundo.
    raio *= s;

    Vector cbv = subtrai_pontos(v, cb);
    Vector cbvEsc = calcula_esc_por_vetor(s, cbv);
    v = Point(cb.x + cbvEsc.i, cb.y + cbvEsc.j, cb.z + cbvEsc.k);

    recalcularDerivados();
}

void Cone::aplicarEscalaNoPivoObjeto(const Vetor& escala, const Point& pivo) {
    // 1) Escala os pontos que definem o cone.
    cb = Objeto::aplicarEscalaNoPivo(cb, escala, pivo);
    v = Objeto::aplicarEscalaNoPivo(v, escala, pivo);

    // 2) Atualiza raio com um fator radial aproximado no plano perpendicular ao eixo.
    Vector eixo = subtrai_pontos(v, cb);
    float alt = calcula_norma(eixo);
    Vector dcLocal(0.0f, 1.0f, 0.0f);
    if (alt > 1e-6f) {
        dcLocal = Vector(eixo.i / alt, eixo.j / alt, eixo.k / alt);
    }

    // Base ortonormal (u,v) perpendicular ao eixo para medir o "alongamento" radial.
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

    // 3) Recalcula parâmetros derivados (dc, altura, teta).
    recalcularDerivados();
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

    naSombraCil = cena.estaNaSombra(Pi_mod, l, dist_Pi_Luz, const_cast<Cone*>(this));

    if (naSombraCil) {
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
    bool naSombraCil = false;

    Vetor kaV = getKa();
    Vetor kdV = getKd();
    Vetor keV = getKe();
    Color Ka(kaV.i, kaV.j, kaV.k);
    Color Kd(kdV.i, kdV.j, kdV.k);
    Color Ke(keV.i, keV.j, keV.k);

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
    Color Ia(cena.luzAmbiente.r * Ka.r, cena.luzAmbiente.g * Ka.g, cena.luzAmbiente.b * Ka.b);

    Point Pi_mod(Pi.x + l.i * 1e-4f, Pi.y + l.j * 1e-4f, Pi.z + l.k * 1e-4f);
    float dist_Pi_Luz = calcula_norma(subtrai_pontos(cena.luz.pos, Pi_mod));

    naSombraCil = cena.estaNaSombra(Pi_mod, l, dist_Pi_Luz, const_cast<Cone*>(this));

    if (naSombraCil) {
        Color I(lidarExcecao(Ia.r), lidarExcecao(Ia.g), lidarExcecao(Ia.b));
        int R = static_cast<int>(roundf(I.r * 255.0f));
        int G = static_cast<int>(roundf(I.g * 255.0f));
        int B = static_cast<int>(roundf(I.b * 255.0f));
        return Color(R, G, B);
    }

    float fd = lidarExcecao(calcula_prod_esc(n, l));
    float cosAlpha = lidarExcecao(calcula_prod_esc(r, vdir));
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

Vector Cone::CalcularNormal(const Point& Pi, const Vector& dir) const {
    // Base do cone: plano em cb com normal dc
    Vector dist_centro = subtrai_pontos(Pi, cb);
    float altura_Pi = calcula_prod_esc(dist_centro, dc);
    bool na_base = fabs(altura_Pi) < 1e-3f;

    return na_base
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
}

// Wrapper legado (mantido para compatibilidade com código antigo)
extern Point Po;
extern Cone cone;
float calcula_t_cone(Point& Pj) {
    Vector dir = calcula_dr(Po, Pj);
    return cone.CalcularIntersecao(Po, dir);
}

bool Cone::verificarIntersecao(Vetor p0, Vetor dr) {
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

Vetor Cone::calcularNormal(Vetor posicao) {
    Point Pi(posicao.i, posicao.j, posicao.k);

    // Detecta base: plano passando em cb com normal dc
    Vector dist_centro = subtrai_pontos(Pi, cb);
    float altura_Pi = calcula_prod_esc(dist_centro, dc);
    bool na_base = std::fabs(altura_Pi) < 1e-3f;

    Vector nLocal = na_base
        ? dc
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

    return Vetor(nLocal.i, nLocal.j, nLocal.k);
}
