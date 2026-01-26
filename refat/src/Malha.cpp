#include "../include/Malha.h"
#include "../include/Utils.h"
#include "../include/Textura.hpp"
#include <cmath>
#include <limits>
#include <Matriz4x4.h>

Malha::Malha() : arestas(), faces(), vertices(), normal(0.0f, 0.0f, 0.0f), al(), plano(Point(0.0f, 0.0f, 0.0f), Vector(0.0f, 1.0f, 0.0f), Material()) {
}

Malha::Malha(const Plano& plano) : Malha() {
    this->plano = plano;
}

void Malha::adicionarAresta(Aresta aresta) {
    arestas.push_back(aresta);
}

void Malha::adicionarVertice(Vertice vertice) {
    vertices.push_back(vertice);
}

void Malha::adicionarFace(Face face) {
    // Tenta calcular normal se a face tiver pelo menos 3 índices válidos.
    if (face.indices.size() >= 3) {
        bool indicesValidos = true;
        for (int idx : face.indices) {
            if (idx < 0 || static_cast<size_t>(idx) >= vertices.size()) {
                indicesValidos = false;
                break;
            }
        }
        if (indicesValidos) {
            face.normal = calcularNormalFace(face);
        }
    }

    faces.push_back(face);
    atualizarNormalGlobal();
}

void Malha::aplicarEscalaUniforme(float s) {
    if (s <= 0.0f) return;

    for (auto& vert : vertices) {
        Vetor& p = vert.pos;
        p = Vetor(
            p.i * static_cast<float>(s),
            p.j * static_cast<float>(s),
            p.k * static_cast<float>(s),
            p.q
        );
    }

    for (auto& face : faces) {
        if (face.indices.size() >= 3) {
            face.normal = calcularNormalFace(face);
        }
    }
    atualizarNormalGlobal();
}

void Malha::rotacionarX(float angulo) {
    Point centro = calcularCentro();
    Matriz4x4 R = Matriz4x4::rotacaoX(angulo);

    for (auto& vertice : vertices) {
        Vetor& v = vertice.pos;
        Point p(v.i, v.j, v.k);
        p = aplicarTranslacao(p, Vector(-centro.x, -centro.y, -centro.z));

        p = R.multiplicarPonto(p);

        p = aplicarTranslacao(p, Vector(centro.x, centro.y, centro.z));

        vertice.pos.i = p.x;
        vertice.pos.j = p.y;
        vertice.pos.k = p.z;
    }

    for (auto& face : faces) {
        if (face.indices.size() >= 3) {
            face.normal = calcularNormalFace(face);
        }
    }
    atualizarNormalGlobal();
}

void Malha::rotacionarY(float angulo) {
    Point centro = calcularCentro();
    Matriz4x4 R = Matriz4x4::rotacaoY(angulo);

    for (auto& vertice : vertices) {
        Vetor& v = vertice.pos;
        Point p(v.i, v.j, v.k);
        p = aplicarTranslacao(p, Vector(-centro.x, -centro.y, -centro.z));

        p = R.multiplicarPonto(p);

        p = aplicarTranslacao(p, Vector(centro.x, centro.y, centro.z));

        vertice.pos.i = p.x;
        vertice.pos.j = p.y;
        vertice.pos.k = p.z;
    }

    for (auto& face : faces) {
        if (face.indices.size() >= 3) {
            face.normal = calcularNormalFace(face);
        }
    }
    atualizarNormalGlobal();
}

void Malha::rotacionarZ(float angulo) {
    Point centro = calcularCentro();
    Matriz4x4 R = Matriz4x4::rotacaoZ(angulo);

    for (auto& vertice : vertices) {
        Vetor& v = vertice.pos;
        Point p(v.i, v.j, v.k);
        p = aplicarTranslacao(p, Vector(-centro.x, -centro.y, -centro.z));

        p = R.multiplicarPonto(p);

        p = aplicarTranslacao(p, Vector(centro.x, centro.y, centro.z));

        vertice.pos.i = p.x;
        vertice.pos.j = p.y;
        vertice.pos.k = p.z;
    }

    for (auto& face : faces) {
        if (face.indices.size() >= 3) {
            face.normal = calcularNormalFace(face);
        }
    }
    atualizarNormalGlobal();
}

void Malha::rotacaoArbitraria(Vector eixo, float angulo) {
    Point centro = calcularCentro();
    Matriz4x4 R = Matriz4x4::rotacaoArbitraria(eixo, angulo);

    for (auto& vertice : vertices) {
        Vetor& v = vertice.pos;
        Point p(v.i, v.j, v.k);
        p = aplicarTranslacao(p, Vector(-centro.x, -centro.y, -centro.z));

        p = R.multiplicarPonto(p);

        p = aplicarTranslacao(p, Vector(centro.x, centro.y, centro.z));

        vertice.pos.i = p.x;
        vertice.pos.j = p.y;
        vertice.pos.k = p.z;
    }

    for (auto& face : faces) {
        if (face.indices.size() >= 3) {
            face.normal = calcularNormalFace(face);
        }
    }
    atualizarNormalGlobal();
}

Point Malha::calcularCentro() {
    if (this->vertices.empty()) {
        return Point(0, 0, 0);
    }

    float somaX = 0, somaY = 0, somaZ = 0;
    for (const auto& vertice : vertices) {
        somaX += vertice.pos.i;
        somaY += vertice.pos.j;
        somaZ += vertice.pos.k;
    }

    int n = vertices.size();
    return Point(somaX / n, somaY / n, somaZ / n);
}

bool Malha::rayTriangleIntersect(const Vetor& origem, const Vetor& dir,
    const Vetor& v0, const Vetor& v1, const Vetor& v2,
    float& t) const {
    // Möller–Trumbore
    constexpr float EPS = 1e-6f;

    Vetor e1 = al.subtrair(v1, v0);
    Vetor e2 = al.subtrair(v2, v0);

    Vetor pvec = al.produtoVetorial(dir, e2);
    float det = al.produtoEscalar(e1, pvec);

    if (std::fabs(det) < EPS) return false;

    float invDet = 1.0f / det;

    Vetor tvec = al.subtrair(origem, v0);
    float u = al.produtoEscalar(tvec, pvec) * invDet;
    if (u < 0.0f || u > 1.0f) return false;

    Vetor qvec = al.produtoVetorial(tvec, e1);
    float v = al.produtoEscalar(dir, qvec) * invDet;
    if (v < 0.0f || (u + v) > 1.0f) return false;

    float tt = al.produtoEscalar(e2, qvec) * invDet;
    if (tt <= EPS) return false;

    t = tt;
    return true;
}

Vetor Malha::calcularNormalFace(const Face& face) const {
    if (face.indices.size() < 3) return Vetor(0.0f, 0.0f, 0.0f);

    int i0 = face.indices[0];
    int i1 = face.indices[1];
    int i2 = face.indices[2];
    if (i0 < 0 || i1 < 0 || i2 < 0) return Vetor(0.0f, 0.0f, 0.0f);
    if (static_cast<size_t>(i0) >= vertices.size() ||
        static_cast<size_t>(i1) >= vertices.size() ||
        static_cast<size_t>(i2) >= vertices.size()) {
        return Vetor(0.0f, 0.0f, 0.0f);
    }

    const Vetor& v0 = vertices[i0].pos;
    const Vetor& v1 = vertices[i1].pos;
    const Vetor& v2 = vertices[i2].pos;

    Vetor e1 = al.subtrair(v1, v0);
    Vetor e2 = al.subtrair(v2, v0);
    Vetor n = al.produtoVetorial(e1, e2);
    return al.normalizar(n);
}

void Malha::atualizarNormalGlobal() {
    Vetor soma(0.0f, 0.0f, 0.0f);

    for (const Face& f : faces) {
        Vetor fn = f.normal;
        // Se não tiver normal pré-computada, tenta calcular agora.
        if (al.norma(fn) <= 1e-8f) {
            fn = calcularNormalFace(f);
        }
        soma = al.somar(soma, fn);
    }

    normal = al.normalizar(soma);
}

bool Malha::verificarIntersecao(Vetor p0, Vetor dr) {
    // Calcula a interseção mais próxima e preenche os campos do Objeto:
    // - temIntersecao
    // - distancia (t)
    // - pontoIntersecao
    // Também atualiza `normal` da malha para a normal da face atingida.

    constexpr float EPS = 1e-6f;

    setTemIntersecao(false);
    setDistancia(std::numeric_limits<double>::infinity());
    setPontoIntersecao(Vetor(0.0f, 0.0f, 0.0f));

    float melhorT = std::numeric_limits<float>::infinity();
    Vetor melhorNormal(0.0f, 0.0f, 0.0f);

    for (const Face& face : faces) {
        if (face.indices.size() < 3) continue;

        int i0 = face.indices[0];
        if (i0 < 0 || static_cast<size_t>(i0) >= vertices.size()) continue;
        const Vetor& v0 = vertices[i0].pos;

        for (size_t k = 1; k + 1 < face.indices.size(); ++k) {
            int i1 = face.indices[k];
            int i2 = face.indices[k + 1];
            if (i1 < 0 || i2 < 0) continue;
            if (static_cast<size_t>(i1) >= vertices.size() || static_cast<size_t>(i2) >= vertices.size()) continue;

            const Vetor& v1 = vertices[i1].pos;
            const Vetor& v2 = vertices[i2].pos;

            float t = 0.0f;
            if (!rayTriangleIntersect(p0, dr, v0, v1, v2, t) || t <= EPS) continue;

            if (t < melhorT) {
                melhorT = t;

                // Normal da face (se já tiver, usa; senão calcula)
                Vetor nFace = face.normal;
                if (al.norma(nFace) <= 1e-8f) {
                    Face tmp = face;
                    tmp.normal = Vetor(0.0f, 0.0f, 0.0f);
                    nFace = calcularNormalFace(tmp);
                }

                // Orienta a normal contra o raio (opcional, mas útil para shading)
                if (al.produtoEscalar(nFace, dr) > 0.0f) {
                    nFace = al.escalar(-1.0f, nFace);
                }
                melhorNormal = nFace;
            }
        }
    }

    if (melhorT < std::numeric_limits<float>::infinity()) {
        setTemIntersecao(true);
        setDistancia(static_cast<double>(melhorT));
        Vetor Pi = al.somar(p0, al.escalar(melhorT, dr));
        setPontoIntersecao(Pi);
        normal = melhorNormal;
        return true;
    }

    return false;
}

Vetor Malha::calcularNormal(Vetor /*posicao*/) {
    // Se houve interseção recentemente, `normal` já foi atualizada para a normal da face atingida.
    if (getTemIntersecao() && al.norma(normal) > 1e-8f) {
        return normal;
    }

    // Caso contrário, retorna a normal global (média das faces).
    if (al.norma(normal) <= 1e-8f) {
        atualizarNormalGlobal();
    }
    return normal;
}

void Malha::aplicarEscalaNoPivoObjeto(const Vetor& escala, const Point& pivo) {
    const float sx = escala.i;
    const float sy = escala.j;
    const float sz = escala.k;

    for (auto& vert : vertices) {
        Vetor& p = vert.pos;
        Point pt(p.i, p.j, p.k);
        Point pt2 = Objeto::aplicarEscalaNoPivo(pt, Vetor(sx, sy, sz), pivo);
        p = Vetor(
            pt2.x,
            pt2.y,
            pt2.z,
            p.q
        );
    }

    for (auto& face : faces) {
        if (face.indices.size() >= 3) {
            face.normal = calcularNormalFace(face);
        }
    }
    atualizarNormalGlobal();
}


Color Malha::CalcularCor(const Cena& cena, float t, const Vector& dir) const {
    Vetor kaV = getKa();
    Vetor kdV = getKd();
    Vetor keV = getKe();
    Color K_a(kaV.i, kaV.j, kaV.k);
    Color K_d(kdV.i, kdV.j, kdV.k);
    Color K_e(keV.i, keV.j, keV.k);

    Point P = calcula_eq_ray(cena.observador, t, dir);
    Vetor normal_vec = normal; // Usa a normal da face atingida, calculada na interseção.
    Vector normal_v(normal_vec.i, normal_vec.j, normal_vec.k);

    Color Kd_textured = K_d;
    if (plano.material.usarTextura && plano.material.textura != nullptr) {
        Vector nLocal = normal_v;
        Vector arbitrary = (fabs(nLocal.i) < 0.9f) ? Vector(1.0f, 0.0f, 0.0f) : Vector(0.0f, 1.0f, 0.0f);
        Vector u_axis = cross(arbitrary, nLocal);
        float nu = calcula_norma(u_axis);
        if (nu == 0.0f) nu = 1.0f;
        u_axis = calcula_esc_por_vetor(1.0f / nu, u_axis);
        Vector v_axis = cross(nLocal, u_axis);

        Vector vecPi = subtrai_pontos(P, plano.p_pi);
        float u = calcula_prod_esc(vecPi, u_axis) * plano.material.texturaScale;
        float v = calcula_prod_esc(vecPi, v_axis) * plano.material.texturaScale;
        u = u - floor(u);
        if (u < 0) u += 1.0f;
        v = v - floor(v);
        if (v < 0) v += 1.0f;

        size_t tx = static_cast<size_t>(u * plano.material.textura->get_largura_pixels()) % plano.material.textura->get_largura_pixels();
        size_t ty = static_cast<size_t>(v * plano.material.textura->get_altura_pixels()) % plano.material.textura->get_altura_pixels();

        rgb px = plano.material.textura->get_cor_pixel(ty, tx);
        Color texCol(px[0] / 255.0f, px[1] / 255.0f, px[2] / 255.0f);
        Kd_textured = Color(texCol.r, texCol.g, texCol.b);
    }

    bool naSombraSpot = false;
    Color I_spot(0.0f, 0.0f, 0.0f);

    // --- Luz Spot ---
    if (cena.luzSpotAtiva) {
        Vector l_spot = calcula_l(cena.luzSpot.pos, P);
        Point P_mod_spot(P.x + l_spot.i * 1e-4f, P.y + l_spot.j * 1e-4f, P.z + l_spot.k * 1e-4f);
        float dist_Pi_Luz_spot = calcula_norma(subtrai_pontos(cena.luzSpot.pos, P_mod_spot));

        naSombraSpot = cena.estaNaSombra(P_mod_spot, l_spot, dist_Pi_Luz_spot, const_cast<Malha*>(this));

        if (!naSombraSpot) {
            Vector r1_spot = calcula_esc_por_vetor(2.0f, calcula_esc_por_vetor(calcula_prod_esc(normal_v, l_spot), normal_v));
            Vector r_spot(r1_spot.i - l_spot.i, r1_spot.j - l_spot.j, r1_spot.k - l_spot.k);

            float fd_spot = lidarExcecao(calcula_prod_esc(normal_v, l_spot));
            float cosAlpha_spot = lidarExcecao(calcula_prod_esc(r_spot, Vector(-dir.i, -dir.j, -dir.k)));
            float fe_spot = pow(cosAlpha_spot, static_cast<float>(getShininess()));

            Vector l_spot_neg = calcula_esc_por_vetor(-1.0f, l_spot);
            float cos_dr_l = calcula_prod_esc(cena.luzSpot.direcao, l_spot_neg);
            float cos_ang = cosf(cena.luzSpot.angulo * (3.14159265f / 180.0f));

            if (cos_dr_l >= cos_ang) {
                Color Id_spot(cena.luzSpot.intensidade.r * cos_dr_l * Kd_textured.r * fd_spot,
                    cena.luzSpot.intensidade.g * cos_dr_l * Kd_textured.g * fd_spot,
                    cena.luzSpot.intensidade.b * cos_dr_l * Kd_textured.b * fd_spot);

                Color Ie_spot(cena.luzSpot.intensidade.r * cos_dr_l * K_e.r * fe_spot,
                    cena.luzSpot.intensidade.g * cos_dr_l * K_e.g * fe_spot,
                    cena.luzSpot.intensidade.b * cos_dr_l * K_e.b * fe_spot);

                I_spot = Color(Id_spot.r + Ie_spot.r, Id_spot.g + Ie_spot.g, Id_spot.b + Ie_spot.b);
            }
        }
    }

    // --- Luz Pontual ---
    Vector l_pontual = calcula_l(cena.luz.pos, P);
    Point P_mod_pontual(P.x + l_pontual.i * 1e-4f, P.y + l_pontual.j * 1e-4f, P.z + l_pontual.k * 1e-4f);
    float dist_Pi_Luz_pontual = calcula_norma(subtrai_pontos(cena.luz.pos, P_mod_pontual));

    bool naSombraPontual = cena.estaNaSombra(P_mod_pontual, l_pontual, dist_Pi_Luz_pontual, const_cast<Malha*>(this));

    Color I_pontual(0.0f, 0.0f, 0.0f);
    if (!naSombraPontual) {
        Vector v(-dir.i, -dir.j, -dir.k);
        Vector r1_pontual = calcula_esc_por_vetor(2.0f, calcula_esc_por_vetor(calcula_prod_esc(normal_v, l_pontual), normal_v));
        Vector r_pontual(r1_pontual.i - l_pontual.i, r1_pontual.j - l_pontual.j, r1_pontual.k - l_pontual.k);

        float fd_pontual = lidarExcecao(calcula_prod_esc(normal_v, l_pontual));
        float cosAlpha_pontual = lidarExcecao(calcula_prod_esc(r_pontual, v));
        float fe_pontual = pow(cosAlpha_pontual, static_cast<float>(getShininess()));

        Color Id_pontual(cena.luz.intensidade.r * Kd_textured.r * fd_pontual,
            cena.luz.intensidade.g * Kd_textured.g * fd_pontual,
            cena.luz.intensidade.b * Kd_textured.b * fd_pontual);

        Color Ie_pontual(cena.luz.intensidade.r * K_e.r * fe_pontual,
            cena.luz.intensidade.g * K_e.g * fe_pontual,
            cena.luz.intensidade.b * K_e.b * fe_pontual);

        I_pontual = Color(Id_pontual.r + Ie_pontual.r, Id_pontual.g + Ie_pontual.g, Id_pontual.b + Ie_pontual.b);
    }

    // --- Cor Final ---
    Color Ia(cena.luzAmbiente.r * K_a.r, cena.luzAmbiente.g * K_a.g, cena.luzAmbiente.b * K_a.b);

    Color I(lidarExcecao(I_pontual.r + I_spot.r + Ia.r),
        lidarExcecao(I_pontual.g + I_spot.g + Ia.g),
        lidarExcecao(I_pontual.b + I_spot.b + Ia.b));

    int R = static_cast<int>(roundf(I.r * 255.0f));
    int G = static_cast<int>(roundf(I.g * 255.0f));
    int B = static_cast<int>(roundf(I.b * 255.0f));
    return Color(R, G, B);
}
