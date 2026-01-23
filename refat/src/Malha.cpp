#include "../include/Malha.h"
#include <cmath>
#include <limits>

Malha::Malha() : arestas(), faces(), vertices(), normal(0.0f, 0.0f, 0.0f), al() {
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
