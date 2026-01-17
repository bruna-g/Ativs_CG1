#ifndef MALHA_H
#define MALHA_H

#include "Objeto.h"
#include "AlgebraLinear.h"

#include <vector>
#include <initializer_list>

struct Vertice {
    Vetor pos;

    explicit Vertice(const Vetor& p) : pos(p) {}
};

struct Aresta {
    int v0;
    int v1;

    Aresta(int a, int b) : v0(a), v1(b) {}
};

struct Face {
    std::vector<int> indices;
    Vetor normal;

    Face() : indices(), normal(0.0f, 0.0f, 0.0f) {}
    Face(int a, int b, int c) : indices{ a, b, c }, normal(0.0f, 0.0f, 0.0f) {}
    Face(std::initializer_list<int> ids) : indices(ids), normal(0.0f, 0.0f, 0.0f) {}
};

class Malha : public Objeto {
public:
    std::vector<Aresta> arestas;
    std::vector<Face> faces;
    std::vector<Vertice> vertices;

    Vetor normal;

    AlgebraLinear al;

    Malha();

    void adicionarAresta(Aresta aresta);
    void adicionarFace(Face face);
    void adicionarVertice(Vertice vertice);

    bool verificarIntersecao(Vetor p0, Vetor dr) override;
    Vetor calcularNormal(Vetor posicao) override;

private:
    bool rayTriangleIntersect(const Vetor& origem, const Vetor& dir,
        const Vetor& v0, const Vetor& v1, const Vetor& v2,
        float& t) const;

    Vetor calcularNormalFace(const Face& face) const;
    void atualizarNormalGlobal();
};

#endif // MALHA_H
