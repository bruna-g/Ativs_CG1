#include "../include/Cubo.h"
#include "../include/Utils.h"
#include "../include/Plano.h"

Malha Cubo::criarCubo(Vetor centroBase, double comprimentoAresta, Vetor Ka, Vetor Kd, Vetor Ke, double shininess) {
    Malha malha;

    malha.setShininess(shininess);
    malha.setKa(Ka);
    malha.setKd(Kd);
    malha.setKe(Ke);

    float half = static_cast<float>(comprimentoAresta / 2.0);
    float a = static_cast<float>(comprimentoAresta);

    // 8 vértices (base em y=centroBase.y e topo em y+aresta)
    malha.adicionarVertice(Vertice(Vetor(centroBase.i - half, centroBase.j, centroBase.k - half, 0.0f)));
    malha.adicionarVertice(Vertice(Vetor(centroBase.i - half, centroBase.j, centroBase.k + half, 1.0f)));
    malha.adicionarVertice(Vertice(Vetor(centroBase.i + half, centroBase.j, centroBase.k + half, 1.0f)));
    malha.adicionarVertice(Vertice(Vetor(centroBase.i + half, centroBase.j, centroBase.k - half, 1.0f)));

    malha.adicionarVertice(Vertice(Vetor(centroBase.i - half, centroBase.j + a, centroBase.k - half, 1.0f)));
    malha.adicionarVertice(Vertice(Vetor(centroBase.i - half, centroBase.j + a, centroBase.k + half, 1.0f)));
    malha.adicionarVertice(Vertice(Vetor(centroBase.i + half, centroBase.j + a, centroBase.k + half, 1.0f)));
    malha.adicionarVertice(Vertice(Vetor(centroBase.i + half, centroBase.j + a, centroBase.k - half, 1.0f)));

    // Arestas (12 básicas)
    malha.adicionarAresta(Aresta(0, 1));
    malha.adicionarAresta(Aresta(1, 2));
    malha.adicionarAresta(Aresta(2, 3));
    malha.adicionarAresta(Aresta(3, 0));

    malha.adicionarAresta(Aresta(4, 5));
    malha.adicionarAresta(Aresta(5, 6));
    malha.adicionarAresta(Aresta(6, 7));
    malha.adicionarAresta(Aresta(7, 4));

    malha.adicionarAresta(Aresta(0, 4));
    malha.adicionarAresta(Aresta(1, 5));
    malha.adicionarAresta(Aresta(2, 6));
    malha.adicionarAresta(Aresta(3, 7));

    // Faces trianguladas (12 triângulos = 6 faces * 2)
    // Base (normal -Y)
    malha.adicionarFace(Face(0, 2, 1));
    malha.adicionarFace(Face(0, 3, 2));

    // Topo (normal +Y)
    malha.adicionarFace(Face(4, 5, 6));
    malha.adicionarFace(Face(4, 6, 7));

    // Lado -X
    malha.adicionarFace(Face(0, 1, 5));
    malha.adicionarFace(Face(0, 5, 4));

    // Lado +X
    malha.adicionarFace(Face(3, 6, 2));
    malha.adicionarFace(Face(3, 7, 6));

    // Lado -Z
    malha.adicionarFace(Face(0, 7, 3));
    malha.adicionarFace(Face(0, 4, 7));

    // Lado +Z
    malha.adicionarFace(Face(1, 2, 6));
    malha.adicionarFace(Face(1, 6, 5));

    return malha;
}
