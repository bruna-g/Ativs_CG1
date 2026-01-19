#ifndef CENA_H
#define CENA_H

#include "Point.h"
#include "Color.h"

class Textura;
class Cone;
class Cilindro;
class Malha;

struct LuzPontual {
    Point pos;
    Color intensidade;
};

// Contexto mínimo para cálculo de cor/sombra.
// A ideia é evitar dependências de variáveis globais dentro dos objetos.
struct Cena {
    Point observador;
    LuzPontual luz;
    Color luzAmbiente;

    // Esfera (usada hoje para teste de sombra)
    Point centroEsfera;
    float raioEsfera;

    // Objetos que participam dos testes de sombra
    const Cone* cone = nullptr;
    const Cilindro* cilindro = nullptr;
    Malha* cubo = nullptr;

    // Recursos
    Textura* texturaMadeira = nullptr;

    // Parâmetros de shading
    float expoenteEspecular = 10.0f;

    Cena()
        : observador(0.0f, 0.0f, 0.0f),
        luz{ Point(0.0f, 0.0f, 0.0f), Color(0.0f, 0.0f, 0.0f) },
        luzAmbiente(0.0f, 0.0f, 0.0f),
        centroEsfera(0.0f, 0.0f, 0.0f),
        raioEsfera(0.0f) {
    }
};

#endif // CENA_H
