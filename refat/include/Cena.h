#ifndef CENA_H
#define CENA_H

#include "Point.h"
#include "Vector.h"
#include "Color.h"

#include <vector>

class Textura;
class Objeto;

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

    // Objetos que participam dos testes de sombra (ex.: cones/cilindros/esferas/malhas).
    // Obs: não inclui planos por padrão.
    std::vector<Objeto*> objetosSombra;

    // Recursos
    Textura* texturaMadeira = nullptr;

    // Parâmetros de shading
    float expoenteEspecular = 10.0f;

    // Retorna true se algum objeto em objetosSombra bloquear o raio (Pi_mod + s*l)
    // antes de chegar na luz. Use "ignorar" para evitar auto-sombra.
    bool estaNaSombra(const Point& Pi_mod, const Vector& l, float dist_Pi_Luz, Objeto* ignorar = nullptr) const;

    Cena()
        : observador(0.0f, 0.0f, 0.0f),
        luz{ Point(0.0f, 0.0f, 0.0f), Color(0.0f, 0.0f, 0.0f) },
        luzAmbiente(0.0f, 0.0f, 0.0f) {
    }
};

#endif // CENA_H
