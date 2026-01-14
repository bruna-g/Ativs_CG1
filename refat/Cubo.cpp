#include "Cubo.h"
#include "Utils.h"
#include "Plano.h"
#include <cmath>

struct Face {
    Point p;
    Vector n;
};

Cubo::Cubo(const Point& centro_c, const float aresta_c)
    : centro_base(centro_c), aresta(aresta_c) {
}

IntersecaoCubo Cubo::CalcularIntersecaoCompleta(const Point& origem, const Vector& dir) const {
    IntersecaoCubo resultado = { false, INFINITY, Vector(0, 0, 0), Point(0, 0, 0) };

    float half = aresta / 2.0f;

    Face faces[6] = {
        {Point(centro_base.x, centro_base.y, centro_base.z - half), Vector(0, 0, -1)},
        {Point(centro_base.x, centro_base.y, centro_base.z + half), Vector(0, 0, 1)},
        {Point(centro_base.x - half, centro_base.y, centro_base.z), Vector(-1, 0, 0)},
        {Point(centro_base.x + half, centro_base.y, centro_base.z), Vector(1, 0, 0)},
        {Point(centro_base.x, centro_base.y, centro_base.z), Vector(0, -1, 0)},
        {Point(centro_base.x, centro_base.y + aresta, centro_base.z), Vector(0, 1, 0)}
    };

    for (int i = 0; i < 6; i++) {
        float denom = calcula_prod_esc(dir, faces[i].n);
        if (fabs(denom) < 1e-6f) continue;

        Vector w = subtrai_pontos(origem, faces[i].p);
        float t_aux = calcula_t_plano(w, faces[i].n, dir);
        if (t_aux <= 1e-4f || t_aux >= resultado.t) continue;

        Point Pi = calcula_eq_ray(origem, t_aux, dir);

        bool dentro_x = (Pi.x >= centro_base.x - half - 1e-4f) &&
            (Pi.x <= centro_base.x + half + 1e-4f);
        bool dentro_y = (Pi.y >= centro_base.y - 1e-4f) &&
            (Pi.y <= centro_base.y + aresta + 1e-4f);
        bool dentro_z = (Pi.z >= centro_base.z - half - 1e-4f) &&
            (Pi.z <= centro_base.z + half + 1e-4f);

        if (dentro_x && dentro_y && dentro_z) {
            resultado.intercepta = true;
            resultado.t = t_aux;
            resultado.normal = faces[i].n;
            resultado.ponto = Pi;
        }
    }

    return resultado;
}

IntersecaoCubo calcula_intersecao_cubo_completa(const Cubo& cubo, const Vector& dr, Point& Po) {
    return cubo.CalcularIntersecaoCompleta(Po, dr);
}
