#include "Cubo.h"
#include "Utils.h"
#include "Plano.h"
#include <cmath>

struct Face {
    Point p;
    Vector n;
};

Cubo::Cubo(const Point& centro_c, const float aresta_c) 
    : centro_base(centro_c), aresta(aresta_c) {}

IntersecaoCubo calcula_intersecao_cubo_completa(const Cubo& cubo, const Vector& dr, Point& Po) {
    IntersecaoCubo resultado = {false, INFINITY, Vector(0, 0, 0), Point(0, 0, 0)};
    
    float half = cubo.aresta / 2.0f;
    
    // Define as 6 faces do cubo
    Face faces[6] = {
        // Face frontal (z = centro.z - half)
        {Point(cubo.centro_base.x, cubo.centro_base.y, cubo.centro_base.z - half), Vector(0, 0, -1)},
        // Face traseira (z = centro.z + half)
        {Point(cubo.centro_base.x, cubo.centro_base.y, cubo.centro_base.z + half), Vector(0, 0, 1)},
        // Face esquerda (x = centro.x - half)
        {Point(cubo.centro_base.x - half, cubo.centro_base.y, cubo.centro_base.z), Vector(-1, 0, 0)},
        // Face direita (x = centro.x + half)
        {Point(cubo.centro_base.x + half, cubo.centro_base.y, cubo.centro_base.z), Vector(1, 0, 0)},
        // Face inferior (y = centro.y)
        {Point(cubo.centro_base.x, cubo.centro_base.y, cubo.centro_base.z), Vector(0, -1, 0)},
        // Face superior (y = centro.y + aresta)
        {Point(cubo.centro_base.x, cubo.centro_base.y + cubo.aresta, cubo.centro_base.z), Vector(0, 1, 0)}
    };
    
    // Testa interseção com cada face
    for (int i = 0; i < 6; i++) {
        float denom = calcula_prod_esc(dr, faces[i].n);
        
        // Se o denominador for próximo de zero, o raio é paralelo ao plano
        if (fabs(denom) < 1e-6f) continue;
        
        Vector w = subtrai_pontos(Po, faces[i].p);
        float t_aux = calcula_t_plano(w, faces[i].n, dr);
        
        // Verifica se t é positivo e menor que o melhor t encontrado
        if (t_aux <= 1e-4f || t_aux >= resultado.t) continue;
        
        // Calcula o ponto de interseção
        Point Pi = calcula_eq_ray(Po, t_aux, dr);
        
        // Verifica se o ponto está dentro dos limites do cubo
        bool dentro_x = (Pi.x >= cubo.centro_base.x - half - 1e-4f) && 
                       (Pi.x <= cubo.centro_base.x + half + 1e-4f);
        bool dentro_y = (Pi.y >= cubo.centro_base.y - 1e-4f) && 
                       (Pi.y <= cubo.centro_base.y + cubo.aresta + 1e-4f);
        bool dentro_z = (Pi.z >= cubo.centro_base.z - half - 1e-4f) && 
                       (Pi.z <= cubo.centro_base.z + half + 1e-4f);
        
        if (dentro_x && dentro_y && dentro_z) {
            resultado.intercepta = true;
            resultado.t = t_aux;
            resultado.normal = faces[i].n;
            resultado.ponto = Pi;
        }
    }
    
    return resultado;
}
