#include "Cone.h"
#include "Vector.h"
#include "Matriz3x3.h"
#include "Utils.h"
#include <cmath>
#include <algorithm>

// Forward declarations e variáveis externas
extern Point Po;
extern Cone cone;
extern Vector dc;
extern float altura_cone;
extern float raio_cone;
extern float teta;



extern Matriz3x3 transpostaVetor(const Vector& v);

Cone::Cone(const Point& cb_c, const Point& v_c, const float raio_c) 
    : cb(cb_c), v(v_c), raio(raio_c) {}

float calcula_t_cone(Point& Pj) {
    Vector w_cone = subtrai_pontos(Po, cone.cb);
    Matriz3x3 w_cone_transposto = transpostaVetor(w_cone);
    Vector dr = calcula_dr(Po, Pj);
    Matriz3x3 dr_transposto = transpostaVetor(dr);
    Vector v(cone.v.x - Po.x, cone.v.y - Po.y, cone.v.z - Po.z);

    float a = pow((calcula_prod_esc(dr, dc)), 2) - pow(cos(teta), 2) * calcula_prod_esc(dr, dr);
    float b = pow(cos(teta), 2) * calcula_prod_esc(v, dr) - calcula_prod_esc(v, dc) * calcula_prod_esc(dr, dc);
    float c = pow(calcula_prod_esc(v, dc), 2) - pow(cos(teta), 2) * calcula_prod_esc(v, v);

    float dn = calcula_prod_esc(dc, dr);
    float tEscolhido = 1000000.0f;
    bool temIntersecao = false;

    if (a == 0 && b != 0) {
        float t = -c / (2.0f * b);

        if (t > 0.0001f) {
            Point Pi = calcula_eq_ray(Po, t, dr);
            Vector w_pi = subtrai_pontos(cone.v, Pi);
            float h = calcula_prod_esc(w_pi, dc);
            if (h >= 0 && h <= altura_cone) {
                tEscolhido = std::min(tEscolhido, t);
                temIntersecao = true;
            }
        }
    }

    if (dn != 0) {
        Vector v_base = w_cone;
        float t_base = -(calcula_prod_esc(v_base, dc) / dn);

        if (t_base > 0.0001f) {
            Point p_base = calcula_eq_ray(Po, t_base, dr);
            Vector dist_centro = subtrai_pontos(p_base, cone.cb);
            float dist_quadrada = calcula_prod_esc(dist_centro, dist_centro);

            if (dist_quadrada <= cone.raio * cone.raio) {
                tEscolhido = std::min(tEscolhido, t_base);
                temIntersecao = true;
            }
        }
    }

    float delta = b * b - a * c;

    if (delta > 0) {
        float t1 = (-b - sqrt(delta)) / (a);
        float t2 = (-b + sqrt(delta)) / (a);

        // Função auxiliar para verificar se t está nos limites do cone
        auto validarT = [&](float t) -> bool {
            if (t <= 0.0f) return false;

            Point Pi = calcula_eq_ray(Po, t, dr);
            Vector w_pi = subtrai_pontos(cone.v, Pi);
            float h = calcula_prod_esc(w_pi, dc);

            // Verifica se está entre a base (h=0) e o vértice (h=altura_cone)
            return (h >= 0.0f && h <= altura_cone);
            };

        // Testa t1 e t2 na ordem correta
        bool t1_valido = validarT(t1);
        bool t2_valido = validarT(t2);

        if (t1_valido && t2_valido) {
            float t_menor = std::min(t1, t2);
            tEscolhido = std::min(tEscolhido, t_menor);
            temIntersecao = true;
        }
        if (t1_valido) {
            tEscolhido = std::min(tEscolhido, t1);
            temIntersecao = true;
        }
        if (t2_valido) {
            tEscolhido = std::min(tEscolhido, t2);
            temIntersecao = true;
        }
    }

    if (delta < 0.0f) return -1.0f;

    if (temIntersecao) return tEscolhido;

    return -1.0f;
}
