#include "Color.h"
#include "Plano.h"
#include "Vector.h"
#include "Cilindro.h"
#include "Point.h"
#include "Cone.h"
#include "Textura.hpp"
#include <cmath>
#include <algorithm>

// Forward declarations e variáveis externas necessárias
extern Point Po;
extern Point P_F;
extern Point centroEsfera;
extern float rEsfera;
extern Color I_F;
extern Color I_A;
extern float m_e;
extern Color KCil_a;
extern Color KCil_d;
extern Color KCil_e;
extern Color KCone;
extern Textura* texturaMadeira;
extern Cone cone;
extern Cilindro cilindro;
extern float altura_cone;
extern float raio_cone;
extern Vector dc;

// Funções auxiliares externas
extern Vector subtrai_pontos(Point& p1, Point& p2);
extern float calcula_norma(const Vector& v);
extern float calcula_prod_esc(const Vector& v1, const Vector& v2);
extern Vector calcula_l(Point PF, Point Pi);
extern Point calcula_eq_ray(Point Po, float t, Vector dr);
extern Vector calcula_esc_por_vetor(float v1, const Vector& v2);
extern float calcula_t_plano(Vector w, Vector n, Vector dr);
extern float calcula_t_cone(Point& Pj);
extern float calcula_t_Cilindro(Cilindro cilindro, Vector dr, Point Po);
extern Vector cross(const Vector& a, const Vector& b);
extern Vector calcula_n_cilindro(Point P, Cilindro cilindro);
extern Vector calcula_dr(Point& Po, Point& Pj);

class Matriz3x3 {
public:
    float m[3][3];
    Matriz3x3(float m00, float m01, float m02,
        float m10, float m11, float m12,
        float m20, float m21, float m22) {
        m[0][0] = m00; m[0][1] = m01; m[0][2] = m02;
        m[1][0] = m10; m[1][1] = m11; m[1][2] = m12;
        m[2][0] = m20; m[2][1] = m21; m[2][2] = m22;
    }
};

extern Matriz3x3 M_id;
extern Matriz3x3 outerProduto(const Vector& d);
extern Matriz3x3 matrizSubtrai(const Matriz3x3& A, const Matriz3x3& B);
extern Vector matrizVetorProduto(const Matriz3x3& M, const Vector& v);

static inline float lidarExcecao(float v) {
    if (v < 0.0f) return 0.0f;
    if (v > 1.0f) return 1.0f;
    return v;
}

Color::Color(float r_c, float g_c, float b_c) {
    this->r = r_c;
    this->g = g_c;
    this->b = b_c;
}

// Nota: As implementações das funções Color estão em ativ5.cpp devido a dependências globais
// Para uso completo, estas funções devem ser mantidas em ativ5.cpp ou refatoradas
