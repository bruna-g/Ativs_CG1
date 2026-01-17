#include <iostream>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <SDL2/SDL.h>
#include "Textura.hpp"
using namespace std;


class Point {
public:
    float x, y, z, p;
    Point(float x_p, float y_p, float z_p) {
        this->x = x_p;
        this->y = y_p;
        this->z = z_p;
        this->p = 1;
    }
};

class Vector {
public:
    float i, j, k, q;
    Vector(float i_v, float j_v, float k_v) {
        this->i = i_v;
        this->j = j_v;
        this->k = k_v;
        this->q = 0;
    }
};

class Color {
public:
    float r, g, b;
    Color(float r_c, float g_c, float b_c) {
        this->r = r_c;
        this->g = g_c;
        this->b = b_c;
    }
};

class Plano {
public:
    Point p_pi;
    Vector n;
    bool tem_textura;
    Plano(const Point& p_pi_p, const Vector& n_v, bool tem_textura_p = false) : p_pi(p_pi_p), n(n_v), tem_textura(tem_textura_p) {}
};

class Cilindro {
public:
    Point cb;
    float raio;
    float altura;
    Vector dc;
    Cilindro(const Point& cb_c, const float raio_c, const float altura_c, const Vector& dc_c) : cb(cb_c), raio(raio_c), altura(altura_c), dc(dc_c) {}
};

class Cone {
public:
    Point cb;
    Point v;
    float raio;
    Cone(const Point& cb_c, const Point& v_c, const float raio_c) : cb(cb_c), v(v_c), raio(raio_c) {}
};

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

class Cubo {
public:
    Point centro_base;
    float aresta;
    Cubo(const Point& centro_c, const float aresta_c) : centro_base(centro_c), aresta(aresta_c) {}
};

Vector subtrai_pontos(Point& p1, Point& p2);
float calcula_prod_esc(const Vector& v1, const Vector& v2);
Vector calcula_dr(Point& Po, Point& Pj);
Color calculaCone(float t, Point p);

float wJanela = 60.0f;
float hJanela = 60.0f;
int nCol = 500;
int nLin = 500;

float dJanela = 10.0f;
float z = -dJanela;

float Dx = wJanela / nCol;
float Dy = hJanela / nLin;

Point Po(50.f, 0.f, 0.f);

//esfera
float rEsfera = 5.0f;
Point centroEsfera(0.f, 95.f, -200.f);
Color K_e(0.854f, 0.647f, 0.125f);  // Material com canal vermelho
float m_e = 10.0f;            // Expoente especular


//textura de madeira
Textura* texturaMadeira = nullptr;

//chão
Point P_pi_chao(0.f, -150.f, 0.f);
Vector n_chao(0.f, 1.0f, 0.f);
Plano plano_chao(P_pi_chao, n_chao, true);

Color KC_d(0.2f, 0.7f, 0.2f);
Color KC_e(0.0f, 0.0f, 0.0f);
float m_c = 1.0f;

//fundo
Point P_pi_fundo(200.f, -150.f, -400.f);
Vector n_fundo(0.f, 0.f, 1.f);
Plano plano_fundo(P_pi_fundo, n_fundo);

Color KF_d(0.686f, 0.933f, 0.933f);
Color KF_e(0.686f, 0.933f, 0.933f);
float m_f = 1.0f;

//parede esquerda
Point P_pi_esq(-200.f, -150.f, 0.f);
Vector n_esq(1.f, 0.f, 0.f);
Plano plano_esq(P_pi_esq, n_esq);

Color KE_d(0.686f, 0.933f, 0.933f);
Color KE_e(0.686f, 0.933f, 0.933f);

//parede direita
Point P_pi_dir(200.f, -150.f, 0.f);
Vector n_dir(-1.f, 0.f, 0.f);
Plano plano_dir(P_pi_dir, n_dir);

Color KD_d(0.686f, 0.933f, 0.933f);
Color KD_e(0.686f, 0.933f, 0.933f);

//teto
Point P_pi_teto(0.f, 150.f, 0.f);
Vector n_teto(0.f, -1.f, 0.f);
Plano plano_teto(P_pi_teto, n_teto);

Color KT_d(0.933f, 0.933f, 0.933f);
Color KT_e(0.933f, 0.933f, 0.933f);

//fonte luminosa
Color I_F(1.5f, 1.5f, 1.5f); // Intensidade da luz (branca)
Point P_F(-100.f, 140.f, -20.f);    // Posição da fonte de luz

//luz ambiente
Color I_A(0.8f, 0.8f, 0.8f); // Intensidade da luz ambiente

//cilindro
Point centroCilindro(0.f, -150.f, -200.f);
float raio_cil = 5.f;
float altura_cil = 90.f;
Vector dc(0.f, 1.0f, 0.f);
Cilindro cilindro(centroCilindro, raio_cil, altura_cil, dc);

Color KCil_d(0.824f, 0.706f, 0.549f);
Color KCil_e(0.824f, 0.706f, 0.549f);
Color KCil_a(0.824f, 0.706f, 0.549f);

//cone
Point aux(altura_cil* dc.i, altura_cil* dc.j, altura_cil* dc.k);
Point c_topo_cilindro(cilindro.cb.x + aux.x, cilindro.cb.y + aux.y, cilindro.cb.z + aux.z);
Point cb_cone(0.f, -60.f, -200.f);
float raio_cone = 90.f;
float altura_cone = 150.f;
Point aux_v(altura_cone* dc.i, altura_cone* dc.j, altura_cone* dc.k);
Point v_cone(cb_cone.x + aux_v.x, cb_cone.y + aux_v.y, cb_cone.z + aux_v.z);
Cone cone(cb_cone, v_cone, raio_cone);
float calcula_t_cone(Point& Pj);
float teta = atan(raio_cone / altura_cone);

Color KCone(0.f, 1.f, 0.498f);

//cubo
Cubo cubo(Point(0.f, -150.f, -165.f), 40.f);



Vector cross(const Vector& a, const Vector& b) {
    return Vector(a.j * b.k - a.k * b.j,
                  a.k * b.i - a.i * b.k,
                  a.i * b.j - a.j * b.i);
}

Matriz3x3 M_id(1.0f, 0.0f, 0.0f,
    0.0f, 1.0f, 0.0f,
    0.0f, 0.0f, 1.0f);

Matriz3x3 outerProduto(const Vector& d) {
    Matriz3x3 M(0.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 0.0f);
    M.m[0][0] = d.i * d.i;
    M.m[0][1] = d.i * d.j;
    M.m[0][2] = d.i * d.k;

    M.m[1][0] = d.j * d.i;
    M.m[1][1] = d.j * d.j;
    M.m[1][2] = d.j * d.k;

    M.m[2][0] = d.k * d.i;
    M.m[2][1] = d.k * d.j;
    M.m[2][2] = d.k * d.k;
    return M;
}

Matriz3x3 matrizSubtrai(const Matriz3x3& A, const Matriz3x3& B) {
    Matriz3x3 R(0.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 0.0f);
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c)
            R.m[r][c] = A.m[r][c] - B.m[r][c];
    return R;
}

Matriz3x3 escalarMatriz(const float s, const Matriz3x3& M) {
    Matriz3x3 R(0.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 0.0f);
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c)
            R.m[r][c] = s * M.m[r][c];
    return R;
}

Vector matrizVetorProduto(const Matriz3x3& M, const Vector& v) {
    float i_res = M.m[0][0] * v.i + M.m[0][1] * v.j + M.m[0][2] * v.k;
    float j_res = M.m[1][0] * v.i + M.m[1][1] * v.j + M.m[1][2] * v.k;
    float k_res = M.m[2][0] * v.i + M.m[2][1] * v.j + M.m[2][2] * v.k;
    return Vector(i_res, j_res, k_res);
}

Matriz3x3 transpostaVetor(const Vector& v) {
    Matriz3x3 M(0.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 0.0f);
    M.m[0][0] = v.i;
    M.m[1][0] = v.j;
    M.m[2][0] = v.k;
    return M;
}

Vector matrizPorMatriz(const Matriz3x3& A, const Matriz3x3& B) {
    float i_res = A.m[0][0] * B.m[0][0] + A.m[0][1] * B.m[1][0] + A.m[0][2] * B.m[2][0];
    float j_res = A.m[1][0] * B.m[0][0] + A.m[1][1] * B.m[1][0] + A.m[1][2] * B.m[2][0];
    float k_res = A.m[2][0] * B.m[0][0] + A.m[2][1] * B.m[1][0] + A.m[2][2] * B.m[2][0];
    return Vector(i_res, j_res, k_res);
}

float transpostoPorVetor(const Matriz3x3& M, const Vector& v) {
    float res = M.m[0][0] * v.i + M.m[1][0] * v.j + M.m[2][0] * v.k;
    return res;
}

Matriz3x3 dc_transposto = transpostaVetor(dc);
Matriz3x3 M = matrizSubtrai(M_id, outerProduto(dc));
Matriz3x3 M_conjugada = outerProduto(dc);
Matriz3x3 M_estrela = matrizSubtrai(M_conjugada, escalarMatriz(pow((altura_cone / raio_cone), 2), M));
// Utilitário simples para limitar valores no intervalo [0,1]
static inline float lidarExcecao(float v) {
    if (v < 0.0f) return 0.0f;
    if (v > 1.0f) return 1.0f;
    return v;
}


Vector subtrai_pontos(Point& p1, Point& p2) {
    Vector sub(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
    return sub;
}

Vector subtrai_vetores(Vector& v1, Vector& v2) {
    Vector sub(v1.i - v2.i, v1.j - v2.j, v1.k - v2.k);
    return sub;
}

float calcula_norma(const Vector& v) {
    return sqrt(v.i * v.i + v.j * v.j + v.k * v.k);
}

Vector calcula_dr(Point& Po, Point& Pj) {
    Vector Dr = subtrai_pontos(Pj, Po);
    float drNorma = calcula_norma(Dr);
    Vector dr(Dr.i / drNorma, Dr.j / drNorma, Dr.k / drNorma);
    return dr;
}

float calcula_prod_esc(const Vector& v1, const Vector& v2) {
    return (v1.i * v2.i + v1.j * v2.j + v1.k * v2.k);
}

Vector calcula_esc_por_vetor(float v1, const Vector& v2) {
    return Vector(v1 * v2.i, v1 * v2.j, v1 * v2.k);
}

Vector calcula_n(Point Pi, Point C, float R) {
    Vector Pi_C = subtrai_pontos(Pi, C);
    Vector n((Pi_C.i / R), (Pi_C.j / R), (Pi_C.k / R));
    return n;
}

Vector calcula_l(Point PF, Point Pi) {
    Vector PF_Pi = subtrai_pontos(PF, Pi);
    float norma = calcula_norma(PF_Pi);
    Vector l(PF_Pi.i / norma, PF_Pi.j / norma, PF_Pi.k / norma);
    return l;
}

Point calcula_eq_ray(Point Po, float t, Vector dr) {
    Vector t_dr(dr.i * t, dr.j * t, dr.k * t);
    // Corrige a equação do raio: P(t) = Po + t * d_r
    Point Pi(Po.x + t_dr.i, Po.y + t_dr.j, Po.z + t_dr.k);
    return Pi;
}

Vector calcula_n_cilindro(Point P, Cilindro cilindro) {
    Vector P_B = subtrai_pontos(P, cilindro.cb);
    float P_B_u = calcula_prod_esc(P_B, cilindro.dc);
    Vector P_B_u_u = calcula_esc_por_vetor(P_B_u, cilindro.dc);
    Vector n = subtrai_vetores(P_B, P_B_u_u);
    float norma = calcula_norma(n);
    return Vector(n.i / norma, n.j / norma, n.k / norma);
}

float calcula_t_plano(Vector w, Vector n, Vector dr) {
    float result = calcula_prod_esc(calcula_esc_por_vetor(-1, w), n) / calcula_prod_esc(dr, n);
    return result;
}

float calcula_t_Cilindro(Cilindro cilindro, Vector dr, Point Po) {
    Vector Po_B = subtrai_pontos(Po, cilindro.cb);
    Vector Po_B_u_u = calcula_esc_por_vetor(calcula_prod_esc(Po_B, cilindro.dc), cilindro.dc);
    Vector v = subtrai_vetores(Po_B, Po_B_u_u);

    Vector d_u_u = calcula_esc_por_vetor(calcula_prod_esc(dr, cilindro.dc), cilindro.dc);
    Vector w = subtrai_vetores(dr, d_u_u);

    float a_delta = calcula_prod_esc(w, w);
    float b_delta = calcula_prod_esc(v, w);
    float c_delta = calcula_prod_esc(v, v) - cilindro.raio * cilindro.raio;

    float delta = b_delta * b_delta - a_delta * c_delta;
    if (delta >= 0) {
        float t1 = (-b_delta + sqrt(delta)) / (a_delta);
        float t2 = (-b_delta - sqrt(delta)) / (a_delta);

        Point P1 = calcula_eq_ray(Po, t1, dr);
        Point P2 = calcula_eq_ray(Po, t2, dr);

        float P1_B_u = calcula_prod_esc(subtrai_pontos(P1, cilindro.cb), cilindro.dc);
        float P2_B_u = calcula_prod_esc(subtrai_pontos(P2, cilindro.cb), cilindro.dc);


        float t = -1.0f;
        if ((P1_B_u >= 0 && P1_B_u <= cilindro.altura && t1 > 0) &&
            (P2_B_u >= 0 && P2_B_u <= cilindro.altura && t2 > 0))
            t = min(t1, t2);
        else if (P1_B_u >= 0 && P1_B_u <= cilindro.altura && t1 > 0)
            t = t1;
        else if (P2_B_u >= 0 && P2_B_u <= cilindro.altura && t2 > 0)
            t = t2;
        return t;
    }
    else {
        return INFINITY;
    }
}

Color calcula_Plano(Plano P, Vector dr, Color K_e, Color K_d) {
    bool naSombraEsf = false;
    bool naSombraCil = false;

    // descobrir ponto Pi no plano a partir do raio que sai do olho do observador
    Vector w = subtrai_pontos(Po, P.p_pi);
    float ti = calcula_t_plano(w, P.n, dr);
    Point Pi = calcula_eq_ray(Po, ti, dr);

    // Verifica sombra: lança um raio de Pi em direção à fonte de luz e checa interseção com a esfera e o cone
    Vector l = calcula_l(P_F, Pi);
    // origem levemente deslocada para evitar auto-interseção
    Point Pi_mod(Pi.x + l.i * 1e-4f, Pi.y + l.j * 1e-4f, Pi.z + l.k * 1e-4f);

    Color Kd_textured = K_d; // fallback: sem textura
    if (P.tem_textura) {
        Vector n = P.n;
        Vector arbitrary = (fabs(n.i) < 0.9f) ? Vector(1.0f,0.0f,0.0f) : Vector(0.0f,1.0f,0.0f);
        Vector u_axis = cross(arbitrary, n);
        float nu = calcula_norma(u_axis);
        if (nu == 0.0f) nu = 1.0f;
        u_axis = calcula_esc_por_vetor(1.0f/nu, u_axis);
        Vector v_axis = cross(n, u_axis);

        Vector vecPi = subtrai_pontos(Pi, P.p_pi);
        float texScale = 0.02f; // ajuste para controlar tamanho do padrão
        float u = calcula_prod_esc(vecPi, u_axis) * texScale;
        float v = calcula_prod_esc(vecPi, v_axis) * texScale;
        u = u - floor(u); if (u < 0) u += 1.0f;
        v = v - floor(v); if (v < 0) v += 1.0f;

        size_t tx = static_cast<size_t>(u * texturaMadeira->get_largura_pixels()) % texturaMadeira->get_largura_pixels();
        size_t ty = static_cast<size_t>(v * texturaMadeira->get_altura_pixels()) % texturaMadeira->get_altura_pixels();

        // Textura::get_cor_pixel(x,y) retorna pixels[linha][coluna], então passamos (ty, tx)
        rgb px = texturaMadeira->get_cor_pixel(ty, tx);
        Color texCol(px[0] / 255.0f, px[1] / 255.0f, px[2] / 255.0f);

        Kd_textured = Color(texCol.r, texCol.g, texCol.b);
    }

    float dist_Pi_Luz = calcula_norma(subtrai_pontos(P_F, Pi_mod));

    //interseção com a esfera
    Vector w_sombra = subtrai_pontos(Pi_mod, centroEsfera);
    float a_delta = calcula_prod_esc(l, l);
    float b_delta = 2.0f * calcula_prod_esc(l, w_sombra);
    float c_delta = calcula_prod_esc(w_sombra, w_sombra) - rEsfera * rEsfera;

    float delta = b_delta * b_delta - 4.0f * a_delta * c_delta;
    float s1 = INFINITY;
    float s2 = INFINITY;

    if (delta > 0.f) {
        s1 = (-b_delta - sqrt(delta)) / (2.0f * a_delta);
        s2 = (-b_delta + sqrt(delta)) / (2.0f * a_delta);
        // se houver interseção positiva antes da fonte (s entre 0 e distância até a luz), ponto está em sombra
    }
    if ((s1 > 1e-4f && s1 < dist_Pi_Luz) || (s2 > 1e-4f && s2 < dist_Pi_Luz)) naSombraEsf = true;

    // Verificar sombra do cone
    if (!naSombraEsf) {
        Point Po_original = Po;  // salva a posição original do observador
        Po = Pi_mod;  // temporariamente muda Po para o ponto de interseção
        float t_cone = calcula_t_cone(P_F);  // verifica interseção com o cone
        Po = Po_original;  // restaura Po para a posição original

        if (t_cone > 1e-4f && t_cone < calcula_norma(subtrai_pontos(P_F, Pi_mod))) {
            naSombraEsf = true;
        }
    }

    Color Ia(I_A.r * Kd_textured.r, I_A.g * Kd_textured.g, I_A.b * Kd_textured.b);
    if (naSombraEsf || naSombraCil) {
        // Apenas componente ambiente
        Color I(lidarExcecao(Ia.r), lidarExcecao(Ia.g), lidarExcecao(Ia.b));
        int R = static_cast<int>(roundf(I.r * 255.0f));
        int G = static_cast<int>(roundf(I.g * 255.0f));
        int B = static_cast<int>(roundf(I.b * 255.0f));
        return Color(R, G, B);
    }

    Vector v = Vector(-dr.i, -dr.j, -dr.k);
    Vector r1 = calcula_esc_por_vetor(2.0f, calcula_esc_por_vetor(calcula_prod_esc(P.n, l), P.n));
    Vector r(r1.i - l.i, r1.j - l.j, r1.k - l.k);

    float fd = lidarExcecao(calcula_prod_esc(P.n, l));
    float cosAlpha = lidarExcecao(calcula_prod_esc(r, v));
    float fe = pow(cosAlpha, m_e);
    Color Id(I_F.r * Kd_textured.r * fd, I_F.g * Kd_textured.g * fd, I_F.b * Kd_textured.b * fd);
    Color Ie(I_F.r * K_e.r * fe, I_F.g * K_e.g * fe, I_F.b * K_e.b * fe);
    Color I(lidarExcecao(Id.r + Ie.r + Ia.r), lidarExcecao(Id.g + Ie.g + Ia.g), lidarExcecao(Id.b + Ie.b + Ia.b));
    int R = static_cast<int>(roundf(I.r * 255.0f));
    int G = static_cast<int>(roundf(I.g * 255.0f));
    int B = static_cast<int>(roundf(I.b * 255.0f));
    return Color(R, G, B);

}



Color calcula_color_cil(Cilindro cilindro, float t_cil, Vector dr) {
    Point P = calcula_eq_ray(Po, t_cil, dr);

    bool naSombra = false;
    float dist_Pi_Luz = calcula_norma(subtrai_pontos(P_F, P));
    if (t_cil < dist_Pi_Luz) naSombra = true;

    Color Ia(I_A.r * KCil_a.r, I_A.g * KCil_a.g, I_A.b * KCil_a.b);

    if (naSombra) {
        // Apenas componente ambiente
        Color I(lidarExcecao(Ia.r), lidarExcecao(Ia.g), lidarExcecao(Ia.b));
        int R = static_cast<int>(roundf(I.r * 255.0f));
        int G = static_cast<int>(roundf(I.g * 255.0f));
        int B = static_cast<int>(roundf(I.b * 255.0f));
        return Color(R, G, B);
    }

    Vector n = calcula_n_cilindro(P, cilindro);
    Vector l = calcula_l(P_F, P);
    Vector v(-dr.i, -dr.j, -dr.k);
    Vector r1 = calcula_esc_por_vetor(2.0f, calcula_esc_por_vetor(calcula_prod_esc(n, l), n));
    Vector r(r1.i - l.i, r1.j - l.j, r1.k - l.k);

    float fd = lidarExcecao(calcula_prod_esc(n, l));
    float cosAlpha = lidarExcecao(calcula_prod_esc(r, v));
    float fe = pow(cosAlpha, m_e);

    Color Id(I_F.r * KCil_d.r * fd, I_F.g * KCil_d.g * fd, I_F.b * KCil_d.b * fd);
    Color Ie(I_F.r * KCil_e.r * fe, I_F.g * KCil_e.g * fe, I_F.b * KCil_e.b * fe);
    Color I(lidarExcecao(Id.r + Ie.r + Ia.r), lidarExcecao(Id.g + Ie.g + Ia.g), lidarExcecao(Id.b + Ie.b + Ia.b));

    int R = static_cast<int>(roundf(I.r * 255.0f));
    int G = static_cast<int>(roundf(I.g * 255.0f));
    int B = static_cast<int>(roundf(I.b * 255.0f));
    return Color(R, G, B);

}

struct Face {
        Point p;
        Vector n;
    };
struct IntersecaoCubo {
    bool intercepta;
    float t;
    Vector normal;
    Point ponto;
};


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


int main() {
    std::ofstream img("ativ5.ppm");
    img << "P3\n" << nCol << " " << nLin << "\n255\n";
    texturaMadeira = new Textura("madeira","madeira.bmp");

    for (int l = 0; l < nLin; l++) {
        float y = hJanela / 2 - Dy / 2 - l * Dy;
        for (int c = 0; c < nCol; c++) {
            float x = -wJanela / 2 + Dx / 2 + c * Dx;
            Point Pj(x, y, z);

            Vector dr_e = calcula_dr(Po, Pj);

            Vector w_e = subtrai_pontos(Po, centroEsfera);

            float a_delta = calcula_prod_esc(dr_e, dr_e);
            float b_delta = 2.0f * calcula_prod_esc(dr_e, w_e);
            float c_delta = calcula_prod_esc(w_e, w_e) - rEsfera * rEsfera;

            float delta = b_delta * b_delta - 4.0f * a_delta * c_delta;

            float t = -1.0f;
            if (delta >= 0) {
                float t1 = (-b_delta + sqrt(delta)) / (2.0f * a_delta);
                float t2 = (-b_delta - sqrt(delta)) / (2.0f * a_delta);

                if (t1 > 0.0f && t2 > 0.0f) t = min(t1, t2);
                else if (t1 > 0.0f) t = t1;
                else if (t2 > 0.0f) t = t2;
            }

            // Verifica interseção com os planos
            Vector w_p_c = subtrai_pontos(Po, plano_chao.p_pi); // w do chão
            float ti_c = calcula_t_plano(w_p_c, plano_chao.n, dr_e); // ti do chão

            Vector w_p_f = subtrai_pontos(Po, plano_fundo.p_pi); // w do plano de fundo
            float ti_f = calcula_t_plano(w_p_f, plano_fundo.n, dr_e); // ti do plano de fundo

            Vector w_p_e = subtrai_pontos(Po, plano_esq.p_pi); // w da parede esquerda
            float ti_e = calcula_t_plano(w_p_e, plano_esq.n, dr_e); // ti da parede esquerda

            Vector w_p_d = subtrai_pontos(Po, plano_dir.p_pi); // w da parede direita
            float ti_d = calcula_t_plano(w_p_d, plano_dir.n, dr_e); // ti da parede direita

            Vector w_p_t = subtrai_pontos(Po, plano_teto.p_pi); // w do teto
            float ti_t = calcula_t_plano(w_p_t, plano_teto.n, dr_e); // ti do teto

            //cilindro
            float t_cil = calcula_t_Cilindro(cilindro, dr_e, Po);
            float ti_cone = calcula_t_cone(Pj);

            IntersecaoCubo inter_cubo = calcula_intersecao_cubo_completa(cubo, dr_e, Po);
            float t_cubo = inter_cubo.intercepta ? inter_cubo.t : -1.0f;

            if (ti_c > 0.0f && (ti_c < t || t < 0.0f)) t = ti_c;
            if (ti_f > 0.0f && (ti_f < t || t < 0.0f)) t = ti_f;
            if (ti_e > 0.0f && (ti_e < t || t < 0.0f)) t = ti_e;
            if (ti_d > 0.0f && (ti_d < t || t < 0.0f)) t = ti_d;
            if (ti_t > 0.0f && (ti_t < t || t < 0.0f)) t = ti_t;
            if (t_cil > 0.0f && (t_cil < t || t < 0.0f)) t = t_cil;
            if (ti_cone > 0.0f && (ti_cone < t || t < 0.0f))  t = ti_cone;
            if (t_cubo > 0.0f && (t_cubo < t || t < 0.0f))  t = t_cubo;

            //se interceptar fundo primeiro
            if (ti_f == t) {
                Color cor_plano = calcula_Plano(plano_fundo, dr_e, KF_e, KF_d);
                img << static_cast<int>(cor_plano.r) << ' ' << static_cast<int>(cor_plano.g) << ' ' << static_cast<int>(cor_plano.b) << ' ';
                continue;
            } //se interceptar o chão primeiro
            else if (ti_c == t) {
                Color cor_plano = calcula_Plano(plano_chao, dr_e, KC_e, KC_d);
                img << static_cast<int>(cor_plano.r) << ' ' << static_cast<int>(cor_plano.g) << ' ' << static_cast<int>(cor_plano.b) << ' ';
                continue;
            }
            else if (ti_e == t) {
                Color cor_plano = calcula_Plano(plano_esq, dr_e, KE_e, KE_d);
                img << static_cast<int>(cor_plano.r) << ' ' << static_cast<int>(cor_plano.g) << ' ' << static_cast<int>(cor_plano.b) << ' ';
                continue;
            }
            else if (ti_d == t) {
                Color cor_plano = calcula_Plano(plano_dir, dr_e, KD_e, KD_d);
                img << static_cast<int>(cor_plano.r) << ' ' << static_cast<int>(cor_plano.g) << ' ' << static_cast<int>(cor_plano.b) << ' ';
                continue;
            }
            else if (ti_t == t) {
                Color cor_plano = calcula_Plano(plano_teto, dr_e, KT_e, KT_d);
                img << static_cast<int>(cor_plano.r) << ' ' << static_cast<int>(cor_plano.g) << ' ' << static_cast<int>(cor_plano.b) << ' ';
                continue;
            }
            else if (t_cil == t) {
                Color color_Cil = calcula_color_cil(cilindro, t, dr_e);
                img << static_cast<int>(color_Cil.r) << ' ' << static_cast<int>(color_Cil.g) << ' ' << static_cast<int>(color_Cil.b) << ' ';
                continue;
            }
            else if (ti_cone == t) {
                Color cor_cone = calculaCone(t, Pj);
                img << static_cast<int>(cor_cone.r) << ' ' << static_cast<int>(cor_cone.g) << ' ' << static_cast<int>(cor_cone.b) << ' ';
                continue;
            }
            else if (t_cubo == t) {
                // Renderiza o cubo (você pode criar uma função similar a calcula_Plano)
                Color K_cubo(1.f, 0.078f, 0.576f); // Cor vermelha para o cubo
                Plano plano_cubo(inter_cubo.ponto, inter_cubo.normal);
                Color cor_cubo = calcula_Plano(plano_cubo, dr_e, K_cubo, K_cubo);
                img << static_cast<int>(cor_cubo.r) << ' ' << static_cast<int>(cor_cubo.g) << ' ' << static_cast<int>(cor_cubo.b) << ' ';
                continue;
            }
            else if (t <= 0.0f || delta < 0.0f) {
                img << "100 100 100 ";
                continue;
            }

            // esfera
            Point Pi = calcula_eq_ray(Po, t, dr_e);
            Vector n = calcula_n(Pi, centroEsfera, rEsfera);
            Vector X(-dr_e.i, -dr_e.j, -dr_e.k);
            Vector l = calcula_l(P_F, Pi);
            Vector r1 = calcula_esc_por_vetor(2.0f, calcula_esc_por_vetor(calcula_prod_esc(n, l), n));
            Vector r(r1.i - l.i, r1.j - l.j, r1.k - l.k);

            float fd = lidarExcecao(calcula_prod_esc(n, l));
            float cosAlpha = lidarExcecao(calcula_prod_esc(r, X));
            float fe = pow(cosAlpha, m_e);

            Color Id(I_F.r * K_e.r * fd, I_F.g * K_e.g * fd, I_F.b * K_e.b * fd);
            Color Ie(I_F.r * K_e.r * fe, I_F.g * K_e.g * fe, I_F.b * K_e.b * fe);
            Color Ia(I_A.r * K_e.r, I_A.g * K_e.g, I_A.b * K_e.b);
            Color I(lidarExcecao(Id.r + Ie.r + Ia.r), lidarExcecao(Id.g + Ie.g + Ia.g), lidarExcecao(Id.b + Ie.b + Ia.b));

            int R = static_cast<int>(roundf(I.r * 255.0f));
            int G = static_cast<int>(roundf(I.g * 255.0f));
            int B = static_cast<int>(roundf(I.b * 255.0f));
            img << R << ' ' << G << ' ' << B << ' ';
        }
        img << "\n";
    }
    img.close();

    delete texturaMadeira;
    SDL_Quit();

    std::cout << "Imagem gerada: ativ5.ppm\n";
    return 0;
}

float calcula_t_cone(Point& Pj) {
    Vector w_cone = subtrai_pontos(Po, cone.cb);
    Matriz3x3 w_cone_transposto = transpostaVetor(w_cone);
    Vector dr = calcula_dr(Po, Pj);
    Matriz3x3 dr_transposto = transpostaVetor(dr);
    Vector v(cone.v.x - Po.x, cone.v.y - Po.y, cone.v.z - Po.z);

    // float a = calcula_prod_esc(matrizPorMatriz(dr_transposto, M_estrela), dr);
    // float b = 2.0f * calcula_prod_esc(matrizPorMatriz(w_cone_transposto, M_estrela), dr)
    //     - 2.0f * altura_cone * transpostoPorVetor(dr_transposto, dc);
    // float c = calcula_prod_esc(matrizPorMatriz(w_cone_transposto, M_estrela), w_cone)
    //     - 2.0f * altura_cone * transpostoPorVetor(w_cone_transposto, dc)
    //     + altura_cone * altura_cone;

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
                tEscolhido = min(tEscolhido, t);
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
                tEscolhido = min(tEscolhido, t_base);
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
            float t_menor = min(t1, t2);
            tEscolhido = min(tEscolhido, t_menor);
            temIntersecao = true;
        }
        if (t1_valido) {
            tEscolhido = min(tEscolhido, t1);
            temIntersecao = true;
        }
        if (t2_valido) {
            tEscolhido = min(tEscolhido, t2);
            temIntersecao = true;
        }
    }

    if (delta < 0.0f) return -1.0f;

    if (temIntersecao) return tEscolhido;

    return -1.0f;
}


Color calculaCone(float t, Point p) {
    bool naSombraEsf = false;
    bool naSombraCil = false;
    Vector dr = calcula_dr(Po, p);
    Point Pi = calcula_eq_ray(Po, t, dr);

    // Verificar se a interseção está na base do cone
    Vector dist_centro = subtrai_pontos(Pi, cone.cb);
    float altura_Pi = calcula_prod_esc(dist_centro, dc);  // projeção no eixo do cone
    bool na_base = fabs(altura_Pi) < 1e-3f;

    Vector n = na_base
        ? Vector(
            (calcula_prod_esc(dc, dr) > 0 ? -dc.i : dc.i),
            (calcula_prod_esc(dc, dr) > 0 ? -dc.j : dc.j),
            (calcula_prod_esc(dc, dr) > 0 ? -dc.k : dc.k))
        : [&]() {
        Vector V = subtrai_pontos(cone.v, Pi);
        float vNorma = calcula_norma(V);
        Vector s_conjugado(V.i / vNorma, V.j / vNorma, V.k / vNorma);
        Matriz3x3 M_e = matrizSubtrai(M_id, outerProduto(s_conjugado));
        Vector N = matrizVetorProduto(M_e, dc);
        float N_norma = calcula_norma(N);
        return Vector(N.i / N_norma, N.j / N_norma, N.k / N_norma);
        }();

    Vector l = calcula_l(P_F, Pi);
    Vector v(-dr.i, -dr.j, -dr.k);
    Vector r1 = calcula_esc_por_vetor(2.0f, calcula_esc_por_vetor(calcula_prod_esc(n, l), n));
    Vector r(r1.i - l.i, r1.j - l.j, r1.k - l.k);
    Color Ia(I_A.r * KCone.r, I_A.g * KCone.g, I_A.b * KCone.b);

    Point Pi_mod(Pi.x + l.i * 1e-4f, Pi.y + l.j * 1e-4f, Pi.z + l.k * 1e-4f);
    float dist_Pi_Luz = calcula_norma(subtrai_pontos(P_F, Pi_mod));

    Vector w_sombra = subtrai_pontos(Pi_mod, centroEsfera);
    float a_delta = calcula_prod_esc(l, l);
    float b_delta = 2.0f * calcula_prod_esc(l, w_sombra);
    float c_delta = calcula_prod_esc(w_sombra, w_sombra) - rEsfera * rEsfera;

    float delta = b_delta * b_delta - 4.0f * a_delta * c_delta;
    float s1 = INFINITY;
    float s2 = INFINITY;

    if (delta > 0.f) {
        s1 = (-b_delta - sqrt(delta)) / (2.0f * a_delta);
        s2 = (-b_delta + sqrt(delta)) / (2.0f * a_delta);
    }
    if ((s1 > 1e-4f && s1 < dist_Pi_Luz) || (s2 > 1e-4f && s2 < dist_Pi_Luz)) naSombraEsf = true;

    if (!naSombraEsf) {
        float t_cil_shadow = calcula_t_Cilindro(cilindro, l, Pi_mod);
        if (t_cil_shadow > 1e-4f && t_cil_shadow < dist_Pi_Luz) {
            naSombraCil = true;
        }
    }

    if (naSombraCil || naSombraEsf) {
        Color I(lidarExcecao(Ia.r), lidarExcecao(Ia.g), lidarExcecao(Ia.b));
        int R = static_cast<int>(roundf(I.r * 255.0f));
        int G = static_cast<int>(roundf(I.g * 255.0f));
        int B = static_cast<int>(roundf(I.b * 255.0f));
        return Color(R, G, B);
    }

    float fd = lidarExcecao(calcula_prod_esc(n, l));
    float cosAlpha = lidarExcecao(calcula_prod_esc(r, v));
    float fe = pow(cosAlpha, m_e);

    Color Id(I_F.r * KCone.r * fd, I_F.g * KCone.g * fd, I_F.b * KCone.b * fd);
    Color Ie(I_F.r * KCone.r * fe, I_F.g * KCone.g * fe, I_F.b * KCone.b * fe);
    Color I(lidarExcecao(Id.r + Ie.r + Ia.r), lidarExcecao(Id.g + Ie.g + Ia.g), lidarExcecao(Id.b + Ie.b + Ia.b));
    int R = static_cast<int>(roundf(I.r * 255.0f));
    int G = static_cast<int>(roundf(I.g * 255.0f));
    int B = static_cast<int>(roundf(I.b * 255.0f));
    return Color(R, G, B);
}