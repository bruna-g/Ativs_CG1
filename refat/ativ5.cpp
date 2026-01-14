#include <iostream>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <SDL2/SDL.h>
#include "Textura.hpp"
#include "Color.h"
#include "Plano.h"
#include "Cone.h"
#include "Cubo.h"
#include "Point.h"
#include "Vector.h"
#include "Cilindro.h"
#include "Utils.h"
#include "Matriz3x3.h"
using namespace std;


float wJanela = 60.0f;
float hJanela = 60.0f;
int nCol = 500;
int nLin = 500;

float dJanela = 30.0f;
float z = -dJanela;

float Dx = wJanela / nCol;
float Dy = hJanela / nLin;

Point Po(0.f, 0.f, 0.f);

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
Color I_F(0.7f, 0.7f, 0.7f); // Intensidade da luz (branca)
Point P_F(-100.f, 140.f, -20.f);    // Posição da fonte de luz

//luz ambiente
Color I_A(0.3f, 0.3f, 0.3f); // Intensidade da luz ambiente

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

struct Face {
    Point p;
    Vector n;
};


IntersecaoCubo calcula_intersecao_cubo_completa(const Cubo& cubo, const Vector& dr, Point& Po);



Matriz3x3 M_id(1.0f, 0.0f, 0.0f,
    0.0f, 1.0f, 0.0f,
    0.0f, 0.0f, 1.0f);

Matriz3x3 dc_transposto = transpostaVetor(dc);
Matriz3x3 M = matrizSubtrai(M_id, outerProduto(dc));
Matriz3x3 M_conjugada = outerProduto(dc);
Matriz3x3 M_estrela = matrizSubtrai(M_conjugada, escalarMatriz(pow((altura_cone / raio_cone), 2), M));

Color calcula_Plano(Plano P, Vector dr, Color K_e, Color K_d) {
    bool naSombraEsf = false;
    bool naSombraCil = false;
    bool naSombraCone = false;

    // descobrir ponto Pi no plano a partir do raio que sai do olho do observador
    float ti = P.CalcularIntersecao(Po, dr);
    Point Pi = calcula_eq_ray(Po, ti, dr);

    // Verifica sombra: lança um raio de Pi em direção à fonte de luz e checa interseção com a esfera e o cone
    Vector l = calcula_l(P_F, Pi);
    // origem levemente deslocada para evitar auto-interseção
    Point Pi_mod(Pi.x + l.i * 1e-4f, Pi.y + l.j * 1e-4f, Pi.z + l.k * 1e-4f);

    Color Kd_textured = K_d; // fallback: sem textura
    if (P.tem_textura) {
        Vector n = P.n;
        Vector arbitrary = (fabs(n.i) < 0.9f) ? Vector(1.0f, 0.0f, 0.0f) : Vector(0.0f, 1.0f, 0.0f);
        Vector u_axis = cross(arbitrary, n);
        float nu = calcula_norma(u_axis);
        if (nu == 0.0f) nu = 1.0f;
        u_axis = calcula_esc_por_vetor(1.0f / nu, u_axis);
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
        float t_cone = cone.CalcularIntersecao(Pi_mod, l);
        if (t_cone > 1e-4f && t_cone < calcula_norma(subtrai_pontos(P_F, Pi_mod))) naSombraCone = true;
    }

    Color Ia(I_A.r * Kd_textured.r, I_A.g * Kd_textured.g, I_A.b * Kd_textured.b);
    if (naSombraEsf || naSombraCil || naSombraCone) {
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


int main() {
    std::ofstream img("ativ5.ppm");
    img << "P3\n" << nCol << " " << nLin << "\n255\n";
    texturaMadeira = new Textura("madeira", "madeira.bmp");

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

            // Interseção com os planos
            float ti_c = plano_chao.CalcularIntersecao(Po, dr_e);
            float ti_f = plano_fundo.CalcularIntersecao(Po, dr_e);
            float ti_e = plano_esq.CalcularIntersecao(Po, dr_e);
            float ti_d = plano_dir.CalcularIntersecao(Po, dr_e);
            float ti_t = plano_teto.CalcularIntersecao(Po, dr_e);

            // Objetos
            float t_cil = cilindro.CalcularIntersecao(Po, dr_e);
            float ti_cone = cone.CalcularIntersecao(Po, Pj);

            IntersecaoCubo inter_cubo = cubo.CalcularIntersecaoCompleta(Po, dr_e);
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
        float t_cil_shadow = cilindro.CalcularIntersecao(Pi_mod, l);
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

