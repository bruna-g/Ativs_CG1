#include <iostream>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <SDL2/SDL.h>
#include "../include/Textura.hpp"
#include "../include/Color.h"
#include "../include/Plano.h"
#include "../include/Cone.h"
#include "../include/Cubo.h"
#include "../include/Point.h"
#include "../include/Vector.h"
#include "../include/Cilindro.h"
#include "../include/Utils.h"
#include "../include/Matriz3x3.h"
#include "../include/Cena.h"
#include "../include/Material.h"
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

static Cena* gCena = nullptr;

//chão
Point P_pi_chao(0.f, -150.f, 0.f);
Vector n_chao(0.f, 1.0f, 0.f);

Color KC_d(0.2f, 0.7f, 0.2f);
Color KC_e(0.0f, 0.0f, 0.0f);
float m_c = 1.0f;

Material mat_chao;
Plano plano_chao(P_pi_chao, n_chao, mat_chao);

//fundo
Point P_pi_fundo(200.f, -150.f, -400.f);
Vector n_fundo(0.f, 0.f, 1.f);
Plano plano_fundo(P_pi_fundo, n_fundo);

Color KF_d(0.686f, 0.933f, 0.933f);
Color KF_e(0.686f, 0.933f, 0.933f);
float m_f = 1.0f;

Material mat_fundo;

//parede esquerda
Point P_pi_esq(-200.f, -150.f, 0.f);
Vector n_esq(1.f, 0.f, 0.f);
Plano plano_esq(P_pi_esq, n_esq);

Color KE_d(0.686f, 0.933f, 0.933f);
Color KE_e(0.686f, 0.933f, 0.933f);

Material mat_esq;

//parede direita
Point P_pi_dir(200.f, -150.f, 0.f);
Vector n_dir(-1.f, 0.f, 0.f);
Plano plano_dir(P_pi_dir, n_dir);

Color KD_d(0.686f, 0.933f, 0.933f);
Color KD_e(0.686f, 0.933f, 0.933f);

Material mat_dir;

//teto
Point P_pi_teto(0.f, 150.f, 0.f);
Vector n_teto(0.f, -1.f, 0.f);
Plano plano_teto(P_pi_teto, n_teto);

Color KT_d(0.933f, 0.933f, 0.933f);
Color KT_e(0.933f, 0.933f, 0.933f);

Material mat_teto;

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

Material mat_cil;

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

Material mat_cone;

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
    if (gCena == nullptr) return Color(0, 0, 0);
    (void)K_e;
    (void)K_d;
    return P.CalcularCor(*gCena, dr);

}



Color calcula_color_cil(Cilindro cilindro, float t_cil, Vector dr) {
    if (gCena == nullptr) return Color(0, 0, 0);
    return cilindro.CalcularCor(*gCena, t_cil, dr);

}


int main() {
    std::ofstream img("ativ5.ppm");
    img << "P3\n" << nCol << " " << nLin << "\n255\n";
    texturaMadeira = new Textura("madeira", "madeira.bmp");

    // Materiais (inicializa aqui para poder apontar para a textura carregada)
    mat_chao.Ka = KC_d;
    mat_chao.Kd = KC_d;
    mat_chao.Ke = KC_e;
    mat_chao.m = m_e;
    mat_chao.usarTextura = true;
    mat_chao.textura = texturaMadeira;

    mat_fundo.Ka = KF_d;
    mat_fundo.Kd = KF_d;
    mat_fundo.Ke = KF_e;
    mat_fundo.m = m_f;

    mat_esq.Ka = KE_d;
    mat_esq.Kd = KE_d;
    mat_esq.Ke = KE_e;
    mat_esq.m = m_e;

    mat_dir.Ka = KD_d;
    mat_dir.Kd = KD_d;
    mat_dir.Ke = KD_e;
    mat_dir.m = m_e;

    mat_teto.Ka = KT_d;
    mat_teto.Kd = KT_d;
    mat_teto.Ke = KT_e;
    mat_teto.m = m_e;

    mat_cil.Ka = KCil_a;
    mat_cil.Kd = KCil_d;
    mat_cil.Ke = KCil_e;
    mat_cil.m = m_e;

    mat_cone = Material::Solido(KCone, m_e);

    // Aplica os materiais nos objetos globais (já instanciados)
    plano_chao.material = mat_chao;
    plano_fundo.material = mat_fundo;
    plano_esq.material = mat_esq;
    plano_dir.material = mat_dir;
    plano_teto.material = mat_teto;
    cilindro.material = mat_cil;
    cone.material = mat_cone;

    Cena cena;
    cena.observador = Po;
    cena.luz = LuzPontual{ P_F, I_F };
    cena.luzAmbiente = I_A;
    cena.centroEsfera = centroEsfera;
    cena.raioEsfera = rEsfera;
    cena.cone = &cone;
    cena.cilindro = &cilindro;
    cena.texturaMadeira = texturaMadeira;
    cena.expoenteEspecular = m_e;
    gCena = &cena;

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

            Color cor(100, 100, 100);

            //se interceptar fundo primeiro
            if (ti_f == t) {
                Color cor_plano = plano_fundo.CalcularCor(cena, dr_e);
                img << static_cast<int>(cor_plano.r) << ' ' << static_cast<int>(cor_plano.g) << ' ' << static_cast<int>(cor_plano.b) << ' ';
                continue;
            } //se interceptar o chão primeiro
            else if (ti_c == t) {
                Color cor_plano = plano_chao.CalcularCor(cena, dr_e);
                img << static_cast<int>(cor_plano.r) << ' ' << static_cast<int>(cor_plano.g) << ' ' << static_cast<int>(cor_plano.b) << ' ';
                continue;
            }
            else if (ti_e == t) {
                Color cor_plano = plano_esq.CalcularCor(cena, dr_e);
                img << static_cast<int>(cor_plano.r) << ' ' << static_cast<int>(cor_plano.g) << ' ' << static_cast<int>(cor_plano.b) << ' ';
                continue;
            }
            else if (ti_d == t) {
                Color cor_plano = plano_dir.CalcularCor(cena, dr_e);
                img << static_cast<int>(cor_plano.r) << ' ' << static_cast<int>(cor_plano.g) << ' ' << static_cast<int>(cor_plano.b) << ' ';
                continue;
            }
            else if (ti_t == t) {
                Color cor_plano = plano_teto.CalcularCor(cena, dr_e);
                img << static_cast<int>(cor_plano.r) << ' ' << static_cast<int>(cor_plano.g) << ' ' << static_cast<int>(cor_plano.b) << ' ';
                continue;
            }
            else if (t_cil == t) {
                Color color_Cil = cilindro.CalcularCor(cena, t, dr_e);
                img << static_cast<int>(color_Cil.r) << ' ' << static_cast<int>(color_Cil.g) << ' ' << static_cast<int>(color_Cil.b) << ' ';
                continue;
            }
            else if (ti_cone == t) {
                Color cor_cone = cone.CalcularCor(cena, t, dr_e);
                img << static_cast<int>(cor_cone.r) << ' ' << static_cast<int>(cor_cone.g) << ' ' << static_cast<int>(cor_cone.b) << ' ';
                continue;
            }
            else if (t_cubo == t) {
                // Renderiza o cubo (você pode criar uma função similar a calcula_Plano)
                Color K_cubo(1.f, 0.078f, 0.576f); // Cor vermelha para o cubo
                Plano plano_cubo(inter_cubo.ponto, inter_cubo.normal, Material::Solido(K_cubo, m_e));
                Color cor_cubo = plano_cubo.CalcularCor(cena, dr_e);
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
    if (gCena == nullptr) return Color(0, 0, 0);
    Vector dr = calcula_dr(Po, p);
    return cone.CalcularCor(*gCena, t, dr);
}

