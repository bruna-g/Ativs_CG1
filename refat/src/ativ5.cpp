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
#include "../include/Malha.h"
#include "../include/Esfera.h"
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

Esfera esfera(centroEsfera, rEsfera);


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
Malha cuboMalha;

Matriz3x3 M_id(1.0f, 0.0f, 0.0f,
    0.0f, 1.0f, 0.0f,
    0.0f, 0.0f, 1.0f);

Matriz3x3 dc_transposto = transpostaVetor(dc);
Matriz3x3 M = matrizSubtrai(M_id, outerProduto(dc));
Matriz3x3 M_conjugada = outerProduto(dc);
Matriz3x3 M_estrela = matrizSubtrai(M_conjugada, escalarMatriz(pow((altura_cone / raio_cone), 2), M));

Color calcula_Plano(Plano P, Vector dr, Color K_e, Color K_d) {
    if (gCena == nullptr) return Color(0, 0, 0);
    P.setKa(Vetor(K_d.r, K_d.g, K_d.b));
    P.setKd(Vetor(K_d.r, K_d.g, K_d.b));
    P.setKe(Vetor(K_e.r, K_e.g, K_e.b));
    P.setShininess(gCena->expoenteEspecular);
    return P.CalcularCor(*gCena, dr);

}



Color calcula_color_cil(Cilindro cilindro, float t_cil, Vector dr) {
    if (gCena == nullptr) return Color(0, 0, 0);
    cilindro.setKa(Vetor(KCil_a.r, KCil_a.g, KCil_a.b));
    cilindro.setKd(Vetor(KCil_d.r, KCil_d.g, KCil_d.b));
    cilindro.setKe(Vetor(KCil_e.r, KCil_e.g, KCil_e.b));
    cilindro.setShininess(gCena->expoenteEspecular);
    return cilindro.CalcularCor(*gCena, t_cil, dr);

}


int main() {
    std::ofstream img("ativ5.ppm");
    img << "P3\n" << nCol << " " << nLin << "\n255\n";
    texturaMadeira = new Textura("madeira", "madeira.bmp");

    // Material da esfera usando o Objeto base
    esfera.setKa(Vetor(K_e.r, K_e.g, K_e.b));
    esfera.setKd(Vetor(K_e.r, K_e.g, K_e.b));
    esfera.setKe(Vetor(K_e.r, K_e.g, K_e.b));
    esfera.setShininess(m_e);

    // Cubo como malha
    Color K_cubo(1.f, 0.078f, 0.576f);
    cuboMalha = Cubo::criarCubo(
        Vetor(0.f, -150.f, -165.f, 0.0f),
        40.0,
        Vetor(K_cubo.r, K_cubo.g, K_cubo.b, 0.0f),
        Vetor(K_cubo.r, K_cubo.g, K_cubo.b, 0.0f),
        Vetor(K_cubo.r, K_cubo.g, K_cubo.b, 0.0f),
        m_e
    );

    // Materiais (inicializa aqui para poder apontar para a textura carregada)
    mat_chao.usarTextura = true;
    mat_chao.textura = texturaMadeira;

    plano_chao.material = mat_chao;

    // Materiais via Objeto (Ka/Kd/Ke/shininess)
    plano_chao.setKa(Vetor(KC_d.r, KC_d.g, KC_d.b));
    plano_chao.setKd(Vetor(KC_d.r, KC_d.g, KC_d.b));
    plano_chao.setKe(Vetor(KC_e.r, KC_e.g, KC_e.b));
    plano_chao.setShininess(m_e);

    plano_fundo.setKa(Vetor(KF_d.r, KF_d.g, KF_d.b));
    plano_fundo.setKd(Vetor(KF_d.r, KF_d.g, KF_d.b));
    plano_fundo.setKe(Vetor(KF_e.r, KF_e.g, KF_e.b));
    plano_fundo.setShininess(m_e);

    plano_esq.setKa(Vetor(KE_d.r, KE_d.g, KE_d.b));
    plano_esq.setKd(Vetor(KE_d.r, KE_d.g, KE_d.b));
    plano_esq.setKe(Vetor(KE_e.r, KE_e.g, KE_e.b));
    plano_esq.setShininess(m_e);

    plano_dir.setKa(Vetor(KD_d.r, KD_d.g, KD_d.b));
    plano_dir.setKd(Vetor(KD_d.r, KD_d.g, KD_d.b));
    plano_dir.setKe(Vetor(KD_e.r, KD_e.g, KD_e.b));
    plano_dir.setShininess(m_e);

    plano_teto.setKa(Vetor(KT_d.r, KT_d.g, KT_d.b));
    plano_teto.setKd(Vetor(KT_d.r, KT_d.g, KT_d.b));
    plano_teto.setKe(Vetor(KT_e.r, KT_e.g, KT_e.b));
    plano_teto.setShininess(m_e);

    cilindro.setKa(Vetor(KCil_a.r, KCil_a.g, KCil_a.b));
    cilindro.setKd(Vetor(KCil_d.r, KCil_d.g, KCil_d.b));
    cilindro.setKe(Vetor(KCil_e.r, KCil_e.g, KCil_e.b));
    cilindro.setShininess(m_e);

    cone.setKa(Vetor(KCone.r, KCone.g, KCone.b));
    cone.setKd(Vetor(KCone.r, KCone.g, KCone.b));
    cone.setKe(Vetor(KCone.r, KCone.g, KCone.b));
    cone.setShininess(m_e);

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

            float t = -1.0f;
            if (esfera.verificarIntersecao(Vetor(Po.x, Po.y, Po.z), Vetor(dr_e.i, dr_e.j, dr_e.k))) {
                t = static_cast<float>(esfera.getDistancia());
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

            float t_cubo = -1.0f;
            if (cuboMalha.verificarIntersecao(Vetor(Po.x, Po.y, Po.z, 0.0f), Vetor(dr_e.i, dr_e.j, dr_e.k, 0.0f))) {
                t_cubo = static_cast<float>(cuboMalha.getDistancia());
            }

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
                Vetor PiV = cuboMalha.getPontoIntersecao();
                Vetor nV = cuboMalha.calcularNormal(PiV);
                Plano plano_cubo(Point(PiV.i, PiV.j, PiV.k), Vector(nV.i, nV.j, nV.k), Material());

                Vetor ka = cuboMalha.getKa();
                Vetor kd = cuboMalha.getKd();
                Vetor ke = cuboMalha.getKe();
                plano_cubo.setKa(ka);
                plano_cubo.setKd(kd);
                plano_cubo.setKe(ke);
                plano_cubo.setShininess(cuboMalha.getShininess());

                Color cor_cubo = plano_cubo.CalcularCor(cena, dr_e);
                img << static_cast<int>(cor_cubo.r) << ' ' << static_cast<int>(cor_cubo.g) << ' ' << static_cast<int>(cor_cubo.b) << ' ';
                continue;
            }
            else if (t <= 0.0f) {
                img << "100 100 100 ";
                continue;
            }

            // esfera
            Point Pi = calcula_eq_ray(Po, t, dr_e);
            Vetor nV = esfera.calcularNormal(Vetor(Pi.x, Pi.y, Pi.z));
            Vector n(nV.i, nV.j, nV.k);
            Vector X(-dr_e.i, -dr_e.j, -dr_e.k);
            Vector l = calcula_l(P_F, Pi);
            Vector r1 = calcula_esc_por_vetor(2.0f, calcula_esc_por_vetor(calcula_prod_esc(n, l), n));
            Vector r(r1.i - l.i, r1.j - l.j, r1.k - l.k);

            float fd = lidarExcecao(calcula_prod_esc(n, l));
            float cosAlpha = lidarExcecao(calcula_prod_esc(r, X));
            float fe = pow(cosAlpha, esfera.getShininess());

            Vetor keV = esfera.getKe();
            Color K_esf(keV.i, keV.j, keV.k);

            Color Id(I_F.r * K_esf.r * fd, I_F.g * K_esf.g * fd, I_F.b * K_esf.b * fd);
            Color Ie(I_F.r * K_esf.r * fe, I_F.g * K_esf.g * fe, I_F.b * K_esf.b * fe);
            Color Ia(I_A.r * K_esf.r, I_A.g * K_esf.g, I_A.b * K_esf.b);
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
    cone.setKa(Vetor(KCone.r, KCone.g, KCone.b));
    cone.setKd(Vetor(KCone.r, KCone.g, KCone.b));
    cone.setKe(Vetor(KCone.r, KCone.g, KCone.b));
    cone.setShininess(gCena->expoenteEspecular);
    return cone.CalcularCor(*gCena, t, dr);
}

