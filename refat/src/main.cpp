#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <limits>
#include <SDL2/SDL.h>
#include "../include/Textura.hpp"
#include "../include/Color.h"
#include "../include/Esfera.h"
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
#include "../include/Camera.h"
using namespace std;

// Configuração da janela de visualização
int nCol = 800;
int nLin = 800;

Point eye(175.f, 180.f, 300.f);
Point at(175.f, 150.f, 200.f);
Vector up(0.f, 1.f, 0.f);
float distanciaFocal = 30.0f;

// Campo de visão - IMPORTANTE: manter proporção igual a nCol/nLin
float xMin = -50.0f;
float xMax = 50.0f;
float hJanelaDesejada = (xMax - xMin) * nLin / nCol;  // Manter proporção
float yMin = -hJanelaDesejada / 2.0f;
float yMax = hJanelaDesejada / 2.0f;

Camera camera(eye, at, up, distanciaFocal, xMin, xMax, yMin, yMax);

float wJanela = camera.getLarguraJanela();
float hJanela = camera.getAlturaJanela();
float Dx = wJanela / nCol;
float Dy = hJanela / nLin;

Point Po(0.f, 10.f, 0.f);

// Esfera (topoete da árvore)
float rEsfera = 5.0f;
Point centroEsfera(100.f, 150.f + rEsfera, 100.f);
Color K_e(0.8f, 0.6f, 0.2f);
float m_e = 10.0f;
Esfera esfera(centroEsfera, rEsfera);

// Textura de madeira
Textura* texturaMadeira = nullptr;
static Cena* gCena = nullptr;

// Chão
Point P_pi_chao(100.f, 0.f, 100.f);
Vector n_chao(0.f, 1.0f, 0.f);
Color KC_d(0.2f, 0.7f, 0.2f);
Color KC_e(0.0f, 0.0f, 0.0f);
float m_c = 1.0f;
Material mat_chao;
Plano plano_chao(P_pi_chao, n_chao, mat_chao);

// Fundo
Point P_pi_fundo(100.f, 100.f, 400.f);
Vector n_fundo(0.f, 0.f, -1.f);
Color KF_d(0.686f, 0.933f, 0.933f);
Color KF_e(0.686f, 0.933f, 0.933f);
float m_f = 1.0f;
Material mat_fundo;
Plano plano_fundo(P_pi_fundo, n_fundo, mat_fundo);

// Parede esquerda
Point P_pi_esq(0.f, 100.f, 100.f);
Vector n_esq(1.f, 0.f, 0.f);
Color KE_d(0.686f, 0.933f, 0.933f);
Color KE_e(0.686f, 0.933f, 0.933f);
Material mat_esq;
Plano plano_esq(P_pi_esq, n_esq, mat_esq);

// Parede direita
Point P_pi_dir(200.f, 100.f, 100.f);
Vector n_dir(-1.f, 0.f, 0.f);
Color KD_d(0.686f, 0.933f, 0.933f);
Color KD_e(0.686f, 0.933f, 0.933f);
Material mat_dir;
Plano plano_dir(P_pi_dir, n_dir, mat_dir);

// Teto
Point P_pi_teto(100.f, 300.f, 100.f);
Vector n_teto(0.f, -1.f, 0.f);
Color KT_d(0.933f, 0.933f, 0.933f);
Color KT_e(0.933f, 0.933f, 0.933f);
Material mat_teto;
Plano plano_teto(P_pi_teto, n_teto, mat_teto);

// Fonte luminosa
Color I_F(0.7f, 0.7f, 0.7f);
Point P_F(-20.f, 185.f, 350.f);

// Luz ambiente
Color I_A(0.3f, 0.3f, 0.3f);

//-----------------Arvore 1---------------------//
// Cilindro (tronco da árvore)
Point centroCilindro(25.f, 00.f, 150.f);
float raio_cil = 8.f;
float altura_cil = 60.f;
Vector dc(0.f, 1.0f, 0.f);
Cilindro cilindro(centroCilindro, raio_cil, altura_cil, dc);
Color KCil_d(0.545f, 0.270f, 0.075f);
Color KCil_e(0.545f, 0.270f, 0.075f);
Color KCil_a(0.545f, 0.270f, 0.075f);
Material mat_cil;

// Cone (copa da árvore)
Point cb_cone(25.f, 60.f, 150.f);
float raio_cone = 40.f;
float altura_cone = 50.f;
Vector aux_v_cone = calcula_esc_por_vetor(altura_cone, dc);
Point v_cone(cb_cone.x + aux_v_cone.i, cb_cone.y + aux_v_cone.j, cb_cone.z + aux_v_cone.k);
Cone cone(cb_cone, v_cone, raio_cone);
Color KCone(0.f, 0.6f, 0.2f);
Material mat_cone;
//----------------------------------------------//

//-----------------Arvore 2---------------------//
// Cilindro (tronco da árvore)
Point centroCilindro2(175.f, 00.f, 0.f);
float raio_cil2 = 8.f;
float altura_cil2 = 60.f;
Vector dc2(0.f, 1.0f, 0.f);
Cilindro cilindro2(centroCilindro2, raio_cil2, altura_cil2, dc2);
Color KCil_d2(0.545f, 0.270f, 0.075f);
Color KCil_e2(0.545f, 0.270f, 0.075f);
Color KCil_a2(0.545f, 0.270f, 0.075f);
Material mat_cil2;

// Cone (copa da árvore)
Point cb_cone2(175.f, 60.f, 0.f);
float raio_cone2 = 40.f;
float altura_cone2 = 50.f;
Vector aux_v_cone2 = calcula_esc_por_vetor(altura_cone2, dc2);
Point v_cone2(cb_cone2.x + aux_v_cone2.i, cb_cone2.y + aux_v_cone2.j, cb_cone2.z + aux_v_cone2.k);
Cone cone2(cb_cone2, v_cone2, raio_cone2);
Color KCone2(0.f, 0.6f, 0.2f);
Material mat_cone2;
//----------------------------------------------//

//-----------------Arvore 3---------------------//
// Cilindro (tronco da árvore)
Point centroCilindro3(325.f, 00.f, 150.f);
float raio_cil3 = 8.f;
float altura_cil3 = 60.f;
Vector dc3(0.f, 1.0f, 0.f);
Cilindro cilindro3(centroCilindro3, raio_cil3, altura_cil3, dc3);
Color KCil_d3(0.545f, 0.270f, 0.075f);
Color KCil_e3(0.545f, 0.270f, 0.075f);
Color KCil_a3(0.545f, 0.270f, 0.075f);
Material mat_cil3;

// Cone (copa da árvore)
Point cb_cone3(325.f, 60.f, 150.f);
float raio_cone3 = 40.f;
float altura_cone3 = 50.f;
Vector aux_v_cone3 = calcula_esc_por_vetor(altura_cone3, dc3);
Point v_cone3(cb_cone3.x + aux_v_cone3.i, cb_cone3.y + aux_v_cone3.j, cb_cone3.z + aux_v_cone3.k);
Cone cone3(cb_cone3, v_cone3, raio_cone3);
Color KCone3(0.f, 0.6f, 0.2f);
Material mat_cone3;
//----------------------------------------------//
//-----------------Vaca---------------------//
// Cilindro (tronco da árvore)
Point centroCilindro4(150.f, 00.f, 110.f);
float raio_cil4 = 2.f;
float altura_cil4 = 20.f;
Vector dc4(0.f, 1.0f, 0.f);
Cilindro cilindro4(centroCilindro4, raio_cil4, altura_cil4, dc4);
Color KCil_d4(0.545f, 0.270f, 0.075f);
Color KCil_e4(0.545f, 0.270f, 0.075f);
Color KCil_a4(0.545f, 0.270f, 0.075f);
Material mat_cil4;

Point centroCilindro5(150.f, 00.f, 140.f);
float raio_cil5 = 2.f;
float altura_cil5 = 20.f;
Vector dc5(0.f, 1.0f, 0.f);
Cilindro cilindro5(centroCilindro5, raio_cil5, altura_cil5, dc5);
Color KCil_d5(0.545f, 0.270f, 0.075f);
Color KCil_e5(0.545f, 0.270f, 0.075f);
Color KCil_a5(0.545f, 0.270f, 0.075f);
Material mat_cil5;

Point centroCilindro6(190.f, 00.f, 110.f);
float raio_cil6 = 2.f;
float altura_cil6 = 20.f;
Vector dc6(0.f, 1.0f, 0.f);
Cilindro cilindro6(centroCilindro6, raio_cil6, altura_cil6, dc6);
Color KCil_d6(0.545f, 0.270f, 0.075f);
Color KCil_e6(0.545f, 0.270f, 0.075f);
Color KCil_a6(0.545f, 0.270f, 0.075f);
Material mat_cil6;

Point centroCilindro7(190.f, 00.f, 140.f);
float raio_cil7 = 2.f;
float altura_cil7 = 20.f;
Vector dc7(0.f, 1.0f, 0.f);
Cilindro cilindro7(centroCilindro7, raio_cil7, altura_cil7, dc7);
Color KCil_d7(0.545f, 0.270f, 0.075f);
Color KCil_e7(0.545f, 0.270f, 0.075f);
Color KCil_a7(0.545f, 0.270f, 0.075f);
Material mat_cil7;
//--------------------------------------//

// nave
Point cb_nave(175.f, 200.f, 130.f);
float raio_nave = 80.f;
float altura_nave = 60.f;
Vector aux_v_nave = calcula_esc_por_vetor(altura_nave, dc);
Point v_nave(cb_nave.x + aux_v_nave.i, cb_nave.y + aux_v_nave.j, cb_nave.z + aux_v_nave.k);
Cone nave(cb_nave, v_nave, raio_nave);
Color Ka_nave(0.15f, 0.15f, 0.15f); // Ka - Ambiente
Color Kd_nave(0.65f, 0.68f, 0.70f); // Kd - Difuso
Color Ke_nave(0.0f, 0.0f, 0.0f); // Ke - Especular
Material mat_nave;

// Esfera1_nave
float rEsfera1_nave = 4.0f;
Point centroEsfera1_nave(150.f, 205.f + rEsfera1_nave, 190.f);
Color K_e1_nave(1.f, 0.f, 0.f);  // Vermelho
Color K_d1_nave(1.f, 0.f, 0.f);  // Vermelho
Color K_a1_nave(1.f, 0.1f, 0.f);  // Vermelho
float m_e1_nave = 10.0f;
Esfera esfera1_nave(centroEsfera1_nave, rEsfera1_nave);

// Esfera2_nave
float rEsfera2_nave = 4.0f;
Point centroEsfera2_nave(175.f, 205.f + rEsfera2_nave, 196.6f);
Color K_e2_nave(1.f, 0.f, 0.f);  // Vermelho
Color K_d2_nave(1.f, 0.f, 0.f);  // Vermelho
Color K_a2_nave(1.f, 0.1f, 0.f);  // Vermelho
float m_e2_nave = 10.0f;
Esfera esfera2_nave(centroEsfera2_nave, rEsfera2_nave);

// Esfera3_nave
float rEsfera3_nave = 4.0f;
Point centroEsfera3_nave(200.f, 205.f + rEsfera3_nave, 190.f);
Color K_e3_nave(1.f, 0.f, 0.f);  // Vermelho
Color K_d3_nave(1.f, 0.f, 0.f);  // Vermelho
Color K_a3_nave(1.f, 0.1f, 0.f);  // Vermelho
float m_e3_nave = 10.0f;
Esfera esfera3_nave(centroEsfera3_nave, rEsfera3_nave);

// Cubo (como malha)
Malha cuboMalha;

Matriz3x3 M_id(1.0f, 0.0f, 0.0f,
    0.0f, 1.0f, 0.0f,
    0.0f, 0.0f, 1.0f);

int main() {
    SDL_Window* window = nullptr;
    SDL_Renderer* renderer = nullptr;

    SDL_Init(SDL_INIT_VIDEO);
    SDL_CreateWindowAndRenderer(nCol, nLin, 0, &window, &renderer);
    SDL_SetWindowTitle(window, "Trabalho Final - CG");
    texturaMadeira = new Textura("madeira", "madeira.bmp");

    // Materiais / propriedades via Objeto (Ka/Kd/Ke/shininess)
    esfera.setKa(Vetor(K_e.r, K_e.g, K_e.b));
    esfera.setKd(Vetor(K_e.r, K_e.g, K_e.b));
    esfera.setKe(Vetor(K_e.r, K_e.g, K_e.b));
    esfera.setShininess(m_e);

    // esfera 1 da nave
    esfera1_nave.setKa(Vetor(K_a1_nave.r, K_a1_nave.g, K_a1_nave.b));
    esfera1_nave.setKd(Vetor(K_d1_nave.r, K_d1_nave.g, K_d1_nave.b));
    esfera1_nave.setKe(Vetor(K_e1_nave.r, K_e1_nave.g, K_e1_nave.b));
    esfera1_nave.setShininess(m_e1_nave);

    // esfera 2 da nave
    esfera2_nave.setKa(Vetor(K_a2_nave.r, K_a2_nave.g, K_a2_nave.b));
    esfera2_nave.setKd(Vetor(K_d2_nave.r, K_d2_nave.g, K_d2_nave.b));
    esfera2_nave.setKe(Vetor(K_e2_nave.r, K_e2_nave.g, K_e2_nave.b));
    esfera2_nave.setShininess(m_e2_nave);

    // esfera 3 da nave
    esfera3_nave.setKa(Vetor(K_a3_nave.r, K_a3_nave.g, K_a3_nave.b));
    esfera3_nave.setKd(Vetor(K_d3_nave.r, K_d3_nave.g, K_d3_nave.b));
    esfera3_nave.setKe(Vetor(K_e3_nave.r, K_e3_nave.g, K_e3_nave.b));
    esfera3_nave.setShininess(m_e3_nave);

    // Chão com textura
    mat_chao.usarTextura = true;
    mat_chao.textura = texturaMadeira;
    plano_chao.material = mat_chao;
    plano_chao.setKa(Vetor(KC_d.r, KC_d.g, KC_d.b));
    plano_chao.setKd(Vetor(KC_d.r, KC_d.g, KC_d.b));
    plano_chao.setKe(Vetor(KC_e.r, KC_e.g, KC_e.b));
    plano_chao.setShininess(m_e);

    cilindro.setKa(Vetor(KCil_a.r, KCil_a.g, KCil_a.b));
    cilindro.setKd(Vetor(KCil_d.r, KCil_d.g, KCil_d.b));
    cilindro.setKe(Vetor(KCil_e.r, KCil_e.g, KCil_e.b));
    cilindro.setShininess(5);

    cilindro2.setKa(Vetor(KCil_a2.r, KCil_a2.g, KCil_a2.b));
    cilindro2.setKd(Vetor(KCil_d2.r, KCil_d2.g, KCil_d2.b));
    cilindro2.setKe(Vetor(KCil_e2.r, KCil_e2.g, KCil_e2.b));
    cilindro2.setShininess(5);

    cilindro3.setKa(Vetor(KCil_a3.r, KCil_a3.g, KCil_a3.b));
    cilindro3.setKd(Vetor(KCil_d3.r, KCil_d3.g, KCil_d3.b));
    cilindro3.setKe(Vetor(KCil_e3.r, KCil_e3.g, KCil_e3.b));
    cilindro3.setShininess(5);

    cilindro4.setKa(Vetor(KCil_a3.r, KCil_a3.g, KCil_a3.b));
    cilindro4.setKd(Vetor(KCil_d3.r, KCil_d3.g, KCil_d3.b));
    cilindro4.setKe(Vetor(KCil_e3.r, KCil_e3.g, KCil_e3.b));
    cilindro4.setShininess(m_e);

    cilindro5.setKa(Vetor(KCil_a3.r, KCil_a3.g, KCil_a3.b));
    cilindro5.setKd(Vetor(KCil_d3.r, KCil_d3.g, KCil_d3.b));
    cilindro5.setKe(Vetor(KCil_e3.r, KCil_e3.g, KCil_e3.b));
    cilindro5.setShininess(m_e);

    cilindro6.setKa(Vetor(KCil_a3.r, KCil_a3.g, KCil_a3.b));
    cilindro6.setKd(Vetor(KCil_d3.r, KCil_d3.g, KCil_d3.b));
    cilindro6.setKe(Vetor(KCil_e3.r, KCil_e3.g, KCil_e3.b));
    cilindro6.setShininess(m_e);

    cilindro7.setKa(Vetor(KCil_a3.r, KCil_a3.g, KCil_a3.b));
    cilindro7.setKd(Vetor(KCil_d3.r, KCil_d3.g, KCil_d3.b));
    cilindro7.setKe(Vetor(KCil_e3.r, KCil_e3.g, KCil_e3.b));
    cilindro7.setShininess(m_e);

    cone.setKa(Vetor(KCone.r, KCone.g, KCone.b));
    cone.setKd(Vetor(KCone.r, KCone.g, KCone.b));
    cone.setKe(Vetor(KCone.r, KCone.g, KCone.b));
    cone.setShininess(m_e);

    cone2.setKa(Vetor(KCone2.r, KCone2.g, KCone2.b));
    cone2.setKd(Vetor(KCone2.r, KCone2.g, KCone2.b));
    cone2.setKe(Vetor(KCone2.r, KCone2.g, KCone2.b));
    cone2.setShininess(m_e);

    cone3.setKa(Vetor(KCone3.r, KCone3.g, KCone3.b));
    cone3.setKd(Vetor(KCone3.r, KCone3.g, KCone3.b));
    cone3.setKe(Vetor(KCone3.r, KCone3.g, KCone3.b));
    cone3.setShininess(m_e);

    nave.setKa(Vetor(Ka_nave.r, Ka_nave.g, Ka_nave.b));
    nave.setKd(Vetor(Kd_nave.r, Kd_nave.g, Kd_nave.b));
    nave.setKe(Vetor(Ke_nave.r, Ke_nave.g, Ke_nave.b));
    nave.setShininess(80.0f);

    // Cubo como malha
    Color K_cubo(1.f, 0.078f, 0.576f);
    cuboMalha = Cubo::criarCubo(
        Vetor(100.f, 0.f, 100.f, 0.0f),
        40.0,
        Vetor(K_cubo.r, K_cubo.g, K_cubo.b, 0.0f),
        Vetor(K_cubo.r, K_cubo.g, K_cubo.b, 0.0f),
        Vetor(K_cubo.r, K_cubo.g, K_cubo.b, 0.0f),
        m_e
    );

    Cena cena;
    cena.observador = camera.getEye();
    cena.luz = LuzPontual{ P_F, I_F };
    cena.luzAmbiente = I_A;
    cena.objetosSombra = { &cilindro, &cone, &cilindro2, &cone2, &cilindro3, &cone3, &nave,
        &esfera1_nave, &esfera2_nave, &esfera3_nave, &cilindro4, &cilindro5, &cilindro6, &cilindro7 };
    cena.texturaMadeira = texturaMadeira;
    cena.expoenteEspecular = m_e;
    gCena = &cena;

    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);

    // Ray casting
    for (int linha = 0; linha < nLin; linha++) {
        float y = camera.getYMax() - (linha + 0.5f) * Dy;  // CORRIGIDO
        for (int coluna = 0; coluna < nCol; coluna++) {
            float x = camera.getXMin() + (coluna + 0.5f) * Dx;  // TAMBÉM CORRIGIDO PARA CONSISTÊNCIA

            Vector dr_e = camera.gerarRaio(x, y);
            Point Po = camera.getEye();

            float t_best = std::numeric_limits<float>::infinity();
            enum class Hit {
                None,
                Fundo,
                Chao,
                Esq,
                Dir,
                Teto,
                Cilindro,
                Cilindro2,
                Cilindro3,
                Cilindro4,
                Cilindro5,
                Cilindro6,
                Cilindro7,
                Cone,
                Cone2,
                Cone3,
                Cubo,
                Esfera,
                Nave,
                Esfera1_nave,
                Esfera2_nave,
                Esfera3_nave
            };
            Hit hit = Hit::None;

            auto considerar = [&](float t, Hit h) {
                if (t > 1e-4f && std::isfinite(t) && t < t_best) {
                    t_best = t;
                    hit = h;
                }
                };

            // Interseção com os planos
            float ti_c = plano_chao.CalcularIntersecao(Po, dr_e);

            // Objetos
            float t_cil = cilindro.CalcularIntersecao(Po, dr_e);
            float ti_cone = cone.CalcularIntersecao(Po, dr_e);

            float ti_esf1_nave = esfera1_nave.CalcularIntersecao(Po, dr_e);
            float ti_esf2_nave = esfera2_nave.CalcularIntersecao(Po, dr_e);
            float ti_esf3_nave = esfera3_nave.CalcularIntersecao(Po, dr_e);

            float t_cil2 = cilindro2.CalcularIntersecao(Po, dr_e);
            float ti_cone2 = cone2.CalcularIntersecao(Po, dr_e);

            float t_cil3 = cilindro3.CalcularIntersecao(Po, dr_e);
            float ti_cone3 = cone3.CalcularIntersecao(Po, dr_e);

            // float t_cubo = std::numeric_limits<float>::infinity();
            // if (cuboMalha.verificarIntersecao(
            //     Vetor(Po.x, Po.y, Po.z, 0.0f),
            //     Vetor(dr_e.i, dr_e.j, dr_e.k, 0.0f))) {
            //     t_cubo = static_cast<float>(cuboMalha.getDistancia());
            // }

            float ti_nave = nave.CalcularIntersecao(Po, dr_e);

            float t_cil4 = cilindro4.CalcularIntersecao(Po, dr_e);
            float t_cil5 = cilindro5.CalcularIntersecao(Po, dr_e);
            float t_cil6 = cilindro6.CalcularIntersecao(Po, dr_e);
            float t_cil7 = cilindro7.CalcularIntersecao(Po, dr_e);

            considerar(ti_c, Hit::Chao);
            considerar(t_cil, Hit::Cilindro);
            considerar(t_cil2, Hit::Cilindro2);
            considerar(t_cil3, Hit::Cilindro3);
            considerar(t_cil4, Hit::Cilindro4);
            considerar(t_cil5, Hit::Cilindro5);
            considerar(t_cil6, Hit::Cilindro6);
            considerar(t_cil7, Hit::Cilindro7);
            considerar(ti_cone, Hit::Cone);
            considerar(ti_cone2, Hit::Cone2);
            considerar(ti_cone3, Hit::Cone3);
            considerar(ti_nave, Hit::Nave);
            considerar(ti_esf1_nave, Hit::Esfera1_nave);
            considerar(ti_esf2_nave, Hit::Esfera2_nave);
            considerar(ti_esf3_nave, Hit::Esfera3_nave);

            Color cor(46, 68, 130);  // Cor de fundo (azul escuro)

            switch (hit) {
            case Hit::Fundo:
                cor = plano_fundo.CalcularCor(cena, dr_e);
                break;
            case Hit::Chao:
                cor = plano_chao.CalcularCor(cena, dr_e);
                break;
            case Hit::Cilindro:
                cor = cilindro.CalcularCor(cena, t_best, dr_e);
                break;
            case Hit::Cilindro2:
                cor = cilindro2.CalcularCor(cena, t_best, dr_e);
                break;
            case Hit::Cilindro3:
                cor = cilindro3.CalcularCor(cena, t_best, dr_e);
                break;
            case Hit::Cilindro4:
                cor = cilindro4.CalcularCor(cena, t_best, dr_e);
                break;
            case Hit::Cilindro5:
                cor = cilindro5.CalcularCor(cena, t_best, dr_e);
                break;
            case Hit::Cilindro6:
                cor = cilindro6.CalcularCor(cena, t_best, dr_e);
                break;
            case Hit::Cilindro7:
                cor = cilindro7.CalcularCor(cena, t_best, dr_e);
                break;
            case Hit::Cone:
                cor = cone.CalcularCor(cena, t_best, dr_e);
                break;
            case Hit::Cone2:
                cor = cone2.CalcularCor(cena, t_best, dr_e);
                break;
            case Hit::Cone3:
                cor = cone3.CalcularCor(cena, t_best, dr_e);
                break;
            case Hit::Nave:
            {
                // Se o ponto estiver no disco da base da nave, pinta de branco.
                // (Caso contrário, usa o sombreamento normal do cone.)
                Point Pi = calcula_eq_ray(Po, t_best, dr_e);
                Vector dist_centro = subtrai_pontos(Pi, nave.cb);
                float altura_Pi = calcula_prod_esc(dist_centro, nave.dc);
                float dist2 = calcula_prod_esc(dist_centro, dist_centro);

                const float epsBase = 1e-3f;
                if (std::fabs(altura_Pi) < epsBase && dist2 <= nave.raio * nave.raio + epsBase) {
                    cor = Color(45, 50, 60);
                }
                else {
                    cor = nave.CalcularCor(cena, t_best, dr_e);
                }
                break;
            }
            case Hit::Esfera:
                cor = esfera.CalcularCor(cena, t_best, dr_e);
                break;
            case Hit::Esfera1_nave:
                cor = esfera1_nave.CalcularCor(cena, t_best, dr_e);
                break;
            case Hit::Esfera2_nave:
                cor = esfera2_nave.CalcularCor(cena, t_best, dr_e);
                break;
            case Hit::Esfera3_nave:
                cor = esfera3_nave.CalcularCor(cena, t_best, dr_e);
                break;
                // case Hit::Cubo: {
                //     Vetor PiV = cuboMalha.getPontoIntersecao();
                //     Vetor nV = cuboMalha.calcularNormal(PiV);
                //     Plano plano_cubo(Point(PiV.i, PiV.j, PiV.k), Vector(nV.i, nV.j, nV.k), Material());
                //     plano_cubo.setKa(cuboMalha.getKa());
                //     plano_cubo.setKd(cuboMalha.getKd());
                //     plano_cubo.setKe(cuboMalha.getKe());
                //     plano_cubo.setShininess(cuboMalha.getShininess());
                //     cor = plano_cubo.CalcularCor(cena, dr_e);
                //     break;
                // }
            case Hit::None:
            default:
                break;
            }

            SDL_SetRenderDrawColor(renderer,
                static_cast<int>(cor.r),
                static_cast<int>(cor.g),
                static_cast<int>(cor.b),
                255);
            SDL_RenderDrawPoint(renderer, coluna, linha);

        }
    }

    SDL_RenderPresent(renderer);

    // Escuta eventos para manter janela aberta
    SDL_Event event;
    bool isRunning = true;
    while (isRunning) {
        while (SDL_PollEvent(&event) != 0) {
            if (event.type == SDL_QUIT)
                isRunning = false;
        }
    }

    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}
