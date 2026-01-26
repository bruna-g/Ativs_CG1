#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <limits>
#include <sstream>
#include <cctype>
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
#include "../include/Objeto.h"
#include "../include/Picking.h"
#include "../include/Render.h"
using namespace std;

// Configuração da janela de visualização
int nCol = 600;
int nLin = 600;

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
// float rEsfera = 5.0f;
// Point centroEsfera(100.f, 150.f + rEsfera, 100.f);
// Color K_e(0.8f, 0.6f, 0.2f);
float m_e = 10.0f;
// Esfera esfera(centroEsfera, rEsfera);

// Textura de madeira
Textura* texturaMadeira = nullptr;
Textura* texturaVaca = nullptr;
Textura* texturaCeu = nullptr;
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
Point P_F(-20.f, 185.f, 450.f);

// Lua (esfera branca no mesmo ponto da luz)
float r_lua = 25.0f;
Esfera lua(P_F, r_lua);

//luz spot
Color I_Spot(1.0f, 0.7f, 0.0f);
Point P_Spot(165.f, 180.f, 125.f);
Vector D_Spot(0.f, -1.f, 0.f);
float angulo_Spot = 30.f;

// Dica: para deixar as árvores com iluminação mais parecida, coloque a luz bem longe
// e aproximadamente centrada no "grupo" de árvores (vira quase uma luz direcional).
// Colocar a luz abaixo da base das copas (y < 60) evita que o cone projete sombra no tronco.
// E deixar bem longe (z grande) reduz variação de direção entre as 3 árvores.
// Point P_F(175.f, 180.f, 500.f);

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

//-----------------Arvore 4---------------------//
// Cilindro (tronco da árvore)
Point centroCilindro8(25.f, 00.f, 0.f);
float raio_cil8 = 8.f;
float altura_cil8 = 60.f;
Vector dc8(0.f, 1.0f, 0.f);
Cilindro cilindro8(centroCilindro8, raio_cil8, altura_cil8, dc8);
Color KCil_d8(0.545f, 0.270f, 0.075f);
Color KCil_e8(0.545f, 0.270f, 0.075f);
Color KCil_a8(0.545f, 0.270f, 0.075f);
Material mat_cil8;

// Cone (copa da árvore)
Point cb_cone8(25.f, 60.f, 0.f);
float raio_cone8 = 40.f;
float altura_cone8 = 50.f;
Vector aux_v_cone8 = calcula_esc_por_vetor(altura_cone8, dc8);
Point v_cone8(cb_cone8.x + aux_v_cone8.i, cb_cone8.y + aux_v_cone8.j, cb_cone8.z + aux_v_cone8.k);
Cone cone8(cb_cone8, v_cone8, raio_cone8);
Color KCone8(0.f, 0.6f, 0.2f);
Material mat_cone8;
//----------------------------------------------//

//-----------------Arvore 5---------------------//
// Cilindro (tronco da árvore)
Point centroCilindro9(325.f, 00.f, 0.f);
float raio_cil9 = 8.f;
float altura_cil9 = 60.f;
Vector dc9(0.f, 1.0f, 0.f);
Cilindro cilindro9(centroCilindro9, raio_cil9, altura_cil9, dc9);
Color KCil_d9(0.545f, 0.270f, 0.075f);
Color KCil_e9(0.545f, 0.270f, 0.075f);
Color KCil_a9(0.545f, 0.270f, 0.075f);
Material mat_cil9;

// Cone (copa da árvore)
Point cb_cone9(325.f, 60.f, 0.f);
float raio_cone9 = 40.f;
float altura_cone9 = 50.f;
Vector aux_v_cone9 = calcula_esc_por_vetor(altura_cone9, dc9);
Point v_cone9(cb_cone9.x + aux_v_cone9.i, cb_cone9.y + aux_v_cone9.j, cb_cone9.z + aux_v_cone9.k);
Cone cone9(cb_cone9, v_cone9, raio_cone9);
Color KCone9(0.f, 0.6f, 0.2f);
Material mat_cone9;
//----------------------------------------------//

//-----------------Vaca---------------------//
// Cilindro
Point centroCilindro4(150.f, 20.f, 110.f);
float raio_cil4 = 4.f;
float altura_cil4 = 20.f;
Vector dc4(0.f, 1.0f, 0.f);
Cilindro cilindro4(centroCilindro4, raio_cil4, altura_cil4, dc4);
Color KCil_d4(0.8f, 0.8f, 0.8f);
Color KCil_e4(0.8f, 0.8f, 0.8f);
Color KCil_a4(0.8f, 0.8f, 0.8f);
Material mat_cil4;

Point centroCilindro5(150.f, 20.f, 140.f);
float raio_cil5 = 4.f;
float altura_cil5 = 20.f;
Vector dc5(0.f, 1.0f, 0.f);
Cilindro cilindro5(centroCilindro5, raio_cil5, altura_cil5, dc5);
Color KCil_d5(0.8f, 0.8f, 0.8f);
Color KCil_e5(0.8f, 0.8f, 0.8f);
Color KCil_a5(0.8f, 0.8f, 0.8f);
Material mat_cil5;

Point centroCilindro6(125.f, 45.f, 110.f);
float raio_cil6 = 4.f;
float altura_cil6 = 20.f;
Vector dc6(0.f, 1.0f, 0.f);
Cilindro cilindro6(centroCilindro6, raio_cil6, altura_cil6, dc6);
Color KCil_d6(0.8f, 0.8f, 0.8f);
Color KCil_e6(0.8f, 0.8f, 0.8f);
Color KCil_a6(0.8f, 0.8f, 0.8f);
Material mat_cil6;

Point centroCilindro7(125.f, 45.f, 140.f);
float raio_cil7 = 4.f;
float altura_cil7 = 20.f;
Vector dc7(0.f, 1.0f, 0.f);
Cilindro cilindro7(centroCilindro7, raio_cil7, altura_cil7, dc7);
Color KCil_d7(0.8f, 0.8f, 0.8f);
Color KCil_e7(0.8f, 0.8f, 0.8f);
Color KCil_a7(0.8f, 0.8f, 0.8f);
Material mat_cil7;

// Parâmetros do cubo (corpo da vaca)
const Point cubo_centro(170.f, 40.f, 125.f);
const float cubo_lado = 35.0f;
const float cubo_escala_x = 1.8f;
const float cubo_escala_y = 1.0f;
const float cubo_escala_z = 1.2f;

// Esfera cabeça (encostada na face esquerda do cubo após a escala)
float rEsfera_cabeca = 15.0f;
Point centroEsfera_cabeca(135.f, 90.f, 125.f);
Color K_e_cabeca(1.0f, 1.0f, 1.0f);
Color K_d_cabeca(1.0f, 1.0f, 1.0f);
Color K_a_cabeca(1.0f, 1.0f, 1.0f);
float m_e_cabeca = 10.0f;
Esfera esfera_cabeca(centroEsfera_cabeca, rEsfera_cabeca);

// Chifres (2 cones pequenos na cabeça)
const float raio_chifre = 2.5f;
const float altura_chifre = 8.0f;
Vector dc_chifre(0.f, 1.0f, 0.f);

Point cb_chifre_esq(
    centroEsfera_cabeca.x + rEsfera_cabeca * 0.20f,
    centroEsfera_cabeca.y + rEsfera_cabeca * 0.80f,
    centroEsfera_cabeca.z - rEsfera_cabeca * 0.35f);
Point cb_chifre_dir(
    centroEsfera_cabeca.x + rEsfera_cabeca * 0.20f,
    centroEsfera_cabeca.y + rEsfera_cabeca * 0.80f,
    centroEsfera_cabeca.z + rEsfera_cabeca * 0.35f);

Vector aux_v_chifre = calcula_esc_por_vetor(altura_chifre, dc_chifre);
Point v_chifre_esq(cb_chifre_esq.x + aux_v_chifre.i, cb_chifre_esq.y + aux_v_chifre.j, cb_chifre_esq.z + aux_v_chifre.k);
Point v_chifre_dir(cb_chifre_dir.x + aux_v_chifre.i, cb_chifre_dir.y + aux_v_chifre.j, cb_chifre_dir.z + aux_v_chifre.k);

Cone chifre_esq(cb_chifre_esq, v_chifre_esq, raio_chifre);
Cone chifre_dir(cb_chifre_dir, v_chifre_dir, raio_chifre);

Color K_chifre(0.9f, 0.68f, 0.55f);
float m_chifre = 25.0f;

// Cauda (cilindro atrás do corpo, levemente inclinado para baixo)
const float raio_cauda = 1.6f;
const float altura_cauda = 18.0f;
const float cubo_meia_largura_x = (cubo_lado * 0.5f) * cubo_escala_x;
Point cb_cauda(
    cubo_centro.x + cubo_meia_largura_x - 0.5f,
    cubo_centro.y + 10.0f,
    cubo_centro.z);
Vector dc_cauda(0.f, 1.f, 0.0f);
Cilindro cauda(cb_cauda, raio_cauda, altura_cauda, dc_cauda);
Color K_cauda(1.f, 1.f, 1.f);
float m_cauda = 15.0f;

// Cubo (como malha)
Malha cuboMalha;

Matriz3x3 M_id(1.0f, 0.0f, 0.0f,
    0.0f, 1.0f, 0.0f,
    0.0f, 0.0f, 1.0f);
//--------------------------------------//

// nave
Point cb_nave(165.f, 200.f, 125.f);
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
Point centroEsfera1_nave(140.f, 205.f + rEsfera1_nave, 185.f);
Color K_e1_nave(1.f, 0.f, 0.f);  // Vermelho
Color K_d1_nave(1.f, 0.f, 0.f);  // Vermelho
Color K_a1_nave(1.f, 0.1f, 0.f);  // Vermelho
float m_e1_nave = 10.0f;
Esfera esfera1_nave(centroEsfera1_nave, rEsfera1_nave);

// Esfera2_nave
float rEsfera2_nave = 4.0f;
Point centroEsfera2_nave(165.f, 205.f + rEsfera2_nave, 191.6f);
Color K_e2_nave(1.f, 0.f, 0.f);  // Vermelho
Color K_d2_nave(1.f, 0.f, 0.f);  // Vermelho
Color K_a2_nave(1.f, 0.1f, 0.f);  // Vermelho
float m_e2_nave = 10.0f;
Esfera esfera2_nave(centroEsfera2_nave, rEsfera2_nave);

// Esfera3_nave
float rEsfera3_nave = 4.0f;
Point centroEsfera3_nave(190.f, 205.f + rEsfera3_nave, 185.f);
Color K_e3_nave(1.f, 0.f, 0.f);  // Vermelho
Color K_d3_nave(1.f, 0.f, 0.f);  // Vermelho
Color K_a3_nave(1.f, 0.1f, 0.f);  // Vermelho
float m_e3_nave = 10.0f;
Esfera esfera3_nave(centroEsfera3_nave, rEsfera3_nave);



int main() {
    SDL_Window* window = nullptr;
    SDL_Renderer* renderer = nullptr;

    SDL_Init(SDL_INIT_VIDEO);
    SDL_CreateWindowAndRenderer(nCol, nLin, 0, &window, &renderer);
    SDL_SetWindowTitle(window, "Trabalho Final - CG");
    //texturaMadeira = new Textura("madeira", "madeira.bmp");
    texturaVaca = new Textura("vaca", "vaca.bmp");
    texturaCeu = new Textura("ceu", "ceu.bmp");
    texturaMadeira = new Textura("grama", "grama.bmp");

    // Materiais / propriedades via Objeto (Ka/Kd/Ke/shininess)
    // esfera.setKa(Vetor(K_e.r, K_e.g, K_e.b));
    // esfera.setKd(Vetor(K_e.r, K_e.g, K_e.b));
    // esfera.setKe(Vetor(K_e.r, K_e.g, K_e.b));
    // esfera.setShininess(m_e);

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

    // esfera 3 da nave
    esfera_cabeca.setKa(Vetor(K_a_cabeca.r, K_a_cabeca.g, K_a_cabeca.b));
    esfera_cabeca.setKd(Vetor(K_d_cabeca.r, K_d_cabeca.g, K_d_cabeca.b));
    esfera_cabeca.setKe(Vetor(K_e_cabeca.r, K_e_cabeca.g, K_e_cabeca.b));
    esfera_cabeca.setShininess(m_e_cabeca);

    // Lua (branca)
    lua.setKa(Vetor(4.0f, 4.0f, 4.0f));
    lua.setKd(Vetor(1.0f, 1.0f, 1.0f));
    lua.setKe(Vetor(0.0f, 0.0f, 0.0f));
    lua.setShininess(10.0f);

    // Textura da vaca (cabeça)
    // esfera_cabeca.material.usarTextura = true;
    // esfera_cabeca.material.textura = texturaVaca;
    // esfera_cabeca.material.texturaScale = 1.0f;

    // chifres
    chifre_esq.setKa(Vetor(K_chifre.r, K_chifre.g, K_chifre.b));
    chifre_esq.setKd(Vetor(K_chifre.r, K_chifre.g, K_chifre.b));
    chifre_esq.setKe(Vetor(K_chifre.r, K_chifre.g, K_chifre.b));
    chifre_esq.setShininess(m_chifre);

    chifre_dir.setKa(Vetor(K_chifre.r, K_chifre.g, K_chifre.b));
    chifre_dir.setKd(Vetor(K_chifre.r, K_chifre.g, K_chifre.b));
    chifre_dir.setKe(Vetor(K_chifre.r, K_chifre.g, K_chifre.b));
    chifre_dir.setShininess(m_chifre);

    // cauda
    cauda.setKa(Vetor(K_cauda.r, K_cauda.g, K_cauda.b));
    cauda.setKd(Vetor(K_cauda.r, K_cauda.g, K_cauda.b));
    cauda.setKe(Vetor(K_cauda.r, K_cauda.g, K_cauda.b));
    cauda.setShininess(m_cauda);

    // // Textura da vaca (cauda)
    // cauda.material.usarTextura = true;
    // cauda.material.textura = texturaVaca;
    // cauda.material.texturaScale = 1.0f;

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

    cilindro4.setKa(Vetor(KCil_a4.r, KCil_a4.g, KCil_a4.b));
    cilindro4.setKd(Vetor(KCil_d4.r, KCil_d4.g, KCil_d4.b));
    cilindro4.setKe(Vetor(KCil_e4.r, KCil_e4.g, KCil_e4.b));
    cilindro4.setShininess(m_e);

    cilindro5.setKa(Vetor(KCil_a5.r, KCil_a5.g, KCil_a5.b));
    cilindro5.setKd(Vetor(KCil_d5.r, KCil_d5.g, KCil_d5.b));
    cilindro5.setKe(Vetor(KCil_e5.r, KCil_e5.g, KCil_e5.b));
    cilindro5.setShininess(m_e);

    cilindro6.setKa(Vetor(KCil_a6.r, KCil_a6.g, KCil_a6.b));
    cilindro6.setKd(Vetor(KCil_d6.r, KCil_d6.g, KCil_d6.b));
    cilindro6.setKe(Vetor(KCil_e6.r, KCil_e6.g, KCil_e6.b));
    cilindro6.setShininess(m_e);

    cilindro7.setKa(Vetor(KCil_a7.r, KCil_a7.g, KCil_a7.b));
    cilindro7.setKd(Vetor(KCil_d7.r, KCil_d7.g, KCil_d7.b));
    cilindro7.setKe(Vetor(KCil_e7.r, KCil_e7.g, KCil_e7.b));
    cilindro7.setShininess(m_e);

    cilindro8.setKa(Vetor(KCil_a8.r, KCil_a8.g, KCil_a8.b));
    cilindro8.setKd(Vetor(KCil_d8.r, KCil_d8.g, KCil_d8.b));
    cilindro8.setKe(Vetor(KCil_e8.r, KCil_e8.g, KCil_e8.b));
    cilindro8.setShininess(m_e);

    cilindro9.setKa(Vetor(KCil_a9.r, KCil_a9.g, KCil_a9.b));
    cilindro9.setKd(Vetor(KCil_d9.r, KCil_d9.g, KCil_d9.b));
    cilindro9.setKe(Vetor(KCil_e9.r, KCil_e9.g, KCil_e9.b));
    cilindro9.setShininess(m_e);

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

    cone8.setKa(Vetor(KCone8.r, KCone8.g, KCone8.b));
    cone8.setKd(Vetor(KCone8.r, KCone8.g, KCone8.b));
    cone8.setKe(Vetor(KCone8.r, KCone8.g, KCone8.b));
    cone8.setShininess(m_e);

    cone9.setKa(Vetor(KCone9.r, KCone9.g, KCone9.b));
    cone9.setKd(Vetor(KCone9.r, KCone9.g, KCone9.b));
    cone9.setKe(Vetor(KCone9.r, KCone9.g, KCone9.b));
    cone9.setShininess(m_e);

    nave.setKa(Vetor(Ka_nave.r, Ka_nave.g, Ka_nave.b));
    nave.setKd(Vetor(Kd_nave.r, Kd_nave.g, Kd_nave.b));
    nave.setKe(Vetor(Ke_nave.r, Ke_nave.g, Ke_nave.b));
    nave.setShininess(80.0f);

    // Cubo como malha
    Color K_cubo(0.45f, 0.25f, 0.10f);
    cuboMalha = Cubo::criarCubo(
        Vetor(cubo_centro.x, cubo_centro.y, cubo_centro.z, 0.0f),
        cubo_lado,
        Vetor(K_cubo.r, K_cubo.g, K_cubo.b, 0.0f),
        Vetor(K_cubo.r, K_cubo.g, K_cubo.b, 0.0f),
        Vetor(K_cubo.r, K_cubo.g, K_cubo.b, 0.0f),
        m_e
    );
    cuboMalha.aplicarEscalaNoPivoObjeto(
        Vetor(cubo_escala_x, cubo_escala_y, cubo_escala_z),
        cubo_centro);

    cuboMalha.rotacionarZ(-45.0f);

    // Textura da vaca (corpo / cuboMalha) será aplicada no Render.cpp

    // Vetor deslocCone(20.f, 0.f, 0.f); // +X => direita
    // cone.cb = cone.aplicarTranslacao(cone.cb, deslocCone);
    // cone.v = cone.aplicarTranslacao(cone.v, deslocCone);
    // cone.recalcularDerivados(); // opcional, mas seguro

    cilindro4.rotacionarZ(-45.0f);
    cilindro5.rotacionarZ(-45.0f);
    cilindro6.rotacionarZ(-45.0f);
    cilindro7.rotacionarZ(-45.0f);
    cauda.rotacionarZ(-45.0f);

    chifre_esq.rotacionarZ(-45.0f);
    chifre_dir.rotacionarZ(-45.0f);



    // Vetor deslocCone(20.f, 0.f, 0.f); // +X => direita
    // cone.cb = cone.aplicarTranslacao(cone.cb, deslocCone);
    // cone.v = cone.aplicarTranslacao(cone.v, deslocCone);
    // cone.recalcularDerivados(); // opcional, mas seguro


    Cena cena;
    cena.observador = camera.getEye();
    cena.luz = LuzPontual{ P_F, I_F };
    cena.luzSpot = LuzSpot{ I_Spot, P_Spot, D_Spot, angulo_Spot };
    cena.luzSpotAtiva = false;
    cena.luzAmbiente = I_A;
    cena.objetosSombra = { &cilindro, &cone, &cilindro2, &cone2, &cilindro3, &cone3, &nave,
        &esfera1_nave, &esfera2_nave, &esfera3_nave, &cilindro4, &cilindro5, &cilindro6, &cilindro7, &cilindro8, &cuboMalha, &esfera_cabeca,
        &chifre_esq, &chifre_dir, &cauda, &cone8, &cilindro9, &cone8, &cone9 };
    cena.texturaMadeira = texturaMadeira;
    cena.expoenteEspecular = m_e;
    gCena = &cena;

    renderScene(renderer, cena);

    // Escuta eventos para manter janela aberta
    SDL_Event event;
    bool isRunning = true;
    while (isRunning) {
        while (SDL_PollEvent(&event) != 0) {
            if (event.type == SDL_QUIT)
                isRunning = false;

            if (event.type == SDL_MOUSEBUTTONDOWN && event.button.button == SDL_BUTTON_LEFT) {
                gSelecionado = pickFromScreen(event.button.x, event.button.y, nCol, nLin);
                imprimirSelecaoDetalhada(gSelecionado);
                if (gSelecionado.hit != Hit::None) {
                    imprimirAjudaSelecao();
                }
            }

            if (event.type == SDL_KEYDOWN) {
                // Evita repetir comando quando a tecla fica pressionada (auto-repeat do SDL).
                if (event.key.repeat != 0) {
                    continue;
                }

                if (event.key.keysym.sym == SDLK_h) {
                    imprimirAjudaSelecao();
                    continue;
                }

                if (gSelecionado.hit == Hit::None || gSelecionado.objeto == nullptr) {
                    if (event.key.keysym.sym == SDLK_c || event.key.keysym.sym == SDLK_m ||
                        event.key.keysym.sym == SDLK_t || event.key.keysym.sym == SDLK_r ||
                        event.key.keysym.sym == SDLK_s || event.key.keysym.sym == SDLK_v) {
                        std::cout << "Nenhum objeto selecionado. Clique em um objeto primeiro." << std::endl;
                    }
                    continue;
                }

                if (event.key.keysym.sym == SDLK_c) {
                    float r = 0.0f, g = 0.0f, b = 0.0f;
                    if (ler3Floats("Digite cor (r g b) em [0..1] ou [0..255]: ", r, g, b)) {
                        setCorObjeto(gSelecionado.objeto, r, g, b);
                        renderScene(renderer, cena);
                        refocarJanela(window);
                    }
                    else {
                        std::cout << "Entrada invalida. Exemplo: 0.2 0.7 0.2" << std::endl;
                        refocarJanela(window);
                    }
                }
                else if (event.key.keysym.sym == SDLK_m) {
                    float ka_r = 0.0f, ka_g = 0.0f, ka_b = 0.0f;
                    float ke_r = 0.0f, ke_g = 0.0f, ke_b = 0.0f;
                    float shininess = 1.0f;

                    if (!ler3Floats("Ka (r g b) em [0..1] ou [0..255]: ", ka_r, ka_g, ka_b)) {
                        std::cout << "Entrada invalida. Exemplo: 0.2 0.2 0.2" << std::endl;
                        refocarJanela(window);
                        continue;
                    }

                    if (!ler3Floats("Ke (r g b) em [0..1] ou [0..255]: ", ke_r, ke_g, ke_b)) {
                        std::cout << "Entrada invalida. Exemplo: 0.0 0.0 0.0" << std::endl;
                        refocarJanela(window);
                        continue;
                    }

                    if (!lerFloat("Shininess (> 0): ", shininess) || shininess <= 0.0f) {
                        std::cout << "Entrada invalida. Exemplo: 30" << std::endl;
                        refocarJanela(window);
                        continue;
                    }

                    aplicarMaterialCustom(gSelecionado.objeto,
                        ka_r, ka_g, ka_b,
                        ke_r, ke_g, ke_b,
                        shininess);
                    renderScene(renderer, cena);
                    refocarJanela(window);
                }
                else if (event.key.keysym.sym == SDLK_t) {
                    float dx = 0.0f, dy = 0.0f, dz = 0.0f;
                    if (ler3Floats("Translacao (dx dy dz): ", dx, dy, dz)) {
                        aplicarTranslacaoSelecionado(gSelecionado.hit, Vetor(dx, dy, dz, 0.0f));
                        renderScene(renderer, cena);
                        refocarJanela(window);
                    }
                    else {
                        std::cout << "Entrada invalida. Exemplo: 10 0 -5" << std::endl;
                        refocarJanela(window);
                    }
                }
                else if (event.key.keysym.sym == SDLK_r) {
                    char eixo = 'z';
                    float ang = 0.0f;
                    if (lerCharEFloat("Rotacao: eixo (x/y/z) e angulo (graus): ", eixo, ang)) {
                        eixo = static_cast<char>(std::tolower(static_cast<unsigned char>(eixo)));
                        if (eixo != 'x' && eixo != 'y' && eixo != 'z') eixo = 'z';
                        aplicarRotacaoSelecionado(gSelecionado.hit, eixo, ang);
                        renderScene(renderer, cena);
                        refocarJanela(window);
                    }
                    else {
                        std::cout << "Entrada invalida. Exemplo: z 15" << std::endl;
                        refocarJanela(window);
                    }
                }
                else if (event.key.keysym.sym == SDLK_s) {
                    float s = 1.0f;
                    if (lerFloat("Escala uniforme (fator > 0): ", s)) {
                        aplicarEscalaSelecionado(gSelecionado.hit, s);
                        renderScene(renderer, cena);
                        refocarJanela(window);
                    }
                    else {
                        std::cout << "Entrada invalida. Exemplo: 1.25" << std::endl;
                        refocarJanela(window);
                    }
                }
                else if (event.key.keysym.sym == SDLK_v) {
                    float sx = 1.0f, sy = 1.0f, sz = 1.0f;
                    if (!ler3Floats("Escala vetorial (sx sy sz). Use 1 para manter um eixo: ", sx, sy, sz)) {
                        std::cout << "Entrada invalida. Exemplo: 2 1 1" << std::endl;
                        refocarJanela(window);
                        continue;
                    }

                    if (sx <= 0.0f || sy <= 0.0f || sz <= 0.0f) {
                        std::cout << "Escala vetorial invalida: sx, sy, sz precisam ser > 0. Exemplo: 3 1 1" << std::endl;
                        refocarJanela(window);
                        continue;
                    }

                    aplicarEscalaVetorSelecionado(gSelecionado.hit, Vetor(sx, sy, sz, 0.0f));
                    renderScene(renderer, cena);
                    refocarJanela(window);
                }
            }
        }
    }

    gCena = nullptr;

    delete texturaMadeira;
    texturaMadeira = nullptr;
    delete texturaVaca;
    texturaVaca = nullptr;
    delete texturaCeu;
    texturaCeu = nullptr;

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}
