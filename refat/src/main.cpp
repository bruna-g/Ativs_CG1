#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
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
using namespace std;

// Configuração da janela de visualização
float wJanela = 60.0f;
float hJanela = 60.0f;
int nCol = 500;
int nLin = 500;

float dJanela = 30.0f;
float z = -dJanela;

float Dx = wJanela / nCol;
float Dy = hJanela / nLin;

Point Po(0.f, 0.f, 0.f);

// Esfera
float rEsfera = 5.0f;
Point centroEsfera(0.f, 95.f, -200.f);
Color K_e(0.854f, 0.647f, 0.125f);
float m_e = 10.0f;
Esfera esfera(centroEsfera, rEsfera);

// Textura de madeira
Textura* texturaMadeira = nullptr;
static Cena* gCena = nullptr;

// Chão
Point P_pi_chao(0.f, -150.f, 0.f);
Vector n_chao(0.f, 1.0f, 0.f);
Color KC_d(0.2f, 0.7f, 0.2f);
Color KC_e(0.0f, 0.0f, 0.0f);
float m_c = 1.0f;
Material mat_chao;
Plano plano_chao(P_pi_chao, n_chao, mat_chao);

// Fundo
Point P_pi_fundo(200.f, -150.f, -400.f);
Vector n_fundo(0.f, 0.f, 1.f);
Color KF_d(0.686f, 0.933f, 0.933f);
Color KF_e(0.686f, 0.933f, 0.933f);
float m_f = 1.0f;
Material mat_fundo;
Plano plano_fundo(P_pi_fundo, n_fundo, mat_fundo);

// Parede esquerda
Point P_pi_esq(-200.f, -150.f, 0.f);
Vector n_esq(1.f, 0.f, 0.f);
Color KE_d(0.686f, 0.933f, 0.933f);
Color KE_e(0.686f, 0.933f, 0.933f);
Material mat_esq;
Plano plano_esq(P_pi_esq, n_esq, mat_esq);

// Parede direita
Point P_pi_dir(200.f, -150.f, 0.f);
Vector n_dir(-1.f, 0.f, 0.f);
Color KD_d(0.686f, 0.933f, 0.933f);
Color KD_e(0.686f, 0.933f, 0.933f);
Material mat_dir;
Plano plano_dir(P_pi_dir, n_dir, mat_dir);

// Teto
Point P_pi_teto(0.f, 150.f, 0.f);
Vector n_teto(0.f, -1.f, 0.f);
Color KT_d(0.933f, 0.933f, 0.933f);
Color KT_e(0.933f, 0.933f, 0.933f);
Material mat_teto;
Plano plano_teto(P_pi_teto, n_teto, mat_teto);

// Fonte luminosa
Color I_F(0.7f, 0.7f, 0.7f);
Point P_F(-100.f, 140.f, -20.f);

// Luz ambiente
Color I_A(0.3f, 0.3f, 0.3f);

// Cilindro
Point centroCilindro(0.f, -150.f, -200.f);
float raio_cil = 5.f;
float altura_cil = 90.f;
Vector dc(0.f, 1.0f, 0.f);
Cilindro cilindro(centroCilindro, raio_cil, altura_cil, dc);
Color KCil_d(0.824f, 0.706f, 0.549f);
Color KCil_e(0.824f, 0.706f, 0.549f);
Color KCil_a(0.824f, 0.706f, 0.549f);
Material mat_cil;

// Cone
Point cb_cone(0.f, -60.f, -200.f);
float raio_cone = 90.f;
float altura_cone = 150.f;
Vector aux_v_cone = calcula_esc_por_vetor(altura_cone, dc);
Point v_cone(cb_cone.x + aux_v_cone.i, cb_cone.y + aux_v_cone.j, cb_cone.z + aux_v_cone.k);
Cone cone(cb_cone, v_cone, raio_cone);
Color KCone(0.f, 1.f, 0.498f);
Material mat_cone;

// Cubo
Cubo cubo(Point(0.f, -150.f, -165.f), 40.f);

Matriz3x3 M_id(1.0f, 0.0f, 0.0f,
    0.0f, 1.0f, 0.0f,
    0.0f, 0.0f, 1.0f);

int main() {
    SDL_Window *window = nullptr;
    SDL_Renderer *renderer = nullptr;
    
    SDL_Init(SDL_INIT_VIDEO);
    SDL_CreateWindowAndRenderer(nCol, nLin, 0, &window, &renderer);
    SDL_SetWindowTitle(window, "Ray Casting - Cena 3D");
    texturaMadeira = new Textura("madeira", "madeira.bmp");


    // Inicializar materiais
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
    cena.texturaMadeira = nullptr;
    cena.expoenteEspecular = m_e;
    gCena = &cena;

    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);

    // Ray casting
    for (int linha = 0; linha < nLin; linha++) {
        float y = hJanela / 2 - Dy / 2 - linha * Dy;
        for (int coluna = 0; coluna < nCol; coluna++) {
            float x = -wJanela / 2 + Dx / 2 + coluna * Dx;
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
            float ti_esf = esfera.CalcularIntersecao(Po,dr_e);

            IntersecaoCubo inter_cubo = cubo.CalcularIntersecaoCompleta(Po, dr_e);
            float t_cubo = inter_cubo.intercepta ? inter_cubo.t : -1.0f;

            if (ti_c > 0.0f && (ti_c < t || t < 0.0f)) t = ti_c;
            if (ti_f > 0.0f && (ti_f < t || t < 0.0f)) t = ti_f;
            if (ti_e > 0.0f && (ti_e < t || t < 0.0f)) t = ti_e;
            if (ti_d > 0.0f && (ti_d < t || t < 0.0f)) t = ti_d;
            if (ti_t > 0.0f && (ti_t < t || t < 0.0f)) t = ti_t;
            if (t_cil > 0.0f && (t_cil < t || t < 0.0f)) t = t_cil;
            if (ti_cone > 0.0f && (ti_cone < t || t < 0.0f)) t = ti_cone;
            if (t_cubo > 0.0f && (t_cubo < t || t < 0.0f)) t = t_cubo;
            if (ti_esf > 0.0f && (ti_esf < t || t < 0.0f)) t = ti_esf;


            Color cor(100, 100, 100);

            if (ti_f == t) {
                cor = plano_fundo.CalcularCor(cena, dr_e);
            } else if (ti_c == t) {
                cor = plano_chao.CalcularCor(cena, dr_e);
            } else if (ti_e == t) {
                cor = plano_esq.CalcularCor(cena, dr_e);
            } else if (ti_d == t) {
                cor = plano_dir.CalcularCor(cena, dr_e);
            } else if (ti_t == t) {
                cor = plano_teto.CalcularCor(cena, dr_e);
            } else if (t_cil == t) {
                cor = cilindro.CalcularCor(cena, t, dr_e);
            } else if (ti_cone == t) {
                cor = cone.CalcularCor(cena, t, dr_e);
            } else if (t_cubo == t) {
                Color K_cubo(1.f, 0.078f, 0.576f);
                Plano plano_cubo(inter_cubo.ponto, inter_cubo.normal, Material::Solido(K_cubo, m_e));
                cor = plano_cubo.CalcularCor(cena, dr_e);
            } else if (ti_esf == t){
                cor = esfera.CalcularCor(cena, t, dr_e);
            } else if (t > 0.0f && delta >= 0.0f) {
                SDL_SetRenderDrawColor(renderer, 
                static_cast<int>(cor.r ), 
                static_cast<int>(cor.g ), 
                static_cast<int>(cor.b ),
                255);
                SDL_RenderDrawPoint(renderer, coluna, linha);
            }

            SDL_SetRenderDrawColor(renderer, 
                static_cast<int>(cor.r ), 
                static_cast<int>(cor.g ), 
                static_cast<int>(cor.b ),
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
