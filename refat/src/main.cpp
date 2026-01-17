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
    SDL_SetWindowTitle(window, "Ray Casting - Cena 3D");
    texturaMadeira = new Textura("madeira", "madeira.bmp");

    // Materiais / propriedades via Objeto (Ka/Kd/Ke/shininess)
    esfera.setKa(Vetor(K_e.r, K_e.g, K_e.b));
    esfera.setKd(Vetor(K_e.r, K_e.g, K_e.b));
    esfera.setKe(Vetor(K_e.r, K_e.g, K_e.b));
    esfera.setShininess(m_e);

    // Chão com textura
    mat_chao.usarTextura = true;
    mat_chao.textura = texturaMadeira;
    plano_chao.material = mat_chao;
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

    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);

    // Ray casting
    for (int linha = 0; linha < nLin; linha++) {
        float y = hJanela / 2 - Dy / 2 - linha * Dy;
        for (int coluna = 0; coluna < nCol; coluna++) {
            float x = -wJanela / 2 + Dx / 2 + coluna * Dx;
            Point Pj(x, y, z);

            Vector dr_e = calcula_dr(Po, Pj);

            float t_best = std::numeric_limits<float>::infinity();
            enum class Hit {
                None,
                Fundo,
                Chao,
                Esq,
                Dir,
                Teto,
                Cilindro,
                Cone,
                Cubo,
                Esfera
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
            float ti_f = plano_fundo.CalcularIntersecao(Po, dr_e);
            float ti_e = plano_esq.CalcularIntersecao(Po, dr_e);
            float ti_d = plano_dir.CalcularIntersecao(Po, dr_e);
            float ti_t = plano_teto.CalcularIntersecao(Po, dr_e);

            // Objetos
            float t_cil = cilindro.CalcularIntersecao(Po, dr_e);
            float ti_cone = cone.CalcularIntersecao(Po, Pj);
            float ti_esf = esfera.CalcularIntersecao(Po, dr_e);

            float t_cubo = std::numeric_limits<float>::infinity();
            if (cuboMalha.verificarIntersecao(
                Vetor(Po.x, Po.y, Po.z, 0.0f),
                Vetor(dr_e.i, dr_e.j, dr_e.k, 0.0f))) {
                t_cubo = static_cast<float>(cuboMalha.getDistancia());
            }

            considerar(ti_f, Hit::Fundo);
            considerar(ti_c, Hit::Chao);
            considerar(ti_e, Hit::Esq);
            considerar(ti_d, Hit::Dir);
            considerar(ti_t, Hit::Teto);
            considerar(t_cil, Hit::Cilindro);
            considerar(ti_cone, Hit::Cone);
            considerar(t_cubo, Hit::Cubo);
            considerar(ti_esf, Hit::Esfera);


            Color cor(100, 100, 100);

            switch (hit) {
            case Hit::Fundo:
                cor = plano_fundo.CalcularCor(cena, dr_e);
                break;
            case Hit::Chao:
                cor = plano_chao.CalcularCor(cena, dr_e);
                break;
            case Hit::Esq:
                cor = plano_esq.CalcularCor(cena, dr_e);
                break;
            case Hit::Dir:
                cor = plano_dir.CalcularCor(cena, dr_e);
                break;
            case Hit::Teto:
                cor = plano_teto.CalcularCor(cena, dr_e);
                break;
            case Hit::Cilindro:
                cor = cilindro.CalcularCor(cena, t_best, dr_e);
                break;
            case Hit::Cone:
                cor = cone.CalcularCor(cena, t_best, dr_e);
                break;
            case Hit::Esfera:
                cor = esfera.CalcularCor(cena, t_best, dr_e);
                break;
            case Hit::Cubo: {
                Vetor PiV = cuboMalha.getPontoIntersecao();
                Vetor nV = cuboMalha.calcularNormal(PiV);
                Plano plano_cubo(Point(PiV.i, PiV.j, PiV.k), Vector(nV.i, nV.j, nV.k), Material());
                plano_cubo.setKa(cuboMalha.getKa());
                plano_cubo.setKd(cuboMalha.getKd());
                plano_cubo.setKe(cuboMalha.getKe());
                plano_cubo.setShininess(cuboMalha.getShininess());
                cor = plano_cubo.CalcularCor(cena, dr_e);
                break;
            }
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
