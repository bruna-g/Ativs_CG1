#include "Render.h"

#include "Picking.h"

#include <SDL2/SDL.h>

#include "Camera.h"
#include "Plano.h"
#include "Cilindro.h"
#include "Cone.h"
#include "Esfera.h"
#include "Malha.h"
#include "Material.h"
#include "Textura.hpp"
#include "Utils.h"

// Dependências globais (continuam sendo definidas em main.cpp)
extern int nCol;
extern int nLin;
extern float Dx;
extern float Dy;
extern Camera camera;

extern Plano plano_fundo;
extern Plano plano_chao;

extern Cilindro cilindro;
extern Cilindro cilindro2;
extern Cilindro cilindro3;
extern Cilindro cilindro4;
extern Cilindro cilindro5;
extern Cilindro cilindro6;
extern Cilindro cilindro7;
extern Cilindro cauda;

extern Cone cone;
extern Cone cone2;
extern Cone cone3;
extern Cone chifre_esq;
extern Cone chifre_dir;
extern Cone nave;

extern Esfera esfera1_nave;
extern Esfera esfera2_nave;
extern Esfera esfera3_nave;
extern Esfera esfera_cabeca;

extern Malha cuboMalha;

// Textura da vaca (definida em main.cpp)
extern Textura* texturaVaca;

void renderScene(SDL_Renderer* renderer, const Cena& cena) {
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);

    for (int linha = 0; linha < nLin; linha++) {
        float y = camera.getYMax() - (linha + 0.5f) * Dy;
        for (int coluna = 0; coluna < nCol; coluna++) {
            float x = camera.getXMin() + (coluna + 0.5f) * Dx;

            Vector dr_e = camera.gerarRaio(x, y);
            Point Po = camera.getEye();

            PickResult pr = pickRay(Po, dr_e);
            float t_best = pr.t;
            Hit hit = pr.hit;

            Color cor(46, 68, 130);  // fundo

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
            case Hit::Cauda:
                cor = cauda.CalcularCor(cena, t_best, dr_e);
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
            case Hit::Chifre_esq:
                cor = chifre_esq.CalcularCor(cena, t_best, dr_e);
                break;
            case Hit::Chifre_dir:
                cor = chifre_dir.CalcularCor(cena, t_best, dr_e);
                break;
            case Hit::Nave: {
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
            case Hit::Esfera1_nave:
                cor = esfera1_nave.CalcularCor(cena, t_best, dr_e);
                break;
            case Hit::Esfera2_nave:
                cor = esfera2_nave.CalcularCor(cena, t_best, dr_e);
                break;
            case Hit::Esfera3_nave:
                cor = esfera3_nave.CalcularCor(cena, t_best, dr_e);
                break;
            case Hit::Cubo: {
                Vetor PiV = cuboMalha.getPontoIntersecao();
                Vetor nV = cuboMalha.calcularNormal(PiV);

                // Cria um plano auxiliar (mesma face atingida) para reaproveitar o shading do Plano,
                // mas com um ponto de referência estável para o mapeamento de textura.
                // Para um plano com normal n e passando por Pi: d = n·Pi.
                // Escolhe um ponto p_pi no plano resolvendo uma coordenada (a de maior módulo).
                const Vector nFace(nV.i, nV.j, nV.k);
                const float d = nFace.i * PiV.i + nFace.j * PiV.j + nFace.k * PiV.k;
                Point p_pi(0.0f, 0.0f, 0.0f);
                const float ax = std::fabs(nFace.i);
                const float ay = std::fabs(nFace.j);
                const float az = std::fabs(nFace.k);
                if (az >= ax && az >= ay && std::fabs(nFace.k) > 1e-6f) {
                    p_pi = Point(0.0f, 0.0f, d / nFace.k);
                }
                else if (ay >= ax && std::fabs(nFace.j) > 1e-6f) {
                    p_pi = Point(0.0f, d / nFace.j, 0.0f);
                }
                else if (std::fabs(nFace.i) > 1e-6f) {
                    p_pi = Point(d / nFace.i, 0.0f, 0.0f);
                }

                Material mat_vaca;
                mat_vaca.usarTextura = (texturaVaca != nullptr);
                mat_vaca.textura = texturaVaca;
                // Para o cubo, deixa a textura “1x por face” (pode ajustar depois)
                mat_vaca.texturaScale = 0.02f;

                Plano plano_cubo(p_pi, nFace, mat_vaca);
                plano_cubo.setKa(cuboMalha.getKa());
                plano_cubo.setKd(cuboMalha.getKd());
                plano_cubo.setKe(cuboMalha.getKe());
                plano_cubo.setShininess(cuboMalha.getShininess());
                cor = plano_cubo.CalcularCor(cena, dr_e);
                break;
            }
            case Hit::Esfera_cabeca:
                cor = esfera_cabeca.CalcularCor(cena, t_best, dr_e);
                break;
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
}
