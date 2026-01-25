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

#include <cmath>

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
extern Cilindro cilindro8;
extern Cilindro cilindro9;
extern Cilindro cauda;

extern Cone cone;
extern Cone cone2;
extern Cone cone3;
extern Cone cone8;
extern Cone cone9;
extern Cone chifre_esq;
extern Cone chifre_dir;
extern Cone nave;

extern Esfera esfera1_nave;
extern Esfera esfera2_nave;
extern Esfera esfera3_nave;
extern Esfera esfera_cabeca;
extern Esfera lua;

extern Malha cuboMalha;

// Textura da vaca (definida em main.cpp)
extern Textura* texturaVaca;
extern Textura* texturaCeu;

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

            // switch (hit) {
            // case Hit::Fundo:
            //     cor = plano_fundo.CalcularCor(cena, dr_e);
            //     break;
            // case Hit::Chao:
            //     cor = plano_chao.CalcularCor(cena, dr_e);
            //     break;
            // case Hit::Cilindro:
            //     cor = cilindro.CalcularCor(cena, t_best, dr_e);
            //     break;
            // case Hit::Cilindro2:
            //     cor = cilindro2.CalcularCor(cena, t_best, dr_e);
            //     break;
            // case Hit::Cilindro3:
            //     cor = cilindro3.CalcularCor(cena, t_best, dr_e);
            //     break;
            // case Hit::Cilindro4:
            //     cor = cilindro4.CalcularCor(cena, t_best, dr_e);
            //     break;
            // case Hit::Cilindro5:
            //     cor = cilindro5.CalcularCor(cena, t_best, dr_e);
            //     break;
            // case Hit::Cilindro6:
            //     cor = cilindro6.CalcularCor(cena, t_best, dr_e);
            //     break;
            // case Hit::Cilindro7:
            //     cor = cilindro7.CalcularCor(cena, t_best, dr_e);
            //     break;
            // case Hit::Cilindro8:
            //     cor = cilindro8.CalcularCor(cena, t_best, dr_e);
            //     break;
            // case Hit::Cilindro9:
            //     cor = cilindro9.CalcularCor(cena, t_best, dr_e);
            //     break;
            // case Hit::Cauda:
            //     cor = cauda.CalcularCor(cena, t_best, dr_e);
            //     break;
            // case Hit::Cone:
            //     cor = cone.CalcularCor(cena, t_best, dr_e);
            //     break;
            // case Hit::Cone2:
            //     cor = cone2.CalcularCor(cena, t_best, dr_e);
            //     break;
            // case Hit::Cone3:
            //     cor = cone3.CalcularCor(cena, t_best, dr_e);
            //     break;
            // case Hit::Cone8:
            //     cor = cone8.CalcularCor(cena, t_best, dr_e);
            //     break;
            // case Hit::Cone9:
            //     cor = cone9.CalcularCor(cena, t_best, dr_e);
            //     break;
            // case Hit::Chifre_esq:
            //     cor = chifre_esq.CalcularCor(cena, t_best, dr_e);
            //     break;
            // case Hit::Chifre_dir:
            //     cor = chifre_dir.CalcularCor(cena, t_best, dr_e);
            //     break;
            // case Hit::Nave: {
            //     Point Pi = calcula_eq_ray(Po, t_best, dr_e);
            //     Vector dist_centro = subtrai_pontos(Pi, nave.cb);
            //     float altura_Pi = calcula_prod_esc(dist_centro, nave.dc);
            //     float dist2 = calcula_prod_esc(dist_centro, dist_centro);

            //     const float epsBase = 1e-3f;
            //     if (std::fabs(altura_Pi) < epsBase && dist2 <= nave.raio * nave.raio + epsBase) {
            //         cor = Color(45, 50, 60);
            //     }
            //     else {
            //         cor = nave.CalcularCor(cena, t_best, dr_e);
            //     }
            //     break;
            // }
            // case Hit::Esfera1_nave:
            //     cor = esfera1_nave.CalcularCor(cena, t_best, dr_e);
            //     break;
            // case Hit::Esfera2_nave:
            //     cor = esfera2_nave.CalcularCor(cena, t_best, dr_e);
            //     break;
            // case Hit::Esfera3_nave:
            //     cor = esfera3_nave.CalcularCor(cena, t_best, dr_e);
            //     break;
            // case Hit::Cubo: {
            //     Vetor PiV = cuboMalha.getPontoIntersecao();
            //     Vetor nV = cuboMalha.calcularNormal(PiV);

            //     // Cria um plano auxiliar (mesma face atingida) para reaproveitar o shading do Plano,
            //     // mas com um ponto de referência estável para o mapeamento de textura.
            //     // Para um plano com normal n e passando por Pi: d = n·Pi.
            //     // Escolhe um ponto p_pi no plano resolvendo uma coordenada (a de maior módulo).
            //     const Vector nFace(nV.i, nV.j, nV.k);
            //     const float d = nFace.i * PiV.i + nFace.j * PiV.j + nFace.k * PiV.k;
            //     Point p_pi(0.0f, 0.0f, 0.0f);
            //     const float ax = std::fabs(nFace.i);
            //     const float ay = std::fabs(nFace.j);
            //     const float az = std::fabs(nFace.k);
            //     if (az >= ax && az >= ay && std::fabs(nFace.k) > 1e-6f) {
            //         p_pi = Point(0.0f, 0.0f, d / nFace.k);
            //     }
            //     else if (ay >= ax && std::fabs(nFace.j) > 1e-6f) {
            //         p_pi = Point(0.0f, d / nFace.j, 0.0f);
            //     }
            //     else if (std::fabs(nFace.i) > 1e-6f) {
            //         p_pi = Point(d / nFace.i, 0.0f, 0.0f);
            //     }

            //     Material mat_vaca;
            //     mat_vaca.usarTextura = (texturaVaca != nullptr);
            //     mat_vaca.textura = texturaVaca;
            //     // Para o cubo, deixa a textura “1x por face” (pode ajustar depois)
            //     mat_vaca.texturaScale = 0.02f;

            //     Plano plano_cubo(p_pi, nFace, mat_vaca);
            //     plano_cubo.setKa(cuboMalha.getKa());
            //     plano_cubo.setKd(cuboMalha.getKd());
            //     plano_cubo.setKe(cuboMalha.getKe());
            //     plano_cubo.setShininess(cuboMalha.getShininess());
            //     cor = plano_cubo.CalcularCor(cena, dr_e);
            //     break;
            // }
            // case Hit::Esfera_cabeca:
            //     cor = esfera_cabeca.CalcularCor(cena, t_best, dr_e);
            //     break;
            // case Hit::Lua:
            //     cor = lua.CalcularCor(cena, t_best, dr_e);
            //     break;
            // case Hit::None:
            // default:
            //     if (texturaCeu != nullptr) {
            //         const std::size_t w = texturaCeu->get_largura_pixels();
            //         const std::size_t h = texturaCeu->get_altura_pixels();

            //         if (w > 0 && h > 0) {
            //             const float u = (nCol > 1) ? (static_cast<float>(coluna) / static_cast<float>(nCol - 1)) : 0.0f;
            //             const float v = (nLin > 1) ? (static_cast<float>(linha) / static_cast<float>(nLin - 1)) : 0.0f;

            //             std::size_t tx = static_cast<std::size_t>(u * static_cast<float>(w - 1));
            //             std::size_t ty = static_cast<std::size_t>(v * static_cast<float>(h - 1));
            //             if (tx >= w) tx = w - 1;
            //             if (ty >= h) ty = h - 1;

            //             // Observação: get_cor_pixel usa (linha, coluna).
            //             rgb px = texturaCeu->get_cor_pixel(ty, tx);
            //             cor = Color(static_cast<float>(px[0]), static_cast<float>(px[1]), static_cast<float>(px[2]));
            //         }
            //     }
            //     break;
            // }

            Point p_aux = calcula_eq_ray(Po, t_best, dr_e);

            switch (hit) {
            case Hit::Fundo:
                cor = plano_fundo.CalcularCor(cena, dr_e);
                break;
            case Hit::Chao:
                cor = plano_chao.CalcularCor(cena, dr_e);
                break;
            case Hit::Cilindro:
                // cor = cilindro.CalcularCor(cena, t_best, dr_e);
                cor = CalcularCor(cena, t_best, dr_e, cilindro.getKa(), 
                                    cilindro.getKd(), cilindro.getKe(), cilindro.CalcularNormal(p_aux), &cilindro);
                break;
            case Hit::Cilindro2:
                //cor = cilindro2.CalcularCor(cena, t_best, dr_e);
                cor = CalcularCor(cena, t_best, dr_e, cilindro2.getKa(),  
                                     cilindro2.getKd(), cilindro2.getKe(), cilindro2.CalcularNormal(p_aux), &cilindro2);
                break;
            case Hit::Cilindro3:
                // cor = cilindro3.CalcularCor(cena, t_best, dr_e);
                cor = CalcularCor(cena, t_best, dr_e, cilindro3.getKa(), 
                                     cilindro3.getKd(), cilindro3.getKe(), cilindro3.CalcularNormal(p_aux), &cilindro3);
                break;
            case Hit::Cilindro4:
                //cor = cilindro4.CalcularCor(cena, t_best, dr_e);
                cor = CalcularCor(cena, t_best, dr_e, cilindro4.getKa(),
                                cilindro.getKd(), cilindro4.getKe(), cilindro4.CalcularNormal(p_aux), &cilindro4);
                break;
            case Hit::Cilindro5:
                // cor = cilindro5.CalcularCor(cena, t_best, dr_e);
                cor = CalcularCor(cena, t_best, dr_e, cilindro5.getKa(), 
                                    cilindro5.getKd(), cilindro5.getKe(), cilindro5.CalcularNormal(p_aux), &cilindro5);
                break;
            case Hit::Cilindro6:
                // cor = cilindro6.CalcularCor(cena, t_best, dr_e);
                cor = CalcularCor(cena, t_best, dr_e, cilindro6.getKa(), 
                                     cilindro6.getKd(), cilindro6.getKe(), cilindro6.CalcularNormal(p_aux), &cilindro6);
                break;
            case Hit::Cilindro7:
                // cor = cilindro7.CalcularCor(cena, t_best, dr_e);
                cor = CalcularCor(cena, t_best, dr_e, cilindro7.getKa(), 
                                     cilindro7.getKd(), cilindro7.getKe(), cilindro7.CalcularNormal(p_aux), &cilindro7);
                break;
            case Hit::Cilindro8:
                //cor = cilindro8.CalcularCor(cena, t_best, dr_e);
                cor = CalcularCor(cena, t_best, dr_e, cilindro8.getKa(), 
                                     cilindro8.getKd(), cilindro8.getKe(), cilindro8.CalcularNormal(p_aux), &cilindro8);
                break;
            case Hit::Cilindro9:
                // cor = cilindro9.CalcularCor(cena, t_best, dr_e);
                cor = CalcularCor(cena, t_best, dr_e, cilindro9.getKa(), 
                                     cilindro9.getKd(), cilindro9.getKe(), cilindro9.CalcularNormal(p_aux), &cilindro9);
                break;
            case Hit::Cauda:
                // cor = cauda.CalcularCor(cena, t_best, dr_e);
                cor = CalcularCor(cena, t_best, dr_e, cauda.getKa(), 
                                     cauda.getKd(), cauda.getKe(), cauda.CalcularNormal(p_aux), &cauda);
                break;
            case Hit::Cone:
                // cor = cone.CalcularCor(cena, t_best, dr_e);
                cor = CalcularCor(cena, t_best, dr_e, cone.getKa(), 
                                    cone.getKd(), cone.getKe(), cone.CalcularNormal(p_aux, dr_e), &cone);
                break;
            case Hit::Cone2:
                //cor = cone2.CalcularCor(cena, t_best, dr_e);
                cor = CalcularCor(cena, t_best, dr_e, cone2.getKa(), 
                                     cone2.getKd(), cone2.getKe(), cone2.CalcularNormal(p_aux, dr_e), &cone2);
                break;
            case Hit::Cone3:
                // cor = cone3.CalcularCor(cena, t_best, dr_e);
                cor = CalcularCor(cena, t_best, dr_e, cone3.getKa(), 
                                     cone3.getKd(), cone3.getKe(), cone3.CalcularNormal(p_aux, dr_e), &cone3);
                break;
            case Hit::Cone8:
                //cor = cone8.CalcularCor(cena, t_best, dr_e);
                cor = CalcularCor(cena, t_best, dr_e, cone8.getKa(), 
                                     cone8.getKd(), cone8.getKe(), cone8.CalcularNormal(p_aux, dr_e), &cone8);
                break;
            case Hit::Cone9:
                //cor = cone9.CalcularCor(cena, t_best, dr_e);
                cor = CalcularCor(cena, t_best, dr_e, cone9.getKa(), 
                                     cone9.getKd(), cone9.getKe(), cone9.CalcularNormal(p_aux, dr_e), &cone9);
                break;
            case Hit::Chifre_esq:
                //cor = chifre_esq.CalcularCor(cena, t_best, dr_e);
                cor = CalcularCor(cena, t_best, dr_e, chifre_esq.getKa(),
                                        chifre_esq.getKd(), chifre_esq.getKe(), chifre_esq.CalcularNormal(p_aux, dr_e), &chifre_esq);    
                break;
            case Hit::Chifre_dir:
                //cor = chifre_dir.CalcularCor(cena, t_best, dr_e);
                cor = CalcularCor(cena, t_best, dr_e, chifre_dir.getKa(),
                                        chifre_dir.getKd(), chifre_dir.getKe(), chifre_dir.CalcularNormal(p_aux, dr_e), &chifre_dir);
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
                    // cor = nave.CalcularCor(cena, t_best, dr_e);
                    cor = CalcularCor(cena, t_best, dr_e, nave.getKa(),
                                            nave.getKd(), nave.getKe(), nave.CalcularNormal(p_aux, dr_e), &nave); 
                }
                break;
            }
            // case Hit::Esfera:
            //     //cor = esfera.CalcularCor(cena, t_best, dr_e);
            //    break;
            case Hit::Esfera1_nave:
                //cor = esfera1_nave.CalcularCor(cena, t_best, dr_e);
                cor = CalcularCor(cena, t_best, dr_e, esfera1_nave.getKa(), 
                                     esfera1_nave.getKd(), esfera1_nave.getKe(), esfera1_nave.CalcularNormal(p_aux), &esfera1_nave);
                break;
            case Hit::Esfera2_nave:
                //cor = esfera2_nave.CalcularCor(cena, t_best, dr_e);
                cor = CalcularCor(cena, t_best, dr_e, esfera2_nave.getKa(), 
                                     esfera2_nave.getKd(), esfera2_nave.getKe(), esfera2_nave.CalcularNormal(p_aux), &esfera2_nave);
                break;
            case Hit::Esfera3_nave:
                //cor = esfera3_nave.CalcularCor(cena, t_best, dr_e);
                cor = CalcularCor(cena, t_best, dr_e, esfera3_nave.getKa(), 
                                     esfera3_nave.getKd(), esfera3_nave.getKe(), esfera3_nave.CalcularNormal(p_aux), &esfera3_nave);
                break;
            // case Hit::Cubo: {
            //     // Vetor PiV = cuboMalha.getPontoIntersecao();
            //     // Vetor nV = cuboMalha.calcularNormal(PiV);
            //     // Plano plano_cubo(Point(PiV.i, PiV.j, PiV.k), Vector(nV.i, nV.j, nV.k), Material());
            //     // plano_cubo.setKa(cuboMalha.getKa());
            //     // plano_cubo.setKd(cuboMalha.getKd());
            //     // plano_cubo.setKe(cuboMalha.getKe());
            //     // plano_cubo.setShininess(cuboMalha.getShininess());
            //     //cor = plano_cubo.CalcularCor(cena, dr_e);
            //     cor = cuboMalha.CalcularCor(cena, t_best, dr_e);
            //     break;
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
                cuboMalha.plano = plano_cubo;
                // cor = plano_cubo.CalcularCor(cena, dr_e);
                cor = cuboMalha.CalcularCor(cena, t_best, dr_e);

                break;
            }
            case Hit::Esfera_cabeca:
                cor = esfera_cabeca.CalcularCor(cena, t_best, dr_e);
                break;
            case Hit::Lua:
                cor = lua.CalcularCor(cena, t_best, dr_e);
                break;
            case Hit::None:
            default:
                if (texturaCeu != nullptr) {
                    const std::size_t w = texturaCeu->get_largura_pixels();
                    const std::size_t h = texturaCeu->get_altura_pixels();

                    if (w > 0 && h > 0) {
                        const float u = (nCol > 1) ? (static_cast<float>(coluna) / static_cast<float>(nCol - 1)) : 0.0f;
                        const float v = (nLin > 1) ? (static_cast<float>(linha) / static_cast<float>(nLin - 1)) : 0.0f;

                        std::size_t tx = static_cast<std::size_t>(u * static_cast<float>(w - 1));
                        std::size_t ty = static_cast<std::size_t>(v * static_cast<float>(h - 1));
                        if (tx >= w) tx = w - 1;
                        if (ty >= h) ty = h - 1;

                        // Observação: get_cor_pixel usa (linha, coluna).
                        rgb px = texturaCeu->get_cor_pixel(ty, tx);
                        cor = Color(static_cast<float>(px[0]), static_cast<float>(px[1]), static_cast<float>(px[2]));
                    }
                }
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
