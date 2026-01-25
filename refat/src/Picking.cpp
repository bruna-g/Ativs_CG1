#include "Picking.h"

#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <limits>
#include <locale>

#include "Camera.h"
#include "Plano.h"
#include "Cilindro.h"
#include "Cone.h"
#include "Esfera.h"
#include "Malha.h"
#include "Utils.h"

// Dependências globais (continuam sendo definidas em main.cpp)
extern Camera camera;
extern Plano plano_chao;
extern Cilindro cilindro;
extern Cone cone;
extern Esfera esfera1_nave;
extern Esfera esfera2_nave;
extern Esfera esfera3_nave;
extern Cilindro cilindro2;
extern Cone cone2;
extern Cilindro cilindro3;
extern Cone cone3;
extern Cone cone8;
extern Cone cone9;
extern Cone chifre_esq;
extern Cone chifre_dir;
extern Esfera esfera_cabeca;
extern Cone nave;
extern Cilindro cilindro4;
extern Cilindro cilindro5;
extern Cilindro cilindro6;
extern Cilindro cilindro7;
extern Cilindro cilindro8;
extern Cilindro cilindro9;
extern Cilindro cauda;
extern Malha cuboMalha;

PickResult gSelecionado;

const char* hitToString(Hit hit) {
    switch (hit) {
    case Hit::None: return "None";
    case Hit::Fundo: return "Fundo";
    case Hit::Chao: return "Chao";
    case Hit::Esq: return "Esq";
    case Hit::Dir: return "Dir";
    case Hit::Teto: return "Teto";
    case Hit::Cilindro: return "Cilindro";
    case Hit::Cilindro2: return "Cilindro2";
    case Hit::Cilindro3: return "Cilindro3";
    case Hit::Cilindro4: return "Cilindro4";
    case Hit::Cilindro5: return "Cilindro5";
    case Hit::Cilindro6: return "Cilindro6";
    case Hit::Cilindro7: return "Cilindro7";
    case Hit::Cilindro8: return "Cilindro8";
    case Hit::Cilindro9: return "Cilindro9";
    case Hit::Cauda: return "Cauda";
    case Hit::Cone: return "Cone";
    case Hit::Cone2: return "Cone2";
    case Hit::Cone3: return "Cone3";
    case Hit::Cone8: return "Cone8";
    case Hit::Cone9: return "Cone9";
    case Hit::Chifre_esq: return "Chifre_esq";
    case Hit::Chifre_dir: return "Chifre_dir";
    case Hit::Cubo: return "Cubo";
    case Hit::Esfera: return "Esfera";
    case Hit::Nave: return "Nave";
    case Hit::Esfera1_nave: return "Esfera1_nave";
    case Hit::Esfera2_nave: return "Esfera2_nave";
    case Hit::Esfera3_nave: return "Esfera3_nave";
    case Hit::Esfera_cabeca: return "Esfera_cabeca";
    default: return "(desconhecido)";
    }
}

PickResult pickRay(const Point& Po, const Vector& dr_e) {
    PickResult result;

    auto considerar = [&](float t, Hit h, Objeto* obj) {
        if (t > 1e-4f && std::isfinite(t) && t < result.t) {
            result.t = t;
            result.hit = h;
            result.objeto = obj;
        }
        };

    float ti_c = plano_chao.CalcularIntersecao(Po, dr_e);

    float t_cil = cilindro.CalcularIntersecao(Po, dr_e);
    float ti_cone = cone.CalcularIntersecao(Po, dr_e);

    float ti_esf1_nave = esfera1_nave.CalcularIntersecao(Po, dr_e);
    float ti_esf2_nave = esfera2_nave.CalcularIntersecao(Po, dr_e);
    float ti_esf3_nave = esfera3_nave.CalcularIntersecao(Po, dr_e);

    float t_cil2 = cilindro2.CalcularIntersecao(Po, dr_e);
    float ti_cone2 = cone2.CalcularIntersecao(Po, dr_e);

    float t_cil3 = cilindro3.CalcularIntersecao(Po, dr_e);
    float ti_cone3 = cone3.CalcularIntersecao(Po, dr_e);
    float t_cil8 = cilindro8.CalcularIntersecao(Po, dr_e);
    float ti_cone8 = cone8.CalcularIntersecao(Po, dr_e);
    float t_cil9 = cilindro9.CalcularIntersecao(Po, dr_e);
    float ti_cone9 = cone9.CalcularIntersecao(Po, dr_e);

    float ti_chifre_esq = chifre_esq.CalcularIntersecao(Po, dr_e);
    float ti_chifre_dir = chifre_dir.CalcularIntersecao(Po, dr_e);

    float t_cubo = std::numeric_limits<float>::infinity();
    if (cuboMalha.verificarIntersecao(
        Vetor(Po.x, Po.y, Po.z, 0.0f),
        Vetor(dr_e.i, dr_e.j, dr_e.k, 0.0f))) {
        t_cubo = static_cast<float>(cuboMalha.getDistancia());
    }

    float ti_cabeca = esfera_cabeca.CalcularIntersecao(Po, dr_e);
    float ti_nave = nave.CalcularIntersecao(Po, dr_e);

    float t_cil4 = cilindro4.CalcularIntersecao(Po, dr_e);
    float t_cil5 = cilindro5.CalcularIntersecao(Po, dr_e);
    float t_cil6 = cilindro6.CalcularIntersecao(Po, dr_e);
    float t_cil7 = cilindro7.CalcularIntersecao(Po, dr_e);
    float t_cauda = cauda.CalcularIntersecao(Po, dr_e);

    considerar(ti_c, Hit::Chao, &plano_chao);
    considerar(t_cil, Hit::Cilindro, &cilindro);
    considerar(t_cil2, Hit::Cilindro2, &cilindro2);
    considerar(t_cil3, Hit::Cilindro3, &cilindro3);
    considerar(t_cil4, Hit::Cilindro4, &cilindro4);
    considerar(t_cil5, Hit::Cilindro5, &cilindro5);
    considerar(t_cil6, Hit::Cilindro6, &cilindro6);
    considerar(t_cil7, Hit::Cilindro7, &cilindro7);
    considerar(t_cauda, Hit::Cauda, &cauda);
    considerar(ti_cone, Hit::Cone, &cone);
    considerar(ti_cone2, Hit::Cone2, &cone2);
    considerar(ti_cone3, Hit::Cone3, &cone3);
    considerar(ti_cone8, Hit::Cone8, &cone8);
    considerar(ti_cone9, Hit::Cone9, &cone9);
    considerar(t_cil8, Hit::Cilindro8, &cilindro8);
    considerar(t_cil9, Hit::Cilindro9, &cilindro9);
    considerar(ti_chifre_esq, Hit::Chifre_esq, &chifre_esq);
    considerar(ti_chifre_dir, Hit::Chifre_dir, &chifre_dir);
    considerar(ti_nave, Hit::Nave, &nave);
    considerar(ti_esf1_nave, Hit::Esfera1_nave, &esfera1_nave);
    considerar(ti_esf2_nave, Hit::Esfera2_nave, &esfera2_nave);
    considerar(ti_esf3_nave, Hit::Esfera3_nave, &esfera3_nave);
    considerar(t_cubo, Hit::Cubo, &cuboMalha);
    considerar(ti_cabeca, Hit::Esfera_cabeca, &esfera_cabeca);

    return result;
}

PickResult pickFromScreen(int mouseX, int mouseY, int width, int height) {
    const float localDx = camera.getLarguraJanela() / static_cast<float>(width);
    const float localDy = camera.getAlturaJanela() / static_cast<float>(height);

    const float x = camera.getXMin() + (static_cast<float>(mouseX) + 0.5f) * localDx;
    const float y = camera.getYMax() - (static_cast<float>(mouseY) + 0.5f) * localDy;

    const Vector dr_e = camera.gerarRaio(x, y);
    const Point Po = camera.getEye();
    return pickRay(Po, dr_e);
}

static float normalizarCor(float v) {
    if (v > 1.0f) return v / 255.0f;
    return v;
}

void setCorObjeto(Objeto* obj, float r, float g, float b) {
    if (obj == nullptr) return;
    r = normalizarCor(r);
    g = normalizarCor(g);
    b = normalizarCor(b);

    obj->setKa(Vetor(0.2f * r, 0.2f * g, 0.2f * b));
    obj->setKd(Vetor(r, g, b));
    obj->setKe(Vetor(0.0f, 0.0f, 0.0f));
}

void aplicarMaterialCustom(Objeto* obj,
    float ka_r, float ka_g, float ka_b,
    float ke_r, float ke_g, float ke_b,
    float shininess) {
    if (obj == nullptr) return;
    if (!(shininess > 0.0f)) return;

    ka_r = normalizarCor(ka_r);
    ka_g = normalizarCor(ka_g);
    ka_b = normalizarCor(ka_b);
    ke_r = normalizarCor(ke_r);
    ke_g = normalizarCor(ke_g);
    ke_b = normalizarCor(ke_b);

    obj->setKa(Vetor(ka_r, ka_g, ka_b));
    obj->setKe(Vetor(ke_r, ke_g, ke_b));
    obj->setShininess(shininess);
}

void aplicarTranslacaoSelecionado(Hit hit, const Vetor& delta) {
    switch (hit) {
    case Hit::Chao:
        plano_chao.p_pi = plano_chao.aplicarTranslacao(plano_chao.p_pi, delta);
        break;
    case Hit::Esfera1_nave:
        esfera1_nave.centro = esfera1_nave.aplicarTranslacao(esfera1_nave.centro, delta);
        break;
    case Hit::Esfera2_nave:
        esfera2_nave.centro = esfera2_nave.aplicarTranslacao(esfera2_nave.centro, delta);
        break;
    case Hit::Esfera3_nave:
        esfera3_nave.centro = esfera3_nave.aplicarTranslacao(esfera3_nave.centro, delta);
        break;
    case Hit::Esfera_cabeca:
        esfera_cabeca.centro = esfera_cabeca.aplicarTranslacao(esfera_cabeca.centro, delta);
        break;
    case Hit::Cilindro:
        cilindro.cb = cilindro.aplicarTranslacao(cilindro.cb, delta);
        break;
    case Hit::Cilindro2:
        cilindro2.cb = cilindro2.aplicarTranslacao(cilindro2.cb, delta);
        break;
    case Hit::Cilindro3:
        cilindro3.cb = cilindro3.aplicarTranslacao(cilindro3.cb, delta);
        break;
    case Hit::Cilindro4:
        cilindro4.cb = cilindro4.aplicarTranslacao(cilindro4.cb, delta);
        break;
    case Hit::Cilindro5:
        cilindro5.cb = cilindro5.aplicarTranslacao(cilindro5.cb, delta);
        break;
    case Hit::Cilindro6:
        cilindro6.cb = cilindro6.aplicarTranslacao(cilindro6.cb, delta);
        break;
    case Hit::Cilindro7:
        cilindro7.cb = cilindro7.aplicarTranslacao(cilindro7.cb, delta);
        break;
    case Hit::Cilindro8:
        cilindro8.cb = cilindro8.aplicarTranslacao(cilindro8.cb, delta);
        break;
    case Hit::Cilindro9:
        cilindro9.cb = cilindro9.aplicarTranslacao(cilindro9.cb, delta);
        break;
    case Hit::Cauda:
        cauda.cb = cauda.aplicarTranslacao(cauda.cb, delta);
        break;
    case Hit::Cone:
        cone.cb = cone.aplicarTranslacao(cone.cb, delta);
        cone.v = cone.aplicarTranslacao(cone.v, delta);
        cone.recalcularDerivados();
        break;
    case Hit::Cone2:
        cone2.cb = cone2.aplicarTranslacao(cone2.cb, delta);
        cone2.v = cone2.aplicarTranslacao(cone2.v, delta);
        cone2.recalcularDerivados();
        break;
    case Hit::Cone3:
        cone3.cb = cone3.aplicarTranslacao(cone3.cb, delta);
        cone3.v = cone3.aplicarTranslacao(cone3.v, delta);
        cone3.recalcularDerivados();
        break;
    case Hit::Cone8:
        cone8.cb = cone8.aplicarTranslacao(cone8.cb, delta);
        cone8.v = cone8.aplicarTranslacao(cone8.v, delta);
        cone8.recalcularDerivados();
        break;
    case Hit::Cone9:
        cone9.cb = cone9.aplicarTranslacao(cone9.cb, delta);
        cone9.v = cone9.aplicarTranslacao(cone9.v, delta);
        cone9.recalcularDerivados();
        break;
    case Hit::Chifre_esq:
        chifre_esq.cb = chifre_esq.aplicarTranslacao(chifre_esq.cb, delta);
        chifre_esq.v = chifre_esq.aplicarTranslacao(chifre_esq.v, delta);
        chifre_esq.recalcularDerivados();
        break;
    case Hit::Chifre_dir:
        chifre_dir.cb = chifre_dir.aplicarTranslacao(chifre_dir.cb, delta);
        chifre_dir.v = chifre_dir.aplicarTranslacao(chifre_dir.v, delta);
        chifre_dir.recalcularDerivados();
        break;
    case Hit::Nave:
        nave.cb = nave.aplicarTranslacao(nave.cb, delta);
        nave.v = nave.aplicarTranslacao(nave.v, delta);
        nave.recalcularDerivados();
        break;
    case Hit::Cubo: {
        for (auto& vert : cuboMalha.vertices) {
            Point p(vert.pos.i, vert.pos.j, vert.pos.k);
            p = cuboMalha.aplicarTranslacao(p, delta);
            vert.pos = Vetor(p.x, p.y, p.z, vert.pos.q);
        }
        break;
    }
    default:
        break;
    }
}

void aplicarRotacaoSelecionado(Hit hit, char eixo, float graus) {
    switch (hit) {
    case Hit::Cilindro:
    case Hit::Cilindro2:
    case Hit::Cilindro3:
    case Hit::Cilindro4:
    case Hit::Cilindro5:
    case Hit::Cilindro6:
    case Hit::Cilindro7:
    case Hit::Cilindro8:
    case Hit::Cilindro9:
    case Hit::Cauda: {
        Cilindro* c = nullptr;
        if (hit == Hit::Cilindro) c = &cilindro;
        else if (hit == Hit::Cilindro2) c = &cilindro2;
        else if (hit == Hit::Cilindro3) c = &cilindro3;
        else if (hit == Hit::Cilindro4) c = &cilindro4;
        else if (hit == Hit::Cilindro5) c = &cilindro5;
        else if (hit == Hit::Cilindro6) c = &cilindro6;
        else if (hit == Hit::Cilindro7) c = &cilindro7;
        else if (hit == Hit::Cilindro8) c = &cilindro8;
        else if (hit == Hit::Cilindro9) c = &cilindro9;
        else if (hit == Hit::Cauda) c = &cauda;
        if (c == nullptr) break;

        if (eixo == 'x') c->rotacionarX(graus);
        else if (eixo == 'y') c->rotacionarY(graus);
        else c->rotacionarZ(graus);
        break;
    }
    case Hit::Cone:
    case Hit::Cone2:
    case Hit::Cone3:
    case Hit::Cone8:
    case Hit::Cone9:
    case Hit::Chifre_esq:
    case Hit::Chifre_dir:
    case Hit::Nave: {
        Cone* c = nullptr;
        if (hit == Hit::Cone) c = &cone;
        else if (hit == Hit::Cone2) c = &cone2;
        else if (hit == Hit::Cone3) c = &cone3;
        else if (hit == Hit::Cone8) c = &cone8;
        else if (hit == Hit::Cone9) c = &cone9;
        else if (hit == Hit::Chifre_esq) c = &chifre_esq;
        else if (hit == Hit::Chifre_dir) c = &chifre_dir;
        else if (hit == Hit::Nave) c = &nave;
        if (c == nullptr) break;

        if (eixo == 'x') c->rotacionarX(graus);
        else if (eixo == 'y') c->rotacionarY(graus);
        else c->rotacionarZ(graus);
        break;
    }
    case Hit::Cubo:
        if (eixo == 'x') cuboMalha.rotacionarX(graus);
        else if (eixo == 'y') cuboMalha.rotacionarY(graus);
        else cuboMalha.rotacionarZ(graus);
        break;
    default:
        break;
    }
}

void aplicarEscalaSelecionado(Hit hit, float s) {
    if (s <= 0.0f) return;
    aplicarEscalaVetorSelecionado(hit, Vetor(s, s, s, 0.0f));
}

void aplicarEscalaVetorSelecionado(Hit hit, const Vetor& escala) {
    if (escala.i <= 0.0f || escala.j <= 0.0f || escala.k <= 0.0f) return;

    switch (hit) {
    case Hit::Chao:
        if (plano_chao.material.usarTextura && plano_chao.material.textura != nullptr) {
            Vector nLocal = plano_chao.n;
            Vector arbitrary = (std::fabs(nLocal.i) < 0.9f) ? Vector(1.0f, 0.0f, 0.0f) : Vector(0.0f, 1.0f, 0.0f);
            Vector u_axis = cross(arbitrary, nLocal);
            float nu = calcula_norma(u_axis);
            if (nu <= 1e-6f) nu = 1.0f;
            u_axis = calcula_esc_por_vetor(1.0f / nu, u_axis);
            Vector v_axis = cross(nLocal, u_axis);

            const float sx = escala.i;
            const float sy = escala.j;
            const float sz = escala.k;
            auto lenEsc = [&](const Vector& a) {
                Vector svec(sx * a.i, sy * a.j, sz * a.k);
                return calcula_norma(svec);
                };

            const float fu = lenEsc(u_axis);
            const float fv = lenEsc(v_axis);
            const float fator = 0.5f * (fu + fv);
            if (fator > 1e-6f) {
                plano_chao.material.texturaScale /= fator;
            }
        }
        plano_chao.aplicarEscalaNoPivoObjeto(escala, plano_chao.p_pi);
        break;
    case Hit::Esfera1_nave:
        esfera1_nave.aplicarEscalaNoPivoObjeto(escala, esfera1_nave.centro);
        break;
    case Hit::Esfera2_nave:
        esfera2_nave.aplicarEscalaNoPivoObjeto(escala, esfera2_nave.centro);
        break;
    case Hit::Esfera3_nave:
        esfera3_nave.aplicarEscalaNoPivoObjeto(escala, esfera3_nave.centro);
        break;
    case Hit::Esfera_cabeca:
        esfera_cabeca.aplicarEscalaNoPivoObjeto(escala, esfera_cabeca.centro);
        break;
    case Hit::Cilindro:
        cilindro.aplicarEscalaNoPivoObjeto(escala, cilindro.cb);
        break;
    case Hit::Cilindro2:
        cilindro2.aplicarEscalaNoPivoObjeto(escala, cilindro2.cb);
        break;
    case Hit::Cilindro3:
        cilindro3.aplicarEscalaNoPivoObjeto(escala, cilindro3.cb);
        break;
    case Hit::Cilindro4:
        cilindro4.aplicarEscalaNoPivoObjeto(escala, cilindro4.cb);
        break;
    case Hit::Cilindro5:
        cilindro5.aplicarEscalaNoPivoObjeto(escala, cilindro5.cb);
        break;
    case Hit::Cilindro6:
        cilindro6.aplicarEscalaNoPivoObjeto(escala, cilindro6.cb);
        break;
    case Hit::Cilindro7:
        cilindro7.aplicarEscalaNoPivoObjeto(escala, cilindro7.cb);
        break;
    case Hit::Cilindro8:
        cilindro8.aplicarEscalaNoPivoObjeto(escala, cilindro8.cb);
        break;
    case Hit::Cilindro9:
        cilindro9.aplicarEscalaNoPivoObjeto(escala, cilindro9.cb);
        break;
    case Hit::Cauda:
        cauda.aplicarEscalaNoPivoObjeto(escala, cauda.cb);
        break;
    case Hit::Cone:
        cone.aplicarEscalaNoPivoObjeto(escala, cone.cb);
        break;
    case Hit::Cone2:
        cone2.aplicarEscalaNoPivoObjeto(escala, cone2.cb);
        break;
    case Hit::Cone3:
        cone3.aplicarEscalaNoPivoObjeto(escala, cone3.cb);
        break;
    case Hit::Cone8:
        cone8.aplicarEscalaNoPivoObjeto(escala, cone8.cb);
        break;
    case Hit::Cone9:
        cone9.aplicarEscalaNoPivoObjeto(escala, cone9.cb);
        break;
    case Hit::Chifre_esq:
        chifre_esq.aplicarEscalaNoPivoObjeto(escala, chifre_esq.cb);
        break;
    case Hit::Chifre_dir:
        chifre_dir.aplicarEscalaNoPivoObjeto(escala, chifre_dir.cb);
        break;
    case Hit::Nave:
        nave.aplicarEscalaNoPivoObjeto(escala, nave.cb);
        break;
    case Hit::Cubo: {
        Point pivo = cuboMalha.calcularCentro();
        cuboMalha.aplicarEscalaNoPivoObjeto(escala, pivo);
        break;
    }
    default:
        break;
    }
}

void imprimirAjudaSelecao() {
    std::cout << "\nComandos (com um objeto selecionado):\n"
        << "  c = mudar cor (r g b)\n"
        << "  m = mudar material (Ka_r Ka_g Ka_b) (Ke_r Ke_g Ke_b) (shininess)\n"
        << "  t = translacao (dx dy dz)\n"
        << "  r = rotacao (eixo x/y/z e angulo em graus)\n"
        << "  s = escala uniforme (fator)\n"
        << "  v = escala vetorial (sx sy sz)\n"
        << "  h = ajuda\n\n";
}

static bool lerLinha(std::string& line) {
    if (!std::getline(std::cin >> std::ws, line)) {
        return false;
    }
    return true;
}

bool ler3Floats(const char* prompt, float& a, float& b, float& c) {
    std::cout << prompt << std::flush;
    std::string line;
    if (!lerLinha(line)) return false;

    // Aceita tanto '.' quanto ',' como separador decimal.
    for (char& ch : line) {
        if (ch == ',') ch = '.';
    }

    std::istringstream iss(line);
    iss.imbue(std::locale::classic());
    return static_cast<bool>(iss >> a >> b >> c);
}

bool lerInt(const char* prompt, int& v) {
    std::cout << prompt << std::flush;
    std::string line;
    if (!lerLinha(line)) return false;
    std::istringstream iss(line);
    return static_cast<bool>(iss >> v);
}

bool lerCharEFloat(const char* prompt, char& c, float& v) {
    std::cout << prompt << std::flush;
    std::string line;
    if (!lerLinha(line)) return false;

    for (char& ch : line) {
        if (ch == ',') ch = '.';
    }

    std::istringstream iss(line);
    iss.imbue(std::locale::classic());
    return static_cast<bool>(iss >> c >> v);
}

bool lerFloat(const char* prompt, float& v) {
    std::cout << prompt << std::flush;
    std::string line;
    if (!lerLinha(line)) return false;

    for (char& ch : line) {
        if (ch == ',') ch = '.';
    }

    std::istringstream iss(line);
    iss.imbue(std::locale::classic());
    return static_cast<bool>(iss >> v);
}

void refocarJanela(SDL_Window* window) {
    if (window == nullptr) return;
    SDL_RaiseWindow(window);
}
