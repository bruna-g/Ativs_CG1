#pragma once

#include <limits>
#include <SDL2/SDL.h>

#include "Objeto.h"
#include "Point.h"
#include "Vector.h"
#include "tipos.hpp"

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
    Cauda,
    Cone,
    Cone2,
    Cone3,
    Chifre_esq,
    Chifre_dir,
    Cubo,
    Esfera,
    Nave,
    Esfera1_nave,
    Esfera2_nave,
    Esfera3_nave,
    Esfera_cabeca
};

struct PickResult {
    Hit hit = Hit::None;
    float t = std::numeric_limits<float>::infinity();
    Objeto* objeto = nullptr;
};

// Estado global atual de seleção (mesmo comportamento de quando estava no main.cpp).
extern PickResult gSelecionado;

const char* hitToString(Hit hit);

// Núcleo do picking
PickResult pickRay(const Point& Po, const Vector& dr_e);
PickResult pickFromScreen(int mouseX, int mouseY, int width, int height);

// Comandos da seleção
void imprimirAjudaSelecao();

void setCorObjeto(Objeto* obj, float r, float g, float b);
void aplicarMaterialCustom(Objeto* obj,
    float ka_r, float ka_g, float ka_b,
    float ke_r, float ke_g, float ke_b,
    float shininess);

void aplicarTranslacaoSelecionado(Hit hit, const Vetor& delta);
void aplicarRotacaoSelecionado(Hit hit, char eixo, float graus);
void aplicarEscalaSelecionado(Hit hit, float s);
void aplicarEscalaVetorSelecionado(Hit hit, const Vetor& escala);

// Helpers de entrada
bool ler3Floats(const char* prompt, float& a, float& b, float& c);
bool lerInt(const char* prompt, int& v);
bool lerCharEFloat(const char* prompt, char& c, float& v);
bool lerFloat(const char* prompt, float& v);

// Foco da janela SDL (útil após ler stdin)
void refocarJanela(SDL_Window* window);
