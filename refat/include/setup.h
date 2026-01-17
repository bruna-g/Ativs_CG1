#pragma once

#ifndef SETUP_H
#define SETUP_H

#include <iostream>
#include <SDL2/SDL.h>
#include <cmath>

// #include "cenario.h"
// #include "canvas.h"
// #include "janela.h"
// #include "observador.h"
// #include "vetor.h"
// #include "algebraLinear.h"
// #include "cilindro.h"
// #include "cone.h"
// #include "malha.h"
// #include "cubo.h"


class Setup {

public:

	void paintPixel(SDL_Renderer *renderer, int x, int y);

	void update(SDL_Renderer *renderer);

	void setPaintColor(SDL_Renderer *renderer, int r, int g, int b, int a);

	void setWindowBackground(SDL_Renderer *renderer, int r, int g, int b, int a);

	void listenEventQuit(SDL_Window *window);

	void initializeSDLAndWindow(SDL_Window **window, SDL_Renderer **renderer, int width, int height);

	//void pintarCanvas (SDL_Renderer *renderer, Canvas *canvas, Observador observador, Janela janela, Esfera* esfera, Vetor fundo, int d);

	//void renderizar(SDL_Renderer *renderer, Canvas canvas);
};

#endif // SETUP_H