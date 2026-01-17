#include "../include/setup.h"

void Setup::paintPixel(SDL_Renderer *renderer, int x, int y){
    SDL_RenderDrawPoint(renderer, x, y);
}

void Setup::update(SDL_Renderer *renderer){
    SDL_RenderPresent(renderer);
}

void Setup::setPaintColor(SDL_Renderer *renderer, int r, int g, int b, int a){
    SDL_SetRenderDrawColor(renderer, r, g, b, a);
}

void Setup::setWindowBackground(SDL_Renderer *renderer, int r, int g, int b, int a){
    SDL_SetRenderDrawColor(renderer, r, g, b, a);
    SDL_RenderClear(renderer);
}

void Setup::listenEventQuit(SDL_Window *window){
    SDL_Event event;
    bool isRunning = true;

    while (isRunning)
    {
        while (SDL_PollEvent(&event) != 0)
        {
            if (event.type == SDL_QUIT)
                isRunning = false;
        }

        SDL_UpdateWindowSurface(window);
    }

    SDL_DestroyWindow(window);
    SDL_Quit();
}

void Setup::initializeSDLAndWindow(SDL_Window **window, SDL_Renderer **renderer, int width, int height){
    SDL_Init(SDL_INIT_VIDEO);

    SDL_CreateWindowAndRenderer(width, height, 0, window, renderer);
}


