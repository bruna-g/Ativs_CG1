#include "../include/Textura.hpp"
#include <stdexcept>
#include <string>
#include <SDL2/SDL.h>

Textura::Textura(std::string nome, std::string arquivo) {

    this->nome = nome;
    this->carregar_imagem(arquivo.c_str());

}

std::string Textura::get_nome() const {

    return this->nome;

}

std::size_t Textura::get_altura_pixels() const {

    return this->altura_pixels;

}

std::size_t Textura::get_largura_pixels() const {

    return this->largura_pixels;

}

rgb Textura::get_cor_pixel(std::size_t x, std::size_t y) {

    return this->pixels[x][y];

}

void Textura::carregar_imagem(const char* arquivo) {

    SDL_Surface* imagem;
    Uint32 cor_pixel_SDL;
    rgb cor_pixel;

    // Carregando as informações do arquivo BMP.
    imagem = SDL_LoadBMP(arquivo);

    // Checando se a imagem carregou.
    if (imagem != NULL) {

        // Limpando a imagem atual da textura.
        this->pixels.clear();

        // Guardando o tamanho da nova textura.
        this->altura_pixels = imagem->h;
        this->largura_pixels = imagem->w;

        for (std::size_t linha = 0; linha < this->altura_pixels; linha++) {

            // Inserindo linha de pixels.
            this->pixels.push_back(std::vector<rgb>());

            for (std::size_t coluna = 0; coluna < this->largura_pixels; coluna++) {

                // Conseguindo os valores RGB do pixel no formato do SDL.
                cor_pixel_SDL = *((Uint32*) ((Uint8*) imagem->pixels + linha * imagem->pitch + coluna * imagem->format->BytesPerPixel));

                // Guardandos valores RGB do pixel no meu formato.
                SDL_GetRGB(cor_pixel_SDL, imagem->format, &cor_pixel[0], &cor_pixel[1], &cor_pixel[2]);

                // Inserindo um pixel.
                this->pixels[linha].push_back(cor_pixel);

            }

        }

        SDL_FreeSurface(imagem);

    } else {

        throw std::runtime_error(std::string("Erro: Não foi possível carregar a imagem ") + arquivo);

    }

}