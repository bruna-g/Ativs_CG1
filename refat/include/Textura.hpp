#ifndef UTILS_TEXTURA_HPP
#define UTILS_TEXTURA_HPP

#include "tipos.hpp"
#include <string>

class Textura {

    private:

        // Nome da textura
        std::string nome;

        // Altura da textura em pixels.
        std::size_t altura_pixels;
        // Largura da textura em pixels.
        std::size_t largura_pixels;

        // Matriz guardando as cores de cada pixel.
        matriz_cores pixels;

    public:

        // --- CONSTRUTORES ---

            Textura(std::string nome = std::string("madeira"), std::string arquivo = std::string("madeira.bmp"));
        
        // --- GETTERS E SETTERS ---

            std::string get_nome() const;

            std::size_t get_altura_pixels() const;

            std::size_t get_largura_pixels() const;

            rgb get_cor_pixel(std::size_t x, std::size_t y);

        // --- OUTROS MÉTODOS ---

            // Carrega uma imagem da textura no formato BMP.
            void carregar_imagem(const char* arquivo);

};

#endif