#ifndef MATERIAL_H
#define MATERIAL_H

#include "Color.h"

class Textura;

struct Material {
    Color Ka; // ambiente
    Color Kd; // difuso
    Color Ke; // especular
    float m;  // expoente especular

    bool usarTextura;
    Textura* textura;
    float texturaScale;

    Material()
        : Ka(0.0f, 0.0f, 0.0f),
        Kd(0.0f, 0.0f, 0.0f),
        Ke(0.0f, 0.0f, 0.0f),
        m(10.0f),
        usarTextura(false),
        textura(nullptr),
        texturaScale(0.02f) {
    }

    static Material Solido(const Color& K, float expoente = 10.0f) {
        Material mat;
        mat.Ka = K;
        mat.Kd = K;
        mat.Ke = K;
        mat.m = expoente;
        return mat;
    }
};

#endif // MATERIAL_H
