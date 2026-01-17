#ifndef MATERIAL_H
#define MATERIAL_H

class Textura;

struct Material {
    bool usarTextura;
    Textura* textura;
    float texturaScale;

    Material()
        : usarTextura(false),
        textura(nullptr),
        texturaScale(0.02f) {
    }
};

#endif // MATERIAL_H
