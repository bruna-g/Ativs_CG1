#ifndef CUBO_H
#define CUBO_H

#include "Malha.h"
#include "Vetor.h"

class Cubo {
public:
    static Malha criarCubo(Vetor centroBase, double comprimentoAresta, Vetor Ka, Vetor Kd, Vetor Ke, double shininess);
};

#endif // CUBO_H
