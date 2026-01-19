// filepath: /home/bruna/Documentos/ESTUDO/CG/Ativs_CG1/refat/include/Camera.h
#ifndef CAMERA_H
#define CAMERA_H

#include "Point.h"
#include "Vector.h"

class Camera {
private:
    // 2.1 - Especificações obrigatórias
    Point eye;              // 2.1.1 - Posição da câmera
    Point at;               // 2.1.2 - Ponto de visada (para onde a câmera olha)
    Vector up;              // 2.1.3 - Orientação da câmera (vetor "para cima")
    
    // 2.2 - Parâmetros adicionais
    float distanciaFocal;   // 2.2.1 - Distância focal (d)
    
    // 2.2.2 - Campo de visão (coordenadas da janela em coordenadas de câmera)
    float xMin;
    float xMax;
    float yMin;
    float yMax;
    
    // Base ortonormal (sistema de coordenadas da câmera)
    Vector u, v, w;         // u=right, v=up_real, w=-direction
    Vector direcao;         // Direção de visada normalizada
    
    void atualizarBase();

public:
    // Construtores
    Camera();
    Camera(const Point& eye, const Point& at, const Vector& up,
           float d, float xMin, float xMax, float yMin, float yMax);
    
    // Getters - Especificações da câmera
    Point getEye() const { return eye; }
    Point getAt() const { return at; }
    Vector getUp() const { return up; }
    float getDistanciaFocal() const { return distanciaFocal; }
    
    // Getters - Campo de visão
    float getXMin() const { return xMin; }
    float getXMax() const { return xMax; }
    float getYMin() const { return yMin; }
    float getYMax() const { return yMax; }
    
    // Getters - Base ortonormal
    Vector getU() const { return u; }
    Vector getV() const { return v; }
    Vector getW() const { return w; }
    Vector getDirecao() const { return direcao; }
    
    // Setters - Especificações da câmera
    void setEye(const Point& eye);
    void setAt(const Point& at);
    void setUp(const Vector& up);
    void setDistanciaFocal(float d);
    
    // Setters - Campo de visão
    void setCampoVisao(float xMin, float xMax, float yMin, float yMax);
    
    // Método principal: gerar raio para um pixel
    // x, y são coordenadas no plano de visualização
    Vector gerarRaio(float x, float y) const;
    
    // Calcular largura e altura da janela
    float getLarguraJanela() const { return xMax - xMin; }
    float getAlturaJanela() const { return yMax - yMin; }
};

#endif // CAMERA_H