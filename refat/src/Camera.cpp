// filepath: /home/bruna/Documentos/ESTUDO/CG/Ativs_CG1/refat/src/Camera.cpp
#include "../include/Camera.h"
#include "../include/Utils.h"
#include "../include/Vector.h"
#include <cmath>

Camera::Camera() 
    : eye(0.f, 10.f, 0.f),
      at(0.f, 10.f, 30.f),
      up(0.f, 1.f, 0.f),
      distanciaFocal(30.0f),
      xMin(-30.0f),
      xMax(30.0f),
      yMin(-30.0f),
      yMax(30.0f)
{
    atualizarBase();
}

Camera::Camera(const Point& eye, const Point& at, const Vector& up,
               float d, float xMin, float xMax, float yMin, float yMax)
    : eye(eye),
      at(at),
      up(up),
      distanciaFocal(d),
      xMin(xMin),
      xMax(xMax),
      yMin(yMin),
      yMax(yMax)
{
    atualizarBase();
}

void Camera::atualizarBase() {
    // Calcular direção de visada: de eye para at
    direcao = Vector(
        at.x - eye.x,
        at.y - eye.y,
        at.z - eye.z
    );
    
    // Normalizar direção
    float lenDir = sqrt(direcao.i * direcao.i + direcao.j * direcao.j + direcao.k * direcao.k);
    if (lenDir > 1e-6f) {
        direcao.i /= lenDir;
        direcao.j /= lenDir;
        direcao.k /= lenDir;
    }
    
    // w aponta na direção oposta à visada (convenção right-handed)
    w = Vector(-direcao.i, -direcao.j, -direcao.k);
    
    // u é perpendicular a up e w (vetor lateral direito)
    u = produto_vetorial(up, w);
    float lenU = sqrt(u.i * u.i + u.j * u.j + u.k * u.k);
    if (lenU > 1e-6f) {
        u.i /= lenU;
        u.j /= lenU;
        u.k /= lenU;
    }
    
    // v é perpendicular a w e u (vetor "para cima" real da câmera)
    v = produto_vetorial(w, u);
    float lenV = sqrt(v.i * v.i + v.j * v.j + v.k * v.k);
    if (lenV > 1e-6f) {
        v.i /= lenV;
        v.j /= lenV;
        v.k /= lenV;
    }
}

void Camera::setEye(const Point& newEye) {
    eye = newEye;
    atualizarBase();
}

void Camera::setAt(const Point& newAt) {
    at = newAt;
    atualizarBase();
}

void Camera::setUp(const Vector& newUp) {
    up = newUp;
    atualizarBase();
}

void Camera::setDistanciaFocal(float d) {
    distanciaFocal = d;
}

void Camera::setCampoVisao(float newXMin, float newXMax, float newYMin, float newYMax) {
    xMin = newXMin;
    xMax = newXMax;
    yMin = newYMin;
    yMax = newYMax;
}

Vector Camera::gerarRaio(float x, float y) const {
    // x, y são coordenadas no plano de visualização (window coordinates)
    // Transformar para o sistema de coordenadas mundial usando a base ortonormal
    
    // O plano de visualização está a uma distância 'd' (distanciaFocal) do eye
    // na direção de visada
    
    // Ponto no plano de visualização em coordenadas de câmera: (x, y, -d)
    // Transformar para coordenadas mundiais:
    Vector raio(
        x * u.i + y * v.i - distanciaFocal * w.i,
        x * u.j + y * v.j - distanciaFocal * w.j,
        x * u.k + y * v.k - distanciaFocal * w.k
    );
    
    // Normalizar o raio
    float len = sqrt(raio.i * raio.i + raio.j * raio.j + raio.k * raio.k);
    if (len > 1e-6f) {
        raio.i /= len;
        raio.j /= len;
        raio.k /= len;
    }
    
    return raio;
}