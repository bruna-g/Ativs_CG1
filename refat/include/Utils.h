#ifndef UTILS_H
#define UTILS_H

#include "Point.h"
#include "Vector.h"
#include "Color.h"
#include "Cena.h"

float calcula_norma(const Vector& v);
float calcula_prod_esc(const Vector& v1, const Vector& v2);
Vector calcula_esc_por_vetor(float v1, const Vector& v2);
Vector subtrai_pontos(Point& p1, Point& p2);
Point soma_ponto_vetor(const Point& p1, const Vector& v2);
Vector subtrai_vetores(Vector& v1, Vector& v2);
Vector soma_vetores(const Vector& v1, const Vector& v2);
Vector calcula_dr(Point& Po, Point& Pj);
Vector calcula_l(Point PF, Point Pi);
Vector calcula_n(Point Pi, Point C, float R);
Point calcula_eq_ray(Point Po, float t, Vector dr);
float lidarExcecao(float v);
Vector normalizar(const Vector& v);
Vector produto_vetorial(const Vector& a, const Vector& b);

// Overloads const-friendly (mantém compatibilidade com as assinaturas antigas)
Vector subtrai_pontos(const Point& p1, const Point& p2);
Vector subtrai_vetores(const Vector& v1, const Vector& v2);
Vector calcula_dr(const Point& Po, const Point& Pj);
Color CalcularCor(const Cena& cena, float t, const Vector& dir, 
    const Color& K_a, const Color& K_d, const Color& K_e, Vector normal, Objeto* obj);

#endif // UTILS_H
