#ifndef UTILS_H
#define UTILS_H

#include "Point.h"
#include "Vector.h"

float calcula_norma(const Vector& v);
float calcula_prod_esc(const Vector& v1, const Vector& v2);
Vector calcula_esc_por_vetor(float v1, const Vector& v2);
Vector subtrai_pontos(Point& p1, Point& p2);
Vector subtrai_vetores(Vector& v1, Vector& v2);
Vector calcula_dr(Point& Po, Point& Pj);
Vector calcula_l(Point PF, Point Pi);
Vector calcula_n(Point Pi, Point C, float R);
Point calcula_eq_ray(Point Po, float t, Vector dr);
float lidarExcecao(float v);

#endif // UTILS_H
