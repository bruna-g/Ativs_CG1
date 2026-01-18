#include "../include/Utils.h"
#include <cmath>

float calcula_norma(const Vector& v) {
    return std::sqrt(v.i * v.i + v.j * v.j + v.k * v.k);
}

float calcula_prod_esc(const Vector& v1, const Vector& v2) {
    return (v1.i * v2.i + v1.j * v2.j + v1.k * v2.k);
}

Vector calcula_esc_por_vetor(float v1, const Vector& v2) {
    return Vector(v1 * v2.i, v1 * v2.j, v1 * v2.k);
}

Vector subtrai_pontos(Point& p1, Point& p2) {
    Vector sub(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
    return sub;
}

Point soma_ponto_vetor(const Point& p1, const Vector& v2) {
    Point soma(p1.x + v2.i, p1.y + v2.j, p1.z + v2.q);
    return soma;
}

Vector subtrai_pontos(const Point& p1, const Point& p2) {
    Vector sub(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
    return sub;
}

Vector subtrai_vetores(Vector& v1, Vector& v2) {
    Vector sub(v1.i - v2.i, v1.j - v2.j, v1.k - v2.k);
    return sub;
}

Vector subtrai_vetores(const Vector& v1, const Vector& v2) {
    Vector sub(v1.i - v2.i, v1.j - v2.j, v1.k - v2.k);
    return sub;
}

Vector soma_vetores(const Vector& v1, const Vector& v2) {
    Vector sub(v1.i + v2.i, v1.j + v2.j, v1.k + v2.k);
    return sub;
}

Vector calcula_dr(Point& Po, Point& Pj) {
    Vector Dr = subtrai_pontos(Pj, Po);
    float drNorma = calcula_norma(Dr);
    Vector dr(Dr.i / drNorma, Dr.j / drNorma, Dr.k / drNorma);
    return dr;
}

Vector calcula_dr(const Point& Po, const Point& Pj) {
    Vector Dr = subtrai_pontos(Pj, Po);
    float drNorma = calcula_norma(Dr);
    Vector dr(Dr.i / drNorma, Dr.j / drNorma, Dr.k / drNorma);
    return dr;
}

Vector calcula_l(Point PF, Point Pi) {
    Vector PF_Pi = subtrai_pontos(PF, Pi);
    float norma = calcula_norma(PF_Pi);
    Vector l(PF_Pi.i / norma, PF_Pi.j / norma, PF_Pi.k / norma);
    return l;
}

Vector calcula_n(Point Pi, Point C, float R) {
    Vector Pi_C = subtrai_pontos(Pi, C);
    Vector n((Pi_C.i / R), (Pi_C.j / R), (Pi_C.k / R));
    return n;
}

Point calcula_eq_ray(Point Po, float t, Vector dr) {
    Vector t_dr(dr.i * t, dr.j * t, dr.k * t);
    Point Pi(Po.x + t_dr.i, Po.y + t_dr.j, Po.z + t_dr.k);
    return Pi;
}

float lidarExcecao(float v) {
    if (v < 0.0f) return 0.0f;
    if (v > 1.0f) return 1.0f;
    return v;
}

Vector normalizar(const Vector& v) {
    float mag = sqrt(v.i*v.i + v.j*v.j + v.k*v.k);
    return Vector(v.i/mag, v.j/mag, v.k/mag);
}

Vector produto_vetorial(const Vector& a, const Vector& b) {
    return Vector(
        a.j*b.k - a.k*b.j,
        a.k*b.i - a.i*b.k,
        a.i*b.j - a.j*b.i
    );
}
