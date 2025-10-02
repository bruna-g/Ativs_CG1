#include <iostream>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;


class Point {
public:
    float x, y, z, p;
    Point(float x_p, float y_p, float z_p) {
        this->x = x_p;
        this->y = y_p;
        this->z = z_p;
        this->p = 1;
    }
};

class Vector {
public:
    float i, j, k, q;
    Vector(float i_v, float j_v, float k_v) {
        this->i = i_v;
        this->j = j_v;
        this->k = k_v;
        this->q = 0;
    }
};

class Color {
public:
    float r, g, b;
    Color(float r_c, float g_c, float b_c) {
        this->r = r_c;
        this->g = g_c;
        this->b = b_c;
    }
};

// Utilitário simples para limitar valores no intervalo [0,1]
static inline float lidarExcecao(float v) {
    if (v < 0.0f) return 0.0f;
    if (v > 1.0f) return 1.0f;
    return v;
}


Vector subtrai_pontos(Point& p1, Point& p2) {
    Vector sub(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
    return sub;
}

float calcula_norma(Vector& v) {
    return sqrt(v.i * v.i + v.j * v.j + v.k * v.k);
}

Vector calcula_dr(Point& Po, Point& Pj) {
    Vector Dr = subtrai_pontos(Pj, Po);
    float drNorma = calcula_norma(Dr);
    Vector dr(Dr.i / drNorma, Dr.j / drNorma, Dr.k / drNorma);
    return dr;
}

float calcula_prod_esc(const Vector& v1, const Vector& v2) {
    return (v1.i * v2.i + v1.j * v2.j + v1.k * v2.k);
}

Vector calcula_esc_por_vetor(float v1, const Vector& v2) {
    return Vector(v1 * v2.i, v1 * v2.j, v1 * v2.k);
}

Vector calcula_n(Point Pi, Point C, float R) {
    Vector Pi_C = subtrai_pontos(Pi, C);
    Vector n((Pi_C.i / R), (Pi_C.j / R), (Pi_C.k / R));
    return n;
}

Vector calcula_l(Point PF, Point Pi) {
    Vector PF_Pi = subtrai_pontos(PF, Pi);
    float norma = calcula_norma(PF_Pi);
    Vector l(PF_Pi.i / norma, PF_Pi.j / norma, PF_Pi.k / norma);
    return l;
}

Point calcula_eq_ray(Point Po, float t, Vector dr) {
    Vector t_dr(dr.i * t, dr.j * t, dr.k * t);
    // Corrige a equação do raio: P(t) = Po + t * d_r
    Point Pi(Po.x + t_dr.i, Po.y + t_dr.j, Po.z + t_dr.k);
    return Pi;
}

int main() {
    float wJanela = 25.0;
    float hJanela = 25.0;
    int dJanela = 10;

    Point Po(0, 0, 0);

    float rEsfera = 5.0;
    Point centroEsfera(0, 0, -(dJanela));
    int nCol = 500;
    int nLin = 500;


    float Dx = wJanela / nCol;
    float Dy = hJanela / nLin;

    float z = -dJanela;

    Color I_F(0.7f, 0.7f, 0.7f); // Intensidade da luz (branca)
    Point P_F(0, 5, 0);
    Color K(0.7f, 0.0f, 0.0f);  // Material com canal vermelho
    float m = 10.0f;            // Expoente especular

    std::ofstream img("esfera.ppm");
    img << "P3\n" << nCol << " " << nLin << "\n255\n";

    for (int l = 0; l < nLin; l++) {
        float y = hJanela / 2 - Dy / 2 - l * Dy;
        for (int c = 0; c < nCol; c++) {
            float x = -wJanela / 2 + Dx / 2 + c * Dx;
            Point Pj(x, y, z);

            Vector dr = calcula_dr(Po, Pj);

            Vector w = subtrai_pontos(Po, centroEsfera);

            float a_delta = calcula_prod_esc(dr, dr);
            float b_delta = 2.0f * calcula_prod_esc(dr, w);
            float c_delta = calcula_prod_esc(w, w) - rEsfera * rEsfera;

            float delta = b_delta * b_delta - 4.0f * a_delta * c_delta;       

            float t1 = (-b_delta + sqrt(delta)) / (2.0f * a_delta);
            float t2 = (-b_delta - sqrt(delta)) / (2.0f * a_delta);

            // Escolhe o menor t positivo
            float t = -1.0f;
            if (t1 > 0.0f && t2 > 0.0f) t = min(t1, t2);
            else if (t1 > 0.0f) t = t1;
            else if (t2 > 0.0f) t = t2;

            // Se ambos t negativos, pinta fundo e segue
            if (t <= 0.0f || delta < 0.0f) {
                img << "100 100 100 ";
                continue;
            }

            Point Pi = calcula_eq_ray(Po, t, dr);
            Vector n = calcula_n(Pi, centroEsfera, rEsfera);
            Vector X(-dr.i, -dr.j, -dr.k); // direção para o observador
            Vector l = calcula_l(P_F, Pi); // direção para a luz
            Vector r1 = calcula_esc_por_vetor(2.0f, calcula_esc_por_vetor(calcula_prod_esc(n, l), n));
            Vector r(r1.i - l.i, r1.j - l.j, r1.k - l.k); // reflexão

            float fd = lidarExcecao(calcula_prod_esc(n, l));
            float cosAlpha = lidarExcecao(calcula_prod_esc(r, X));
            float fe = pow(cosAlpha, m);

            // Componentes difusa e especular (apenas canal vermelho do material)
            Color Id(I_F.r * K.r * fd, I_F.g * K.g * fd, I_F.b * K.b * fd);
            Color Ie(I_F.r * K.r * fe, I_F.g * K.g * fe, I_F.b * K.b * fe);
            Color I(lidarExcecao(Id.r + Ie.r), lidarExcecao(Id.g + Ie.g), lidarExcecao(Id.b + Ie.b));

            int R = static_cast<int>(roundf(I.r * 255.0f));
            int G = static_cast<int>(roundf(I.g * 255.0f));
            int B = static_cast<int>(roundf(I.b * 255.0f));
            img << R << ' ' << G << ' ' << B << ' ';
        }
        img << "\n";
    }
    img.close();
    std::cout << "Imagem gerada: esfera.ppm\n";
    return 0;
}
