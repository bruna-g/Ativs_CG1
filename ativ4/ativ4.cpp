#include <iostream>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
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

class Plano {
public:
    Point p_pi;
    Vector n;
    Plano(const Point& p_pi_p, const Vector& n_v) : p_pi(p_pi_p), n(n_v) {}
};

class Cilindro {
public:
    Point cb;
    float raio;
    float altura;
    Vector dc;
    Cilindro(const Point& cb_c, const float raio_c, const float altura_c, const Vector& dc_c): cb(cb_c),raio(raio_c),altura(altura_c),dc(dc_c) {}
};

float wJanela = 60.0f;
float hJanela = 60.0f;
int nCol = 500;
int nLin = 500;

float dJanela = 30.0f;
float z = -dJanela;

float Dx = wJanela / nCol;
float Dy = hJanela / nLin;

Point Po(0.f, 0.f, 0.f);

//esfera
float rEsfera = 40.0f;
Point centroEsfera(0.f, 0.f, -100.f);
Color K_e(0.7f, 0.2f, 0.2f);  // Material com canal vermelho
float m_e = 10.0f;            // Expoente especular

//chão
Point P_pi_chao(0.f, -(rEsfera), 0.f);
Vector n_chao(0.f, 1.0f, 0.f);
Plano plano_chao(P_pi_chao, n_chao);

Color KC_d(0.2f, 0.7f, 0.2f);
Color KC_e(0.0f, 0.0f, 0.0f);
float m_c = 1.0f;

//fundo
Point P_pi_fundo(0.f, 0.f, -200.f);
Vector n_fundo(0.f, 0.f, 1.0f);
Plano plano_fundo(P_pi_fundo, n_fundo);

Color KF_d(0.3f, 0.3f, 0.7f);
Color KF_e(0.0f, 0.0f, 0.0f);
float m_f = 1.0f;

//fonte luminosa
Color I_F(0.7f, 0.7f, 0.7f); // Intensidade da luz (branca)
Point P_F(0.f, 60.f, -30.f);    // Posição da fonte de luz

//luz ambiente
Color I_A(0.3f, 0.3f, 0.3f); // Intensidade da luz ambiente

//cilindro
float raio_cil = rEsfera/3;
float altura_cil = 3*rEsfera;
Vector dc(1/sqrt(3), 1/sqrt(3), -1/sqrt(3));
Cilindro cilindro(centroEsfera, raio_cil, altura_cil, dc);

Color KCil_d(0.2, 0.3, 0.8);
Color KCil_e(0.2, 0.3, 0.8);
Color KCil_a(0.2, 0.3, 0.8);


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

float calcula_norma(const Vector& v) {
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

float calcula_t_plano(Vector w, Vector n, Vector dr){
    float result = calcula_prod_esc(calcula_esc_por_vetor(-1, w), n)/calcula_prod_esc(dr,n);
    return result; 
}

Color calcula_Plano(Plano P, Vector dr, Color K_e, Color K_d){
    bool naSombra = false;
    
    // descobrir ponto Pi no plano a partir do raio que sai do olho do observador
    Vector w = subtrai_pontos(Po, P.p_pi);
    float ti = calcula_t_plano(w, P.n, dr);
    Point Pi = calcula_eq_ray(Po, ti, dr);

    // Verifica sombra: lança um raio de Pi em direção à fonte de luz e checa interseção com a esfera
    Vector l = calcula_l(P_F, Pi);
    // origem levemente deslocada para evitar auto-interseção
    Point Pi_mod(Pi.x + l.i * 1e-4f, Pi.y + l.j * 1e-4f, Pi.z + l.k * 1e-4f); 

    Vector w_sombra = subtrai_pontos(Pi_mod, centroEsfera);
    float a_delta = calcula_prod_esc(l, l);
    float b_delta = 2.0f * calcula_prod_esc(l, w_sombra);
    float c_delta = calcula_prod_esc(w_sombra, w_sombra) - rEsfera * rEsfera;

    float delta = b_delta * b_delta - 4.0f * a_delta * c_delta;
    
    if (delta > 0.f){
        float s1 = (-b_delta - sqrt(delta)) / (2.0f * a_delta);
        float s2 = (-b_delta + sqrt(delta)) / (2.0f * a_delta);
        // se houver interseção positiva antes da fonte (s entre 0 e distância até a luz), ponto está em sombra
        float dist_Pi_Luz = calcula_norma(subtrai_pontos(P_F, Pi_mod));
        if ((s1 > 1e-4f && s1 < dist_Pi_Luz) || (s2 > 1e-4f && s2 < dist_Pi_Luz)) naSombra = true;
    }
    
    Vector v = Vector(-dr.i, -dr.j, -dr.k);
    Vector r1 = calcula_esc_por_vetor(2.0f, calcula_esc_por_vetor(calcula_prod_esc(P.n, l), P.n));
    Vector r(r1.i - l.i, r1.j - l.j, r1.k - l.k);

    Color Ia(I_A.r * K_d.r, I_A.g * K_d.g, I_A.b * K_d.b);
    if (naSombra) {
        // Apenas componente ambiente
        Color I(lidarExcecao(Ia.r), lidarExcecao(Ia.g), lidarExcecao(Ia.b));
        int R = static_cast<int>(roundf(I.r * 255.0f));
        int G = static_cast<int>(roundf(I.g * 255.0f));
        int B = static_cast<int>(roundf(I.b * 255.0f));
        return Color(R, G, B);
    }

    float fd = lidarExcecao(calcula_prod_esc(P.n, l));
    float cosAlpha = lidarExcecao(calcula_prod_esc(r, v));
    float fe = pow(cosAlpha, m_c);
    Color Id(I_F.r * K_d.r * fd, I_F.g * K_d.g * fd, I_F.b * K_d.b * fd);
    Color Ie(I_F.r * K_e.r * fe, I_F.g * K_e.g * fe, I_F.b * K_e.b * fe);
    Color I(lidarExcecao(Id.r + Ie.r + Ia.r), lidarExcecao(Id.g + Ie.g + Ia.g), lidarExcecao(Id.b + Ie.b + Ia.b));
    int R = static_cast<int>(roundf(I.r * 255.0f));
    int G = static_cast<int>(roundf(I.g * 255.0f));
    int B = static_cast<int>(roundf(I.b * 255.0f));
    return Color(R, G, B);

}


int main() {
    std::ofstream img("esfera.ppm");
    img << "P3\n" << nCol << " " << nLin << "\n255\n";

    for (int l = 0; l < nLin; l++) {
        float y = hJanela / 2 - Dy / 2 - l * Dy;
        for (int c = 0; c < nCol; c++) {
            float x = -wJanela / 2 + Dx / 2 + c * Dx;
            Point Pj(x, y, z);

            Vector dr_e = calcula_dr(Po, Pj);

            Vector w_e = subtrai_pontos(Po, centroEsfera);

            float a_delta = calcula_prod_esc(dr_e, dr_e);
            float b_delta = 2.0f * calcula_prod_esc(dr_e, w_e);
            float c_delta = calcula_prod_esc(w_e, w_e) - rEsfera * rEsfera;

            float delta = b_delta * b_delta - 4.0f * a_delta * c_delta;

            float t1 = (-b_delta + sqrt(delta)) / (2.0f * a_delta);
            float t2 = (-b_delta - sqrt(delta)) / (2.0f * a_delta);

            // Escolhe o menor t positivo
            float t = -1.0f;
            if (t1 > 0.0f && t2 > 0.0f) t = min(t1, t2);
            else if (t1 > 0.0f) t = t1;
            else if (t2 > 0.0f) t = t2;

            // Verifica interseção com os planos
            Vector w_p_c = subtrai_pontos(Po, plano_chao.p_pi); // w do chão
            float ti_c = calcula_t_plano(w_p_c, plano_chao.n, dr_e); // ti do chão
            
            Vector w_p_f = subtrai_pontos(Po, plano_fundo.p_pi); // w do plano de fundo
            float ti_f = calcula_t_plano(w_p_f, plano_fundo.n, dr_e); // ti do plano de fundo

            if (ti_c > 0.0f && (ti_c < t || t < 0.0f)) t = ti_c;
            if (ti_f > 0.0f && (ti_f < t || t < 0.0f)) t = ti_f;

            //se interceptar fundo primeiro
            if (ti_f == t) {
                //Color cor_plano = calculaPlano_fundo(Pj);
                Color cor_plano = calcula_Plano(plano_fundo, dr_e, KF_e, KF_d);
                img << static_cast<int>(cor_plano.r) << ' ' << static_cast<int>(cor_plano.g) << ' ' << static_cast<int>(cor_plano.b) << ' ';
                continue;
            } //se interceptar o chão primeiro
            else if (ti_c == t) {
                //Color cor_plano = calculaPlano_chao(Pj);
                Color cor_plano = calcula_Plano(plano_chao, dr_e, KC_e, KC_d);
                img << static_cast<int>(cor_plano.r) << ' ' << static_cast<int>(cor_plano.g) << ' ' << static_cast<int>(cor_plano.b) << ' ';
                continue;
            }
            else if (t <= 0.0f || delta < 0.0f) {
                img << "100 100 100 ";
                continue;
            }

            // esfera
            Point Pi = calcula_eq_ray(Po, t, dr_e);
            Vector n = calcula_n(Pi, centroEsfera, rEsfera);
            Vector X(-dr_e.i, -dr_e.j, -dr_e.k);
            Vector l = calcula_l(P_F, Pi);
            Vector r1 = calcula_esc_por_vetor(2.0f, calcula_esc_por_vetor(calcula_prod_esc(n, l), n));
            Vector r(r1.i - l.i, r1.j - l.j, r1.k - l.k);

            float fd = lidarExcecao(calcula_prod_esc(n, l));
            float cosAlpha = lidarExcecao(calcula_prod_esc(r, X));
            float fe = pow(cosAlpha, m_e);

            Color Id(I_F.r * K_e.r * fd, I_F.g * K_e.g * fd, I_F.b * K_e.b * fd);
            Color Ie(I_F.r * K_e.r * fe, I_F.g * K_e.g * fe, I_F.b * K_e.b * fe);
            Color Ia(I_A.r * K_e.r, I_A.g * K_e.g, I_A.b * K_e.b);
            Color I(lidarExcecao(Id.r + Ie.r + Ia.r), lidarExcecao(Id.g + Ie.g + Ia.g), lidarExcecao(Id.b + Ie.b + Ia.b));

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