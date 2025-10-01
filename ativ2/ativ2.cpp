#include <iostream>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;


class Point{
    public:
    float x, y, z, p;
    Point(float x_p, float y_p, float z_p){
        this->x = x_p;
        this->y = y_p;
        this->z = z_p;
        this->p = 1;
    }
};

class Vector{
    public:
    float i, j, k, q;
    Vector(float i_v, float j_v, float k_v){
        this->i = i_v;
        this->j = j_v;
        this->k = k_v;
        this->q = 0;
    }
};

class Color{
    public:
    float r, g, b;
    Color(float r_c, float g_c, float b_c){
        this->r = r_c;
        this->g = g_c;
        this->b = b_c;
    }
};


Vector subtrai_pontos(Point&  p1, Point&  p2){
    Vector sub(p1.x- p2.x, p1.y - p2.y, p1.z - p2.z);
    return sub;
}

float calcula_norma(Vector& v){
    return sqrt(v.i*v.i + v.j*v.j + v.k*v.k);
}

Vector calcula_dr(Point& Po, Point& Pj){
    Vector Dr = subtrai_pontos(Pj, Po);
    float drNorma = calcula_norma(Dr);
    Vector dr(Dr.i/drNorma, Dr.j/drNorma, Dr.k/drNorma);
    return dr;
}

float calcula_prod_esc(Vector& v1, Vector& v2){
    return (v1.i*v2.i + v1.j*v2.j + v1.k*v2.k);
}

Vector calcula_n(Point Pi, Point C, float R){
    Vector Pi_C = subtrai_pontos(Pi, C);
    Vector n((Pi_C.i/R), (Pi_C.j/R), (Pi_C.k/R));
    return n;
}

Vector calcula_l(Point PF, Point Pi){
    Vector PF_Pi = subtrai_pontos(PF, Pi);
    float norma = calcula_norma(PF_Pi);
    Vector l(PF_Pi.i/norma, PF_Pi.j/norma, PF_Pi.k/norma);
    return l;
}

Point calcula_eq_ray(Point Po, float t, Vector dr){
    Vector t_dr(dr.i*t,dr.j*t,dr.k*t);
    Point Pi(Po.p+t_dr.i, Po.x+t_dr.j, Po.y+t_dr.k);
    return Pi;
}

int main(){
    float wJanela = 25.0;
    float hJanela = 25.0;
    int dJanela = 10;

    Point Po(0,0,0);

    float rEsfera = 5.0;
    Point centroEsfera(0,0,-(dJanela));
    int nCol = 500;
    int nLin = 500;


    float Dx = wJanela/nCol;
    float Dy = hJanela/nLin;

    float z = -dJanela;

    Color I_F(0.7,0.7,0.7);
    Point P_F(0,5,0);
    Color K(1,1,1);
    float m = 1;

    std::ofstream img("esfera.ppm");
    img << "P3\n" << nCol << " " << nLin << "\n255\n";
    
    for (int l=0; l<nLin; l++){ 
        float y = hJanela/2 -Dy/2 -l*Dy; 
        for (int c=0; c<nCol; c++){
            float x = -wJanela/2 +Dx/2 +c*Dx;
            Point Pj(x,y,z);

            Vector dr = calcula_dr(Po, Pj);

            Vector w = subtrai_pontos(Po, centroEsfera);

            float a_delta = calcula_prod_esc(dr,dr);
            float b_delta = 2 * calcula_prod_esc(dr, w);
            float c_delta = calcula_prod_esc(w,w) - rEsfera*rEsfera;

            float delta = b_delta*b_delta - 4*a_delta*c_delta;

            float t1 = (-b_delta + sqrt(delta))/(2*a_delta);
            float t2 = (-b_delta - sqrt(delta))/(2*a_delta);
            float t = min(t1, t2);

            Point Pi = calcula_eq_ray(Po, t, dr);
            Vector n = calcula_n(Pi, centroEsfera, rEsfera);
            Vector x(-dr.i, -dr.j, -dr.k);


            if (delta >= 0) {
                img << "255 0 0 ";
            } else {
                img << "100 100 100 ";
            }
        }
        img << "\n";
    }
    img.close();
    std::cout << "Imagem gerada: esfera.ppm\n";
    return 0;
}
