#include <iostream>
#include <vector>
#include <cmath>
#include <netpbm/pam.h>
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

int main(){
    float wJanela = 25.0;
    float hJanela = 25.0;
    int dJanela = 10;

    Point Po(0,0,0);

    float rEsfera = 1.0;
    Point centroEsfera(0,0,-(dJanela));
    int nCol = 500;
    int nLin = 500;


    float Dx = wJanela/nCol;
    float Dy = hJanela/nLin;

    float z = -dJanela;

    pixval maxval = 255;
    pm_init("esfera", 0);
    pixel **pixels = ppm_allocarray(nCol, nLin);
    
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

            if (delta >= 0) {
                PPM_ASSIGN(pixels[l][c], 255, 0, 0);
            } else {
                PPM_ASSIGN(pixels[l][c], 100,100,100);
            }
        }
    }
    FILE *fp = fopen("esfera.ppm", "wb");
    ppm_writeppm(fp, pixels, nCol, nLin, maxval, 0);
    fclose(fp);

    ppm_freearray(pixels, nLin);

    std::cout << "Imagem gerada: esfera.ppm\n";
    return 0;
}
