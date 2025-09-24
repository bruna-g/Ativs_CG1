#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

int main(){
    float wJanela = 5.0;
    float hJanela = 5.0;
    float dJanela = 2.0;

    Point Po(0,0,0);

    float rEsfera = 3.0;
    Point centroEsfera(0,0,-(dJanela + rEsfera));
    int nCol = 15;
    int nLin = 15;

    vector<vector<int>> matCanvas(nLin, vector<int>(nCol,0));

    int esfColor[3] = {255, 0, 0};
    int bgColor[3] = {100,100,100};

    float Dx = wJanela/nCol;
    float Dy = hJanela/nLin;

    int z = -dJanela;
    for (int l=0; l<nLin; l++){ 
        int y = -hJanela/2 -Dy/2 -l*Dx; 
        for (int c=0; c<nCol; c++){
            int x = -wJanela/2 +Dx/2 +c*Dx;
            Point Pj(x,y,z);

            Vector dr = calcula_dr(Po, Pj);

            Vector w = subtrai_pontos(Po, centroEsfera);

            int a = 

            

        }
    }
    return 0;
}

class Point{
    public:
    int x, y, z, p;
    Point(int x_p, int y_p, int z_p){
        x = x_p;
        y = y_p;
        z = z_p;
        p = 1;
    }
};

class Vector{
    public:
    int i, j, k, q;
    Vector(int i_v, int j_v, int k_v){
        i = i_v;
        j = j_v;
        k = k_v;
        q = 0;
    }
};

// class Esfera{
//     public:
//     Vector centro;
//     int raio;
//     Esfera(Vector& centro_esf, int& raio_esf)
//         : centro(centro_esf), raio(raio_esf) {}

//     Vector getCentro() { return centro; }
//     int getRaio() { return raio; }
// };


Vector subtrai_pontos(Point&  p1, Point&  p2){
    Vector sub(p1.x- p2.x, p1.y - p2.y, p1.z - p2.z);
    return sub;
}

Vector calcula_dr(Point& Po, Point& Pj){
    Vector Dr = subtrai_pontos(Po, Pj);
    int drNorma = calcula_norma(Dr);
    Vector dr(Dr.i/drNorma, Dr.j/drNorma, Dr.k/drNorma);
    return dr;
}

int calcula_norma(Vector& v){
    return sqrt(v.i*v.i + v.j*v.j + v.k*v.k);
}


