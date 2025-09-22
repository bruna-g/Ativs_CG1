#include <iostream>
#include <vector>
using namespace std;

int main(){
    float wJanela = 5.0;
    float hJanela = 5.0;
    float dJanela = 2.0;

    float rEsfera = 3.0;
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

        }
    }
    return 0;
}

class Esfera{
    private:
        int centro[3];
        int raio;
//        int p[3];
    public:
        Esfera(int centro_esf[3], int raio_esf){
            centro[0] = centro_esf[0];
            centro[1] = centro_esf[1];
            centro[2] = centro_esf[2];
            raio = raio_esf;
        }

        int* getCentro() {return centro;}
        int getRaio() {return raio;}

        bool intersecta(Raio (int p0[3], int pj[3]), Esfera(int centro[3], int raio)){

            return true;
        }
};

class Raio{
    private:
        int Po[3]; 
        int Pj[3];
    public:
        Raio(int ponto_inicial[3], int ponto_final[3]){
            Po[0] = ponto_inicial[0];
            Po[1] = ponto_inicial[1];
            Po[2] = ponto_inicial[2];
            Pj[0] = ponto_final[0];
            Pj[1] = ponto_final[1];
            Pj[2] = ponto_final[2];
        }
        int* getPo() {return Po;}
        int* getPj() {return Pj;}
        int* getDr(int Po[3], int Pj[3]){
            int Dr[3];
            Dr[0] = Pj[0] - Po[0];
            Dr[1] = Pj[1] - Po[1];
            Dr[2] = Pj[2] - Po[2];
            return Dr;
        }

};

vector<int> substrai(int p1[3], int p2[3]){
    vector<int> sub[3];
    sub[0] = p1[0] - p2[0];
    sub[1] = p1[1] - p2[1];
    sub[2] = p1[2] - p2[2];
    return sub;
}