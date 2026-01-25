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

Color CalcularCor(const Cena& cena, float t, const Vector& dir, 
    const Color& K_a, const Color& K_d, const Color& K_e, Vector normal, Objeto* obj) {
    Point P = calcula_eq_ray(cena.observador, t, dir);
    
    bool naSombraSpot = false;
    Color I_spot(0.0f, 0.0f, 0.0f);
    Color Id_spot(0.0f, 0.0f, 0.0f);
    Color Ie_spot(0.0f, 0.0f, 0.0f);


    //spot**********************************************
    if(cena.luzSpotAtiva) {
        Vector l_spot = calcula_l(cena.luzSpot.pos, P);
        Point P_mod_spot(P.x + l_spot.i * 1e-4f, P.y + l_spot.j * 1e-4f, P.z + l_spot.k * 1e-4f);
        float dist_Pi_Luz_spot = calcula_norma(subtrai_pontos(cena.luzSpot.pos, P_mod_spot));
        
        naSombraSpot = cena.estaNaSombra(P_mod_spot, l_spot, dist_Pi_Luz_spot, obj);
        
        Vector r1_spot = calcula_esc_por_vetor(2.0f, calcula_esc_por_vetor(calcula_prod_esc(normal, l_spot), normal));
        Vector r_spot(r1_spot.i - l_spot.i, r1_spot.j - l_spot.j, r1_spot.k - l_spot.k);
        
        float fd_spot = lidarExcecao(calcula_prod_esc(normal, l_spot));
        float cosAlpha_spot = lidarExcecao(calcula_prod_esc(r_spot, Vector(-dir.i, -dir.j, -dir.k)));
        float fe_spot = pow(cosAlpha_spot, cena.expoenteEspecular);

        float cos_dr_l = calcula_prod_esc(cena.luzSpot.direcao, 
            calcula_esc_por_vetor(-1.0f, l_spot));

        float cos_ang = cosf(cena.luzSpot.angulo * (3.14159265f / 180.0f));
        
        if(cos_dr_l >= cos_ang) {
            Id_spot = Color(cena.luzSpot.intensidade.r*cos_dr_l * K_d.r * fd_spot,
            cena.luzSpot.intensidade.g*cos_dr_l * K_d.g * fd_spot,
            cena.luzSpot.intensidade.b*cos_dr_l * K_d.b * fd_spot);

            Ie_spot = Color(cena.luzSpot.intensidade.r*cos_dr_l * K_e.r * fe_spot,
                cena.luzSpot.intensidade.g*cos_dr_l * K_e.g * fe_spot,
                cena.luzSpot.intensidade.b*cos_dr_l * K_e.b * fe_spot);

            I_spot = Color(Id_spot.r + Ie_spot.r, 
                Id_spot.g + Ie_spot.g, 
                Id_spot.b + Ie_spot.b);
        }        
    }
    
    //pontual**********************************************
    Vector l_pontual = calcula_l(cena.luz.pos, P);
    Point P_mod_pontual(P.x + l_pontual.i * 1e-4f, P.y + l_pontual.j * 1e-4f, P.z + l_pontual.k * 1e-4f);
    float dist_Pi_Luz_pontual = calcula_norma(subtrai_pontos(cena.luz.pos, P_mod_pontual));
    
    bool naSombra = cena.estaNaSombra(P_mod_pontual, l_pontual, dist_Pi_Luz_pontual, obj);

    Vector v(-dir.i, -dir.j, -dir.k);
    Vector r1_pontual = calcula_esc_por_vetor(2.0f, calcula_esc_por_vetor(calcula_prod_esc(normal, l_pontual), normal));
    Vector r_pontual(r1_pontual.i - l_pontual.i, r1_pontual.j - l_pontual.j, r1_pontual.k - l_pontual.k);

    float fd_pontual = lidarExcecao(calcula_prod_esc(normal, l_pontual));
    float cosAlpha_pontual = lidarExcecao(calcula_prod_esc(r_pontual, v));
    float fe_pontual = pow(cosAlpha_pontual, cena.expoenteEspecular);

    Color Ia(cena.luzAmbiente.r * K_a.r, cena.luzAmbiente.g * K_a.g, cena.luzAmbiente.b * K_a.b);
    if (naSombra || naSombraSpot) {
        Color I(lidarExcecao(Ia.r), lidarExcecao(Ia.g), lidarExcecao(Ia.b));
        int R = static_cast<int>(roundf(I.r * 255.0f));
        int G = static_cast<int>(roundf(I.g * 255.0f));
        int B = static_cast<int>(roundf(I.b * 255.0f));
        return Color(R, G, B);
    }

    Color Id_pontual(cena.luz.intensidade.r * K_d.r * fd_pontual,
        cena.luz.intensidade.g * K_d.g * fd_pontual,
        cena.luz.intensidade.b * K_d.b * fd_pontual);

    Color Ie_pontual(cena.luz.intensidade.r * K_e.r * fe_pontual,
        cena.luz.intensidade.g * K_e.g * fe_pontual,
        cena.luz.intensidade.b * K_e.b * fe_pontual);
    
    //calcula da cor final*******************************

    Color I_pontual(Id_pontual.r + Ie_pontual.r, 
        Id_pontual.g + Ie_pontual.g, 
        Id_pontual.b + Ie_pontual.b);
    
    
    Color I(lidarExcecao(I_pontual.r + I_spot.r + Ia.r), 
        lidarExcecao(I_pontual.g + I_spot.g + Ia.g), 
        lidarExcecao(I_pontual.b + I_spot.b + Ia.b));  

    int R = static_cast<int>(roundf(I.r * 255.0f));
    int G = static_cast<int>(roundf(I.g * 255.0f));
    int B = static_cast<int>(roundf(I.b * 255.0f));
    return Color(R, G, B);
}
