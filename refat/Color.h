#ifndef COLOR_H
#define COLOR_H

class Plano;
class Vector;
class Cilindro;
class Point;

class Color {
public:
    float r, g, b;
    
    Color(float r_c, float g_c, float b_c);
};

// Funções que retornam Color
Color calcula_Plano(Plano P, Vector dr, Color K_e, Color K_d);
Color calcula_color_cil(Cilindro cilindro, float t_cil, Vector dr);
Color calculaCone(float t, Point p);

#endif // COLOR_H
