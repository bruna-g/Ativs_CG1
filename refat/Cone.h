#ifndef CONE_H
#define CONE_H

#include "Point.h"

class Cone {
public:
    Point cb;
    Point v;
    float raio;
    
    Cone(const Point& cb_c, const Point& v_c, const float raio_c);
};

// Função relacionada ao Cone
float calcula_t_cone(Point& Pj);

#endif // CONE_H
