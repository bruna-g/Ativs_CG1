#ifndef VECTOR_H
#define VECTOR_H

class Point;
class Matriz3x3;
class Cilindro;

class Vector {
public:
    float i, j, k, q;

    Vector(float i_v, float j_v, float k_v);
    Vector(float i_v, float j_v, float k_v, float q_v);
};

// Funções que retornam Vector
Vector cross(const Vector& a, const Vector& b);
Vector matrizVetorProduto(const Matriz3x3& M, const Vector& v);
Vector matrizPorMatriz(const Matriz3x3& A, const Matriz3x3& B);

#endif // VECTOR_H
