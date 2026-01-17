#ifndef MATRIZ3X3_H
#define MATRIZ3X3_H

class Vector;

class Matriz3x3 {
public:
    float m[3][3];

    Matriz3x3(float m00, float m01, float m02,
              float m10, float m11, float m12,
              float m20, float m21, float m22);
};

// Funções que retornam Matriz3x3
Matriz3x3 outerProduto(const Vector& d);
Matriz3x3 matrizSubtrai(const Matriz3x3& A, const Matriz3x3& B);
Matriz3x3 escalarMatriz(float s, const Matriz3x3& M);
Matriz3x3 transpostaVetor(const Vector& v);

// Função auxiliar envolvendo Matriz3x3
float transpostoPorVetor(const Matriz3x3& M, const Vector& v);

#endif // MATRIZ3X3_H
