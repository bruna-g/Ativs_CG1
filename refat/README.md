como rodar:

g++ -O2 -std=c++17 -Wall -Wextra   ativ5.cpp Vector.cpp Point.cpp Color.cpp Plano.cpp Cilindro.cpp Cone.cpp Cubo.cpp Matriz3x3.cpp Utils.cpp Textura.cpp   $(sdl2-config --cflags --libs) -o ativ

./ativ