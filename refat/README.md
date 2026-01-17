como rodar:

cd src

Para gerar PPM:
g++ -O2 -std=c++17 -Wall -Wextra -I../include \
	ativ5.cpp Vector.cpp Point.cpp Color.cpp Plano.cpp Cilindro.cpp Cone.cpp Cubo.cpp Esfera.cpp \
	Matriz3x3.cpp Utils.cpp Textura.cpp ../include/Malha.cpp \
	$(sdl2-config --cflags --libs) -o ativ

./ativ

Para gerar SDL:
g++ -O2 -std=c++17 -Wall -Wextra -I../include \
	main.cpp Vector.cpp Point.cpp Color.cpp Plano.cpp Cilindro.cpp Cone.cpp Cubo.cpp Esfera.cpp \
	Matriz3x3.cpp Utils.cpp Textura.cpp setup.cpp ../include/Malha.cpp \
	$(sdl2-config --cflags --libs) -o main

./main