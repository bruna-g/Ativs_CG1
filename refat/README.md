Como rodar

Com CMake (recomendado)

Os comandos abaixo já fazem: configurar + compilar + executar.

- Build + executar `ativ` (gera PPM):

```bash
cd refat
```

```bash
cmake -S . -B build && cmake --build build -j && ./ativ
```

- Build + executar `main` (SDL):

```bash
cmake -S . -B build && cmake --build build -j && ./main
```

Observações:
- Os executáveis são gerados diretamente em `refat/`.
- O `main` abre uma janela SDL e fica rodando até você fechar.

Build manual (g++), opcional

Dentro de `refat/src`:

- Para gerar PPM:

```bash
cd refat/src && g++ -O2 -std=c++17 -Wall -Wextra -I../include \
  ativ5.cpp Vector.cpp Point.cpp Color.cpp Plano.cpp Cilindro.cpp Cone.cpp Cubo.cpp Esfera.cpp \
  Matriz3x3.cpp Utils.cpp Textura.cpp ../include/Malha.cpp \
  $(sdl2-config --cflags --libs) -o ../ativ
```

```bash
cd refat && ./ativ
```

- Para gerar SDL:

```bash
cd refat/src && g++ -O2 -std=c++17 -Wall -Wextra -I../include \
  main.cpp Vector.cpp Point.cpp Color.cpp Plano.cpp Cilindro.cpp Cone.cpp Cubo.cpp Esfera.cpp \
  Matriz3x3.cpp Utils.cpp Textura.cpp setup.cpp ../include/Malha.cpp \
  $(sdl2-config --cflags --libs) -o ../main
```

```bash
cd refat && ./main
```