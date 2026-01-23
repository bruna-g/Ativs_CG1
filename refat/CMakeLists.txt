cmake_minimum_required(VERSION 3.16)

project(refat LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Coloca binários diretamente em refat/
# (útil quando você mantém a pasta build dentro do refat)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}")
foreach(OUTPUTCONFIG Debug Release RelWithDebInfo MinSizeRel)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} "${PROJECT_SOURCE_DIR}")
endforeach()

# Flags úteis
if(MSVC)
  add_compile_options(/W4)
else()
  add_compile_options(-Wall -Wextra)
endif()

# SDL2 via pkg-config (Linux)
find_package(PkgConfig REQUIRED)
pkg_check_modules(SDL2 REQUIRED sdl2)

# Paths do projeto
set(REFAT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
set(REFAT_SRC_DIR "${REFAT_DIR}/src")
set(REFAT_INC_DIR "${REFAT_DIR}/include")

# Fontes comuns (compartilhadas pelos dois executáveis)
set(REFAT_COMMON_SOURCES
  "${REFAT_SRC_DIR}/Cena.cpp"
  "${REFAT_SRC_DIR}/Vector.cpp"
  "${REFAT_SRC_DIR}/Point.cpp"
  "${REFAT_SRC_DIR}/Color.cpp"
  "${REFAT_SRC_DIR}/Plano.cpp"
  "${REFAT_SRC_DIR}/Cilindro.cpp"
  "${REFAT_SRC_DIR}/Cone.cpp"
  "${REFAT_SRC_DIR}/Cubo.cpp"
  "${REFAT_SRC_DIR}/Esfera.cpp"
  "${REFAT_SRC_DIR}/Matriz3x3.cpp"
  "${REFAT_SRC_DIR}/Matriz4x4.cpp"
  "${REFAT_SRC_DIR}/Utils.cpp"
  "${REFAT_SRC_DIR}/Textura.cpp"
  "${REFAT_SRC_DIR}/Camera.cpp"
  "${REFAT_SRC_DIR}/Malha.cpp"
)

# Executável SDL (janela)
add_executable(main
  "${REFAT_SRC_DIR}/main.cpp"
  "${REFAT_SRC_DIR}/setup.cpp"
  ${REFAT_COMMON_SOURCES}
)

# Executável PPM (ray casting em arquivo)
# add_executable(ativ
#   "${REFAT_SRC_DIR}/ativ5.cpp"
#   ${REFAT_COMMON_SOURCES}
# )

foreach(target IN ITEMS main ) #ativ
  target_include_directories(${target} PRIVATE
    "${REFAT_INC_DIR}"
    ${SDL2_INCLUDE_DIRS}
  )

  target_compile_options(${target} PRIVATE
    ${SDL2_CFLAGS_OTHER}
  )

  target_link_directories(${target} PRIVATE
    ${SDL2_LIBRARY_DIRS}
  )

  target_link_libraries(${target} PRIVATE
    ${SDL2_LIBRARIES}
  )
endforeach()
