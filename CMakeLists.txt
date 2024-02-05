project(lab1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(lab1 main.cpp functions.cpp functions.h GalerkinMethod.cpp GalerkinMethod.h)