project(lab1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(lab1Galerkin mainGalerkin.cpp functions.cpp functions.h GalerkinMethod.cpp GalerkinMethod.h)
add_executable(lab1Ritz mainRitz.cpp functions.cpp functions.h RitzMethod.cpp RitzMethod.h)
add_executable(lab1Test test.cpp functions.cpp functions.h GalerkinMethod.cpp GalerkinMethod.h RitzMethod.cpp RitzMethod.h)