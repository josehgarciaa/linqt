cmake_minimum_required(VERSION 3.15) # A version that supports most modern CMake features
project(KPM_LIB VERSION 1.0 LANGUAGES CXX) # Assuming your project is C++ based

# Author Information
set(AUTHOR "Jose H. Garcia")

#CREATES KPM_LIB
set(KPM_LIB_SOURCES
dummy.cpp
#  chebyshev_solver.cpp
#  quantum_states.cpp
#  chebyshev_vectors.cpp
#  chebyshev_moments.cpp
#  chebyshev_moments2D.cpp
#  chebyshev_moments1D.cpp
#  chebyshev_momentsTD.cpp
#  sparse_matrix.cpp
#  special_functions.cpp
)



#CREATION OF THE LIBRARY
add_library(kpm_lib STATIC ${KPM_LIB_SOURCES})
target_include_directories(kpm_lib PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(kpm_lib PUBLIC ${BLAS_LIB})

# Function to simplify executable creation
function(create_kpm_executable name)
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} kpm_lib)
endfunction()

# Create the different executables
create_kpm_executable(kpm_spectral_green_function)

#create_kpm_executable(inline_compute-kpm-nonEqOp)
#create_kpm_executable(inline_compute-kpm-CorrOp)
#create_kpm_executable(inline_compute-kpm-spectralOp)
#create_kpm_executable(inline_spectralFunctionFromChebmom)
#create_kpm_executable(inline_kuboGreenwoodFromChebmom)
#create_kpm_executable(inline_kuboBastinIFromChebmom)
#create_kpm_executable(inline_kuboBastinIIFromChebmom)
#create_kpm_executable(inline_kuboBastinFromChebmom)
#create_kpm_executable(inline_kuboBastinSeaFromChebmom)
#create_kpm_executable(inline_kuboBastinSurfFromChebmom)
#create_kpm_executable(inline_timeCorrelationsFromChebmom)
