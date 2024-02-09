# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)

# Set the system name
set(CMAKE_SYSTEM_NAME Linux)

# No need to set a compiler vendor for generic setup

# Use GCC for C compiler
set(CMAKE_C_COMPILER gcc)
# Update C flags: Removed MKL references, ensure C++17 standard, and retain optimization and warning flags
set(CMAKE_C_FLAGS "-O3 -Wall -Wextra -fopenmp" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_RELEASE " ")
set(CMAKE_C_FLAGS_DEBUG "-g -fopenmp")

# Use G++ for C++ compiler
set(CMAKE_CXX_COMPILER g++)
# Update C++ flags: Removed MKL references, explicitly set to C++17 standard, and retain optimization and warning flags
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -fopenmp" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "")
set(CMAKE_CXX_FLAGS_DEBUG "-g -fopenmp")

# Removed INTEL_MKL setting as we're not using MKL

# For BLAS, you may need to link against the system-provided BLAS library.
# This is often provided by packages like libblas-dev or similar on many Linux distributions.
# Adjust the BLAS_LIB variable as necessary for your specific BLAS library or use FindBLAS CMake module.
find_package(BLAS REQUIRED)
if(BLAS_FOUND)
  set(BLAS_LIB ${BLAS_LIBRARIES})
endif()

# For Eigen, since it's a header-only library, you typically only need to include its path.
# Adjust the path as necessary based on where Eigen is located on your system or use FindEigen3 CMake module.
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Example of how to include directories for your target
# include_directories(${EIGEN3_INCLUDE_DIR})

# When linking your target, use the variables for libraries as needed.
# For example:
# target_link_libraries(my_target ${BLAS_LIB})
