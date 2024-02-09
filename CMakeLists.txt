# Based on https://github.com/Barthelemy/CppProjectTemplate
# Assuming the usage of a more modern CMake version for better features and compatibility
cmake_minimum_required(VERSION 3.15) # Bumping this up for better modern CMake features
project(LinQT VERSION 2.0.0 LANGUAGES CXX)

# Note: To use a specific toolchain file, pass it via the command line to cmake with 
# '-DCMAKE_TOOLCHAIN_FILE=path/to/your/toolchain.cmake'
# Example: cmake -DCMAKE_TOOLCHAIN_FILE=../toolchain-gcc.cmake ..

# Set the C++ standard explicitly
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Installation variables
set(PYTHON_VERSION "3.6")
set(ENV_PYTHON_ROOT $ENV{CONDA_PREFIX})

# CMake useful variables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib") 
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# Output directories
message("Project Source Directory: ${PROJECT_SOURCE_DIR}")

# Including the project's header files
# Note: Prefer to use target_include_directories on specific targets rather than include_directories
# include_directories(${PROJECT_SOURCE_DIR}/include) 

# Adding subdirectories
add_subdirectory(src)
#add_subdirectory(utilities)
# If you have these directories and relevant CMakeLists.txt in them, you can uncomment these lines
# add_subdirectory(doc)
# add_subdirectory(test)
# add_subdirectory(samples)

# Example of defining a library or executable target and setting its include directories
# add_library(mylib SHARED src/mylib.cpp)
# target_include_directories(mylib PUBLIC ${PROJECT_SOURCE_DIR}/include)

# Example of defining an executable target
# add_executable(myapp src/main.cpp)
# target_link_libraries(myapp PRIVATE mylib)
