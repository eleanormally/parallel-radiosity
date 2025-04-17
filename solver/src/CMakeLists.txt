cmake_minimum_required (VERSION 3.10)
project(radiosity)
set(executable radiosity)
set (CMAKE_EXPORT_COMPILE_COMMANDS ON)

#including source .cpp files
file(GLOB_RECURSE PROJECT_FILES *.cpp)
add_executable(${executable}
  ${PROJECT_FILES}
)

add_custom_command(
  OUTPUT ${PROJECT_BINARY_DIR}/calculate.o
  COMMAND nvcc ${PROJECT_SOURCE_DIR}/calculate.cu -c -o ${PROJECT_BINARY_DIR}/calculate.o
  DEPENDS ${PROJECT_SOURCE_DIR}/calculate.cu
)
add_custom_target(buildCuda DEPENDS ${PROJECT_BINARY_DIR}/calculate.o)
add_library(Calculate STATIC IMPORTED GLOBAL)
add_dependencies(Calculate buildCuda)
set_target_properties(Calculate PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/calculate.o)

#preparing gcc
set_target_properties(${executable} PROPERTIES CXX_STANDARD 11)
find_package(MPI REQUIRED)
find_package(CUDA 11.2 REQUIRED)
target_link_libraries(Calculate INTERFACE ${CUDA_LIBRARIES})
target_link_libraries(${executable} MPI::MPI_CXX)
target_link_libraries(${executable} Calculate)
