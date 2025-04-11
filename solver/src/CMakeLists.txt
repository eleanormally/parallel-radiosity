
cmake_minimum_required (VERSION 3.10)
project(radiosity)
set(executable radiosity)
set (CMAKE_EXPORT_COMPILE_COMMANDS ON)

#including source .cpp files
set(PROJECT_FILES
  ${PROJECT_SOURCE_DIR}/main.cpp
  ${PROJECT_SOURCE_DIR}/mesh/*.cpp
)
add_executable(${executable}
  ${PROJECT_FILES}
)

#preparing gcc
set_target_properties(${executable} PROPERTIES CXX_STANDARD 11)
find_package(MPI REQUIRED)
target_link_libraries(${executable} MPI::MPI_CXX)
set_target_properties(${executable} PROPERTIES COMPILE_FLAGS "-g -Wall -Wextra")
