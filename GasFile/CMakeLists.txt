# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(Generate)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()

# ---Build executables----------------------------------------------------------
add_executable(generate generate.C)
target_link_libraries(generate Garfield::Garfield)

add_executable(read read.C)
target_link_libraries(read Garfield::Garfield)
