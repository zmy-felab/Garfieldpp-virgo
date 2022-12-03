cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(x-ray)
find_package(ROOT)

if(NOT TARGET Garfield::Garfield)
    find_package(Garfield REQUIRED)
endif()

add_executable(${PROJECT_NAME} ${PROJECT_NAME}.C)
target_link_libraries(${PROJECT_NAME} Garfield::Garfield ${ROOT_LIBRARIES})

#find_package(OpenMP)
#if(OpenMP_CXX_FOUND)
#    target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)
#endif()
