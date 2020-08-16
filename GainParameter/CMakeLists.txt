#---Setup the example project---------------------------------------------------
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(CeramicGEM)
#---------------------------------------------------------------------------------------------------
# ROOT
find_package(ROOT)
if(NOT ROOT_FOUND)
        message(STATUS "Could not find ROOT")
else()
        message(STATUS "Found package ROOT in ${ROOT_SOURCE_DIR}")
        include_directories(${ROOT_INCLUDE_DIRS})
        link_directories($ENV{ROOTSYS}/lib)
endif()
#--------------------------------------------------------------------------------------------------
# Garfield
find_package(Garfield REQUIRED)
if(NOT Garfield_FOUND)
        message(STATUS "Could not find cmake configuration for package Garfield.")
else()
        message(STATUS "Found package Garfield++ in ${GARFIELD_SOURCE_DIR}")
        include_directories(${GARFIELD_INCLUDE_DIRS})
endif()
#--------------------------------------------------------------------------------------------------
# Add the executable, and link it to the dependent libraries.
add_executable(${PROJECT_NAME} ${PROJECT_NAME}.C )
target_link_libraries(${PROJECT_NAME} Garfield -lgfortran -lGeom -lGdml ${ROOT_LIBRARIES})
#--------------------------------------------------------------------------------------------------
# Copy ansys script to the build directory.
file(GLOB ANSYS_SCRIPTS RELATIVE ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/*.lis)
foreach(_script ${ANSYS_SCRIPTS})
        configure_file(
        ${PROJECT_SOURCE_DIR}/${_script}
        ${PROJECT_BINARY_DIR}/${_script}
        COPYONLY
        )
endforeach()