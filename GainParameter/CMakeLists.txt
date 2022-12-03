#---Setup the example project---------------------------------------------------
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(CeramicGEM)
#-------------------------------------------------------------------------------
# ROOT
find_package(ROOT)
if(NOT ROOT_FOUND)
    message(STATUS "Could not find ROOT")
endif()
#--------------------------------------------------------------------------------
# Garfield
if(NOT TARGET Garfield::Garfield)
    find_package(Garfield REQUIRED)
endif()
#--------------------------------------------------------------------------------
# Add the executable, and link it to the dependent libraries.
add_executable(${PROJECT_NAME} ${PROJECT_NAME}.C )
target_link_libraries(${PROJECT_NAME} Garfield::Garfield ${ROOT_LIBRARIES})
#---------------------------------------------------------------------------------
