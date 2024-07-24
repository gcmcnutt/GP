cmake_minimum_required(VERSION 3.10)
project(autoc)

find_package(VTK REQUIRED)
find_package(Eigen3 REQUIRED)
include(${VTK_USE_FILE} )
include_directories(../include ${EIGEN3_INCLUDE_DIR} ${VTK_INCLUDE_DIRS})
link_directories(../lib)
add_compile_options(-Wno-write-strings)

add_executable(autoc autoc.cc autoc.h minisim.h minisim.cc pathgen.cc pathgen.h renderer.cc renderer.h)

target_link_libraries(autoc gp ${VTK_LIBRARIES})