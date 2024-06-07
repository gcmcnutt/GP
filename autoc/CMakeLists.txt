cmake_minimum_required(VERSION 3.10)
project(autoc)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE} )
include_directories(../include ${VTK_INCLUDE_DIRS})
link_directories(../lib)
add_compile_options(-Wno-write-strings)

add_executable(autoc autoc.cc minisim.h minisim.cc pathgen.cc pathgen.h )

target_link_libraries(autoc gp ${VTK_LIBRARIES})