cmake_minimum_required(VERSION 3.10)
project(autoc)

find_package(VTK REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED COMPONENTS thread system serialization log)
find_package(AWSSDK REQUIRED COMPONENTS sqs s3)

option(USE_SANITIZER "Build with Address Sanitizer" OFF)
if(USE_SANITIZER)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
endif()

include(${VTK_USE_FILE} )
include_directories(../include ${EIGEN3_INCLUDE_DIR} ${VTK_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${AWSSDK_INCLUDE_DIRS})
link_directories(../lib)
add_compile_options(-Wno-write-strings)

add_executable(autoc autoc.cc autoc.h minisim.h pathgen.cc pathgen.h threadpool.h logger.cc logger.h)
add_executable(minisim minisim.cc minisim.h logger.cc logger.h)
add_executable(renderer renderer.cc renderer.h minisim.h autoc.h)

target_link_libraries(autoc gp ${Boost_LIBRARIES} Eigen3::Eigen ${AWSSDK_LIBRARIES})
target_link_libraries(minisim ${Boost_LIBRARIES} Eigen3::Eigen ${AWSSDK_LIBRARIES})
target_link_libraries(renderer ${VTK_LIBRARIES} ${Boost_LIBRARIES} Eigen3::Eigen ${AWSSDK_LIBRARIES})
