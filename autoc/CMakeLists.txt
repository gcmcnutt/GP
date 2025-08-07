cmake_minimum_required(VERSION 3.10)
project(autoc)

find_package(VTK REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED COMPONENTS thread system serialization log)
find_package(AWSSDK REQUIRED COMPONENTS core s3)

option(USE_SANITIZER "Build with Address Sanitizer" OFF)

set(CMAKE_EXE_LINKER_FLAGS "-no-pie")

if(USE_SANITIZER)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
endif()

include(${VTK_USE_FILE})
include_directories(../include ${EIGEN3_INCLUDE_DIR} ${VTK_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${AWSSDK_INCLUDE_DIRS})
link_directories(../lib)
add_compile_options(-Wno-write-strings)
add_definitions(-DGP_BUILD)

add_executable(autoc autoc.cc autoc.h autoc-eval.cc minisim.h aircraft_state.h pathgen.cc pathgen.h threadpool.h logger.cc logger.h gp_bytecode.cc gp_bytecode.h gp_evaluator_portable.cc gp_evaluator_portable.h gp_evaluator_desktop.cc gp_evaluator_desktop.h config_manager.cc config_manager.h)
add_executable(minisim minisim.cc minisim.h aircraft_state.h autoc.h autoc-eval.cc gp_bytecode.cc gp_bytecode.h gp_evaluator_portable.cc gp_evaluator_portable.h gp_evaluator_desktop.cc gp_evaluator_desktop.h)
add_executable(renderer renderer.cc renderer.h minisim.h aircraft_state.h autoc.h config_manager.cc config_manager.h)
add_executable(gpextractor gpextractor.cc minisim.h aircraft_state.h config_manager.cc config_manager.h autoc-eval.cc gp_evaluator_portable.cc gp_evaluator_portable.h gp_evaluator_desktop.cc gp_evaluator_desktop.h)
add_executable(bytecode2cpp bytecode2cpp.cc gp_bytecode.cc gp_bytecode.h gp_evaluator_portable.cc gp_evaluator_portable.h gp_evaluator_desktop.cc gp_evaluator_desktop.h config_manager.cc config_manager.h)

target_link_libraries(autoc gp ${Boost_LIBRARIES} Eigen3::Eigen ${AWSSDK_LIBRARIES})
target_link_libraries(minisim gp ${Boost_LIBRARIES} Eigen3::Eigen)
target_link_libraries(renderer gp ${VTK_LIBRARIES} ${Boost_LIBRARIES} Eigen3::Eigen ${AWSSDK_LIBRARIES})
target_link_libraries(gpextractor gp ${Boost_LIBRARIES} Eigen3::Eigen ${AWSSDK_LIBRARIES})
target_link_libraries(bytecode2cpp gp ${Boost_LIBRARIES} Eigen3::Eigen ${AWSSDK_LIBRARIES})
