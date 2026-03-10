cmake_minimum_required(VERSION 3.10)
project(autoc)

find_package(VTK REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED COMPONENTS thread system serialization log)
find_package(AWSSDK REQUIRED COMPONENTS core s3)

option(ENABLE_AUTOC_TESTS "Build autoc unit tests" ON)
if(ENABLE_AUTOC_TESTS)
  include(GoogleTest)
  include(FetchContent)
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_Declare(
      googletest
      URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
  )
  FetchContent_MakeAvailable(googletest)
endif()

option(USE_SANITIZER "Build with Address Sanitizer" OFF)
option(ENABLE_PROFILING "Build with frame pointers and debug info for perf/FlameGraph" OFF)
option(PERFORMANCE_BUILD "Build with high performance flags (keeps debug symbols for core dumps)" OFF)

set(CMAKE_EXE_LINKER_FLAGS "-no-pie")

if(USE_SANITIZER)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
endif()

if(ENABLE_PROFILING)
    message(STATUS "Enabling profiling-friendly compiler flags for autoc")
    # -O2: Optimize but keep code structure recognizable (avoid -O3 which over-inlines)
    # -g: Debug symbols for perf to resolve symbols
    # -fno-omit-frame-pointer: Essential for perf call graph
    # -fno-optimize-sibling-calls: Preserve tail call frames in stack traces
    # Note: -O2 gives representative percentages vs PERFORMANCE_BUILD's -O3
    #       (difference is typically <5% in time distribution)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g -fno-omit-frame-pointer -fno-optimize-sibling-calls")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -g -fno-omit-frame-pointer -fno-optimize-sibling-calls")
endif()

if(PERFORMANCE_BUILD)
    message(STATUS "Enabling performance-optimized build for autoc (with debug symbols)")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -march=native -mtune=native -funroll-loops -flto=auto")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -g -march=native -mtune=native -funroll-loops -flto=auto")
endif()

include(${VTK_USE_FILE})
include_directories(../include ${EIGEN3_INCLUDE_DIR} ${VTK_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${AWSSDK_INCLUDE_DIRS})
link_directories(../lib)
add_compile_options(-Wno-write-strings)
add_definitions(-DGP_BUILD)

# Shared library for common sources (compiled once, linked by all executables)
add_library(autoc_common STATIC
    gp_bytecode.cc
    gp_bytecode.h
    gp_evaluator_portable.cc
    gp_evaluator_portable.h
    gp_evaluator_desktop.cc
    gp_evaluator_desktop.h
    config_manager.cc
    config_manager.h
    logger.cc
    logger.h
    pathgen.cc
    pathgen.h
    autoc-eval.cc
    autoc.h
    minisim.h
    aircraft_state.h
    threadpool.h
    gp_types.h
    variation_generator.h
)
target_link_libraries(autoc_common gp ${Boost_LIBRARIES} Eigen3::Eigen)

add_executable(autoc autoc.cc)
add_executable(minisim minisim.cc)
add_executable(renderer renderer.cc)
add_executable(gpextractor gpextractor.cc)
add_executable(bytecode2cpp bytecode2cpp.cc)

target_link_libraries(autoc autoc_common ${AWSSDK_LIBRARIES})
target_link_libraries(minisim autoc_common)
target_link_libraries(renderer autoc_common ${VTK_LIBRARIES} ${AWSSDK_LIBRARIES})
target_link_libraries(gpextractor autoc_common ${AWSSDK_LIBRARIES})
target_link_libraries(bytecode2cpp autoc_common ${AWSSDK_LIBRARIES})

if(ENABLE_AUTOC_TESTS)
  # Tests need separate compilation with GP_TEST define, can't use autoc_common
  add_executable(autoc_tests tests/gp_evaluator_tests.cc gp_evaluator_portable.cc)
  target_compile_definitions(autoc_tests PRIVATE GP_BUILD GP_TEST)
  target_link_libraries(autoc_tests GTest::gtest GTest::gtest_main Eigen3::Eigen)
  enable_testing()
  add_test(NAME autoc_tests COMMAND autoc_tests)
  # --gtest_brief=1 only prints failures, not every test run
  add_custom_target(run_autoc_tests ALL DEPENDS autoc_tests COMMAND autoc_tests --gtest_brief=1)
endif()
