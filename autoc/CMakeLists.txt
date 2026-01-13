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
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -fno-omit-frame-pointer -fno-optimize-sibling-calls")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -fno-omit-frame-pointer -fno-optimize-sibling-calls")
endif()

if(PERFORMANCE_BUILD)
    message(STATUS "Enabling performance-optimized build for autoc (with debug symbols)")
    # -O3: Aggressive optimization
    # -g: Debug symbols for core dump analysis (doesn't affect runtime)
    # -march=native: Use AVX-512/AVX2 and other CPU-specific instructions
    # -mtune=native: Optimize instruction scheduling for this CPU
    # -funroll-loops: Better SIMD utilization
    # -flto: Link-time optimization for whole-program analysis
    # -ffast-math: Aggressive FP optimizations (safe now that simTime is deterministic)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -march=native -mtune=native -funroll-loops -flto -ffast-math")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -g -march=native -mtune=native -funroll-loops -flto -ffast-math")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -flto")
endif()

include(${VTK_USE_FILE})
include_directories(../include ${EIGEN3_INCLUDE_DIR} ${VTK_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${AWSSDK_INCLUDE_DIRS})
link_directories(../lib)
add_compile_options(-Wno-write-strings)
add_definitions(-DGP_BUILD)

add_executable(autoc autoc.cc autoc.h autoc-eval.cc minisim.h aircraft_state.h pathgen.cc pathgen.h threadpool.h logger.cc logger.h gp_bytecode.cc gp_bytecode.h gp_evaluator_portable.cc gp_evaluator_portable.h gp_evaluator_desktop.cc gp_evaluator_desktop.h config_manager.cc config_manager.h)
add_executable(minisim minisim.cc minisim.h aircraft_state.h autoc.h autoc-eval.cc gp_bytecode.cc gp_bytecode.h gp_evaluator_portable.cc gp_evaluator_portable.h gp_evaluator_desktop.cc gp_evaluator_desktop.h)
add_executable(renderer renderer.cc renderer.h minisim.h aircraft_state.h autoc.h config_manager.cc config_manager.h pathgen.cc pathgen.h gp_types.h)
add_executable(gpextractor gpextractor.cc minisim.h aircraft_state.h config_manager.cc config_manager.h autoc-eval.cc gp_evaluator_portable.cc gp_evaluator_portable.h gp_evaluator_desktop.cc gp_evaluator_desktop.h gp_bytecode.cc gp_bytecode.h autoc.h threadpool.h logger.cc logger.h pathgen.cc pathgen.h gp_types.h)
add_executable(bytecode2cpp bytecode2cpp.cc gp_bytecode.cc gp_bytecode.h gp_evaluator_portable.cc gp_evaluator_portable.h gp_evaluator_desktop.cc gp_evaluator_desktop.h config_manager.cc config_manager.h)

target_link_libraries(autoc gp ${Boost_LIBRARIES} Eigen3::Eigen ${AWSSDK_LIBRARIES})
target_link_libraries(minisim gp ${Boost_LIBRARIES} Eigen3::Eigen)
target_link_libraries(renderer gp ${VTK_LIBRARIES} ${Boost_LIBRARIES} Eigen3::Eigen ${AWSSDK_LIBRARIES})
target_link_libraries(gpextractor gp ${Boost_LIBRARIES} Eigen3::Eigen ${AWSSDK_LIBRARIES})
target_link_libraries(bytecode2cpp gp ${Boost_LIBRARIES} Eigen3::Eigen ${AWSSDK_LIBRARIES})

if(ENABLE_AUTOC_TESTS)
  add_executable(autoc_tests tests/gp_evaluator_tests.cc gp_evaluator_portable.cc)
  target_compile_definitions(autoc_tests PRIVATE GP_BUILD GP_TEST)
  target_link_libraries(autoc_tests GTest::gtest GTest::gtest_main Eigen3::Eigen)
  enable_testing()
  add_test(NAME autoc_tests COMMAND autoc_tests)
  add_custom_target(run_autoc_tests ALL DEPENDS autoc_tests COMMAND autoc_tests)
endif()
