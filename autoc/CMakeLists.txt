cmake_minimum_required(VERSION 3.10)
project(autoc)

find_package(VTK REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(AWSSDK REQUIRED COMPONENTS core s3)

include(FetchContent)

# inih config parser
FetchContent_Declare(
    inih
    GIT_REPOSITORY https://github.com/benhoyt/inih.git
    GIT_TAG r58
)
FetchContent_Populate(inih)
add_library(inih STATIC
    ${inih_SOURCE_DIR}/ini.c
    ${inih_SOURCE_DIR}/cpp/INIReader.cpp
)
target_include_directories(inih PUBLIC ${inih_SOURCE_DIR} ${inih_SOURCE_DIR}/cpp)

# cereal serialization (header-only, replaces Boost.Serialization)
FetchContent_Declare(
    cereal
    GIT_REPOSITORY https://github.com/USCiLab/cereal.git
    GIT_TAG v1.3.2
    GIT_SHALLOW TRUE
)
set(SKIP_PERFORMANCE_COMPARISON ON CACHE BOOL "" FORCE)
set(BUILD_SANDBOX OFF CACHE BOOL "" FORCE)
set(BUILD_DOC OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(cereal)

option(ENABLE_AUTOC_TESTS "Build autoc unit tests" ON)
if(ENABLE_AUTOC_TESTS)
  include(GoogleTest)
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
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g -fno-omit-frame-pointer -fno-optimize-sibling-calls")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -g -fno-omit-frame-pointer -fno-optimize-sibling-calls")
endif()

if(PERFORMANCE_BUILD)
    message(STATUS "Enabling performance-optimized build for autoc (with debug symbols)")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -march=native -mtune=native -funroll-loops -flto=auto")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -g -march=native -mtune=native -funroll-loops -flto=auto")
endif()

include(${VTK_USE_FILE})
include_directories(${EIGEN3_INCLUDE_DIR} ${VTK_INCLUDE_DIRS} ${AWSSDK_INCLUDE_DIRS})
add_compile_options(-Wno-write-strings)

# Shared library for common sources (compiled once, linked by all executables)
add_library(autoc_common STATIC
    config_manager.cc
    config_manager.h
    logger.cc
    logger.h
    pathgen.cc
    pathgen.h
    sensor_math.cc
    sensor_math.h
    autoc.h
    minisim.h
    aircraft_state.h
    threadpool.h
    gp_types.h
    variation_generator.h
    eval_backend.h
    fitness_computer.cc
    fitness_computer.h
    eval_logger.cc
    eval_logger.h
    nn_evaluator_portable.cc
    nn_evaluator_portable.h
    nn_serialization.cc
    nn_serialization.h
    nn_population.cc
    nn_population.h
)
target_link_libraries(autoc_common inih cereal::cereal Eigen3::Eigen)

add_executable(autoc autoc.cc)
add_executable(minisim minisim.cc)
add_executable(renderer renderer.cc)
add_executable(nnextractor nnextractor.cc)
add_executable(nn2cpp nn2cpp.cc)

target_link_libraries(autoc autoc_common ${AWSSDK_LIBRARIES})
target_link_libraries(minisim autoc_common)
target_link_libraries(renderer autoc_common ${VTK_LIBRARIES} ${AWSSDK_LIBRARIES})
target_link_libraries(nnextractor autoc_common ${AWSSDK_LIBRARIES})
target_link_libraries(nn2cpp autoc_common)

if(ENABLE_AUTOC_TESTS)
  # NN evaluator tests
  add_executable(nn_evaluator_tests tests/nn_evaluator_tests.cc nn_evaluator_portable.cc)
  target_link_libraries(nn_evaluator_tests GTest::gtest GTest::gtest_main Eigen3::Eigen cereal::cereal)

  # NN serialization tests
  add_executable(nn_serialization_tests tests/nn_serialization_tests.cc nn_serialization.cc nn_evaluator_portable.cc)
  target_link_libraries(nn_serialization_tests GTest::gtest GTest::gtest_main Eigen3::Eigen cereal::cereal)

  # NN population tests
  add_executable(nn_population_tests tests/nn_population_tests.cc nn_population.cc nn_evaluator_portable.cc)
  target_link_libraries(nn_population_tests GTest::gtest GTest::gtest_main Eigen3::Eigen cereal::cereal)

  # Fitness computer tests
  add_executable(fitness_computer_tests tests/fitness_computer_tests.cc fitness_computer.cc)
  target_link_libraries(fitness_computer_tests GTest::gtest GTest::gtest_main Eigen3::Eigen)

  # Contract tests (014-nn-training-signal)
  add_executable(contract_evaluator_tests tests/contract_evaluator_tests.cc nn_evaluator_portable.cc)
  target_link_libraries(contract_evaluator_tests GTest::gtest GTest::gtest_main Eigen3::Eigen cereal::cereal)

  add_executable(contract_config_tests tests/contract_config_tests.cc)
  target_link_libraries(contract_config_tests GTest::gtest GTest::gtest_main)

  add_executable(contract_evolution_tests tests/contract_evolution_tests.cc nn_population.cc nn_evaluator_portable.cc)
  target_link_libraries(contract_evolution_tests GTest::gtest GTest::gtest_main Eigen3::Eigen cereal::cereal)

  add_executable(rpc_transport_tests tests/rpc_transport_tests.cc)
  target_link_libraries(rpc_transport_tests GTest::gtest GTest::gtest_main)

  enable_testing()
  add_test(NAME nn_evaluator_tests COMMAND nn_evaluator_tests)
  add_test(NAME nn_serialization_tests COMMAND nn_serialization_tests)
  add_test(NAME nn_population_tests COMMAND nn_population_tests)
  add_test(NAME fitness_computer_tests COMMAND fitness_computer_tests)
  add_test(NAME contract_evaluator_tests COMMAND contract_evaluator_tests)
  add_test(NAME contract_config_tests COMMAND contract_config_tests)
  add_test(NAME contract_evolution_tests COMMAND contract_evolution_tests)
  add_test(NAME rpc_transport_tests COMMAND rpc_transport_tests)

  # --gtest_brief=1 only prints failures, not every test run
  add_custom_target(run_autoc_tests ALL
    DEPENDS nn_evaluator_tests nn_serialization_tests nn_population_tests fitness_computer_tests
            contract_evaluator_tests contract_config_tests contract_evolution_tests rpc_transport_tests
    COMMAND nn_evaluator_tests --gtest_brief=1
    COMMAND nn_serialization_tests --gtest_brief=1
    COMMAND nn_population_tests --gtest_brief=1
    COMMAND fitness_computer_tests --gtest_brief=1
    COMMAND contract_evaluator_tests --gtest_brief=1
    COMMAND contract_config_tests --gtest_brief=1
    COMMAND contract_evolution_tests --gtest_brief=1
    COMMAND rpc_transport_tests --gtest_brief=1
  )
endif()
