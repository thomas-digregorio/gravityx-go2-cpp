cmake_minimum_required(VERSION 3.20)
project(gravityx_go2_cpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(GRAVITY_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/.deps/Gravity" CACHE PATH "Pinned Gravity source tree")
set(IPOPT_ROOT "$ENV{IPOPT_ROOT_DIR}" CACHE PATH "Ipopt installation prefix")

if(NOT EXISTS "${GRAVITY_ROOT}/include/gravity/model.h")
    message(FATAL_ERROR "Gravity was not found at ${GRAVITY_ROOT}; run scripts/bootstrap_wsl.sh")
endif()
if(NOT EXISTS "${GRAVITY_ROOT}/lib/libgravity.a")
    message(FATAL_ERROR "Pinned Gravity library is not built; run scripts/bootstrap_wsl.sh")
endif()
if(NOT EXISTS "${IPOPT_ROOT}/lib/libipopt.so")
    message(FATAL_ERROR "Ipopt was not found below IPOPT_ROOT=${IPOPT_ROOT}")
endif()

find_package(nlohmann_json CONFIG REQUIRED HINTS "${IPOPT_ROOT}")
find_package(Eigen3 CONFIG REQUIRED HINTS "${IPOPT_ROOT}")
find_path(HIGHS_INCLUDE_DIR highs/Highs.h HINTS "${IPOPT_ROOT}/include" REQUIRED)
find_library(HIGHS_LIBRARY highs HINTS "${IPOPT_ROOT}/lib" REQUIRED)

add_library(gravity_pinned STATIC IMPORTED GLOBAL)
set_target_properties(gravity_pinned PROPERTIES
    IMPORTED_LOCATION "${GRAVITY_ROOT}/lib/libgravity.a"
    INTERFACE_INCLUDE_DIRECTORIES "${GRAVITY_ROOT}/include;${IPOPT_ROOT}/include/coin;${IPOPT_ROOT}/include"
)

add_executable(gravityx_go2
    src/main.cpp
    src/case_data.cpp
    src/ac_model.cpp
    src/algorithm.cpp
    src/fast_power_flow.cpp
    src/linearized_ac_seed.cpp
    src/validation.cpp
    src/state_io.cpp
)
target_include_directories(gravityx_go2 PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${HIGHS_INCLUDE_DIR}"
    "${HIGHS_INCLUDE_DIR}/highs"
)
target_compile_definitions(gravityx_go2 PRIVATE USE_IPOPT)
target_link_libraries(gravityx_go2 PRIVATE
    gravity_pinned
    "${IPOPT_ROOT}/lib/libipopt.so"
    Eigen3::Eigen
    nlohmann_json::nlohmann_json
    "${HIGHS_LIBRARY}"
    m
    pthread
)
set_target_properties(gravityx_go2 PROPERTIES
    BUILD_RPATH "${IPOPT_ROOT}/lib"
)

enable_testing()
find_package(Python3 COMPONENTS Interpreter REQUIRED)
add_test(NAME gravity_smoke COMMAND gravityx_go2 smoke)
add_test(NAME gravity_component_tests COMMAND gravityx_go2 component-tests)
add_test(NAME gravity_parallel_circuit_regression COMMAND gravityx_go2 parallel-circuit-test)
add_test(
    NAME solution_writer_tests
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_solution_writer.py"
)
