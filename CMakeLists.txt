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
    src/validation.cpp
    src/state_io.cpp
)
target_include_directories(gravityx_go2 PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_compile_definitions(gravityx_go2 PRIVATE USE_IPOPT)
target_link_libraries(gravityx_go2 PRIVATE
    gravity_pinned
    "${IPOPT_ROOT}/lib/libipopt.so"
    nlohmann_json::nlohmann_json
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
