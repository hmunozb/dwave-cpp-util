cmake_minimum_required(VERSION 3.9)
project(dw_examples)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
include_directories(include)

find_package(dwave_sapi REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options filesystem)
include_directories(${dwave_sapi_INCLUDE_DIR})

add_library(dwave_cpp SHARED
        src/dwave_cpp/dwave_cpp_core.cpp
        src/dwave_cpp/dwave_cpp_problems.cpp
        src/dwave_cpp/dwave_cpp_schedules.cpp)

target_link_libraries(dwave_cpp ${dwave_sapi_LIBRARY})

add_executable(sapi_connect_example c_examples/sapi_connectingToSolver.c)
target_link_libraries(sapi_connect_example ${dwave_sapi_LIBRARY})

add_executable(connection_example test/examples/connectingToSolver.cpp)
target_link_libraries(connection_example ${dwave_sapi_LIBRARY} dwave_cpp)

add_executable(dwave_cpp_test test/test_dwave_cpp_problems.cpp)
target_link_libraries(dwave_cpp_test ${dwave_sapi_LIBRARY} dwave_cpp)

add_executable(bound_canc_0989 src/dw_prog.cpp src/bound_canc.cpp)
target_link_libraries(bound_canc_0989 ${dwave_sapi_LIBRARY} dwave_cpp ${Boost_LIBRARIES})

add_executable(chain_qac src/dw_prog.cpp src/chain_qac.cpp)
target_link_libraries(bound_canc_0989 ${dwave_sapi_LIBRARY} dwave_cpp ${Boost_LIBRARIES})