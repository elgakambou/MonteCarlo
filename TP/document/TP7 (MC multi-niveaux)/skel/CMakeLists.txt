cmake_minimum_required(VERSION 2.8)
project(Discretization CXX)

# Release or Debug
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Debug' as none was specified.")
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
endif ()

set(CMAKE_PREFIX_PATH /home/elga/Documents/ENSIMAG/MONTE-CARLO_SIMULATION-DE-LOI/MonteCarlo/pnl/build ${CMAKE_PREFIX_PATH})
# Detect PNL
find_package(Pnl REQUIRED)
set(LIBS ${LIBS} ${PNL_LIBRARIES})
include_directories(${PNL_INCLUDE_DIRS})
if(${CMAKE_BUILD_TYPE} STREQUAL "Release")
    add_definitions(-DPNL_RANGE_CHECK_OFF)
endif()

# Testing the parser
set(MC_SOURCES AsianOption.cpp  Heston.cpp Model.cpp)
set(MC_HEADERS AsianOption.hpp  Heston.hpp  Model.hpp  Option.hpp)
add_executable(mc ${MC_SOURCES} ${MC_HEADERS} MonteCarlo.cpp MonteCarlo.hpp mcpricer.cpp)
target_link_libraries(mc ${LIBS})
pnl_add_postbuild(mc)
add_executable(mlmc ${MC_SOURCES} ${MC_HEADERS} MultiLevelMonteCarlo.cpp MultiLevelMonteCarlo.hpp mlmcpricer.cpp)
target_link_libraries(mlmc ${LIBS})
pnl_add_postbuild(mlmc)

# Print compiler flags
get_directory_property(DirDefs COMPILE_DEFINITIONS)
message(STATUS "COMPILE_DEFINITIONS = ${DirDefs}")
if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
    message(STATUS "C++ Flags: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_LINKER_EXE_FLAGS}" )
else()
    message(STATUS "C++ Flags: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_LINKER_EXE_FLAGS}" )
endif()

