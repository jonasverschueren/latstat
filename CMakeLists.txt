cmake_minimum_required (VERSION 3.4.2)
project (LatStat)

# Version
set (LatStat_VERSION_MAJOR 1)
set (LatStat_VERSION_MINOR 0)

# Download cubature
configure_file(CMakeLists.txt.in ${CMAKE_BINARY_DIR}/cubature-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/cubature-download )
if(result)
  message(FATAL_ERROR "CMake step for cubature failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/cubature-download )
if(result)
  message(FATAL_ERROR "Build step for cubature failed: ${result}")
endif()

add_subdirectory(${CMAKE_BINARY_DIR}/cubature-src
	${CMAKE_BINARY_DIR}/cubature-build
                 EXCLUDE_FROM_ALL)


configure_file (
	"${PROJECT_SOURCE_DIR}/LatStatConfig.h.in"
	"${PROJECT_BINARY_DIR}/LatStatConfig.h"
)

include_directories (
	"${PROJECT_BINARY_DIR}"
	"${PROJECT_SOURCE_DIR}/include"
	"${PROJECT_SOURCE_DIR}/thirdparty/eigen3.3.4"
	"${CMAKE_BINARY_DIR}/cubature-src"
)

add_subdirectory (test)
add_subdirectory (src)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O2")

add_executable(ex_main.exe ${PROJECT_SOURCE_DIR}/ex_main.cpp)
target_link_libraries( ex_main.exe
		m
		latstat
		cubature
)
