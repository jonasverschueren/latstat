# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in ${CMAKE_BINARY_DIR}/googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
	${CMAKE_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
  include_directories("${gtest_SOURCE_DIR}/include")
endif()

include_directories (
	${PROJECT_SOURCE_DIR}/include
	${PROJECT_SOURCE_DIR}/test
	"${PROJECT_SOURCE_DIR}/thirdparty/eigen3.3.4"
)

add_executable(test_suite.exe ${PROJECT_SOURCE_DIR}/test/test_suite.cpp)

target_link_libraries(test_suite.exe 
	m
	gtest_main
	latstat
)

# Copy potential files to build directory
file(GLOB POTENTIAL_FILES "${PROJECT_SOURCE_DIR}/pot/*")
file(COPY ${POTENTIAL_FILES} DESTINATION ${CMAKE_BINARY_DIR})
#configure_file(${PROJECT_SOURCE_DIR}/pot/* ${CMAKE_BINARY_DIR} COPYONLY)

set_target_properties(test_suite.exe PROPERTIES LINK_FLAGS -Wl,-lpthread)

find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -O2")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -O2")
	set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
	message(WARNING "omp not found")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O2")
endif()
