include_directories (
	"${PROJECT_SOURCE_DIR}/include"
	"${PROJECT_SOURCE_DIR}/thirdparty/eigen3.3.4"
)

file(GLOB SRC_FILES "${PROJECT_SOURCE_DIR}/src/*.cpp")

add_library (latstat ${SRC_FILES})

find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -O2")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -O2")
	set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
	target_link_libraries (latstat
		gomp
	)
else()
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
	message(WARNING "omp not found")
endif()
