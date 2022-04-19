FetchContent_Declare(pybind11
  URL https://github.com/pybind/pybind11/archive/v2.6.2.tar.gz
)



if(${CMAKE_VERSION} VERSION_LESS "3.14.0") 
	FetchContent_GetProperties(pybind11)
	if(NOT pybind11_POPULATED)
		FetchContent_Populate(pybind11)
		add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR} EXCLUDE_FROM_ALL)
	endif()
else()
	FetchContent_MakeAvailable(pybind11)
endif()

pybind11_add_module(mormon)

target_compile_options(
	mormon
        PRIVATE
	-std=c++20
        -O3
        -ffast-math
	-march=native
	-D_DEBUG_
	-DTIMERS
	-Wall
)
target_link_options(
	mormon
	PRIVATE
	-lstdc++
	-O3
	-ffast-math
)

add_definitions(-DPY11)

target_sources(
	mormon
	PRIVATE
	io.cpp
	sampler.cpp
	Faddeeva.cpp
	mormon.cpp
)
target_include_directories(mormon
			PRIVATE
			${CMAKE_SOURCE_DIR}/src)

if(${CMAKE_VERSION} VERSION_LESS "3.16.0") 
	message("NOT using precompiled header, CMAKE version too old (< 3.16)")
endif()

target_precompile_headers(
			mormon
			PRIVATE
			pch.h)

include(FortranCInterface)
FortranCInterface_VERIFY(CXX)


if(APPLE)
	target_link_libraries(
		mormon
		PRIVATE
		OpenMP::OpenMP_CXX
		Eigen3::Eigen
		xdr
		random
		fmt::fmt
	)
else()
	target_link_libraries(
		mormon
		PRIVATE
		OpenMP::OpenMP_CXX
		xdr
		random
		fmt::fmt
	)
endif()

