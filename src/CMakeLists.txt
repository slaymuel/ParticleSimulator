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

pybind11_add_module(particlesimulator)

if(APPLE)
	target_compile_options(
		particlesimulator
			PRIVATE
		-std=c++20
			-O3
			-ffast-math
		-march=native
		-D_DEBUG_
		-D_TIMERS_
		-Wall
	)
	target_link_options(
		particlesimulator
		PRIVATE
		-lstdc++
		-O3
		-ffast-math
	)
else()
target_compile_options(
	particlesimulator
        PRIVATE
	-std=c++2a
        -O3
        -ffast-math
	-march=native
	-D_DEBUG_
	-D_TIMERS_
	-Wall
)
target_link_options(
	particlesimulator
	PRIVATE
	-lstdc++
	-O3
	-ffast-math
	-fopenmp=libomp
)
endif()


add_definitions(-DPY11)

target_sources(
	particlesimulator
	PRIVATE
	io.cpp
	sampler.cpp
	move.cpp
	state.cpp
	Faddeeva.cpp
	constants.cpp
	particlesimulator.cpp
)
target_include_directories(particlesimulator
			PRIVATE
			${CMAKE_SOURCE_DIR}/src)

if(${CMAKE_VERSION} VERSION_LESS "3.16.0") 
	message("NOT using precompiled header, CMAKE version too old (< 3.16)")
else()
	target_precompile_headers(
		particlesimulator
		PRIVATE
		pch.h
	)
endif()



include(FortranCInterface)
FortranCInterface_VERIFY(CXX)


if(APPLE)
	target_link_libraries(
		particlesimulator
		PRIVATE
		OpenMP::OpenMP_CXX
		Eigen3::Eigen
		xdr
		random
		fmt::fmt
	)
else()
	target_link_libraries(
		particlesimulator
		PRIVATE
		xdr
		random
		fmt::fmt
	)
endif()

