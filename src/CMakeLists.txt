FetchContent_Declare(pybind11
  URL https://github.com/pybind/pybind11/archive/v2.6.2.tar.gz
)

FetchContent_MakeAvailable(pybind11)


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
target_precompile_headers(
			mormon
			PRIVATE
			pch.h)

include(FortranCInterface)
FortranCInterface_VERIFY(CXX)



target_link_libraries(
    mormon
    PRIVATE
    OpenMP::OpenMP_CXX
    Eigen3::Eigen
    xdr
    random
    fmt::fmt
)
