add_library(
	random
	STATIC
	ran2.f
	random.cpp
	ran2_lib.cpp
)

target_compile_options(
	random
	PRIVATE
	-O3
	-ffast-math
)

target_link_options(
	random
	PRIVATE
	-O3
	-ffast-math
	)

target_link_libraries(random 
                    PRIVATE 
                    Eigen3::Eigen)