cmake_minimum_required(VERSION 3.13 FATAL_ERROR)

include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

project(mormon LANGUAGES CXX C Fortran)

set(CMAKE_CXX_STANDARD 20)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Python 3.6 REQUIRED COMPONENTS Interpreter Development)
find_package(OpenMP REQUIRED COMPONENTS CXX)

Include(FetchContent)
include(CTest)

enable_testing()

add_subdirectory(lib/fmt)
add_subdirectory(src/xdr)
add_subdirectory(src/random)
add_subdirectory(src)

#################################       CATCH2 UNIT TESTING    ##################################
#FetchContent_Declare(
#  Catch2
#  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
#  GIT_TAG        v3.0.0-preview3
#)

#FetchContent_MakeAvailable(Catch2)
add_subdirectory(lib/Catch2)

add_executable(particles_test tests/particles_test.cpp)
target_link_libraries(particles_test 
                      PRIVATE
                      Eigen3::Eigen
                      Catch2::Catch2WithMain
                      random
                      )
#target_link_options(particles_test
#                    PRIVATE
#                    -force_load Catch2)
#list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
#list(APPEND CMAKE_MODULE_PATH lib/Catch2/extras)

include(lib/Catch2/extras/Catch.cmake)
catch_discover_tests(particles_test
                    REPORTER console)
#################################################################################################



#################################       INTEGRATION TESTS    ##################################
add_executable(ran2_test ${CMAKE_CURRENT_LIST_DIR}/tests/ran2_test.cpp)

target_link_libraries(ran2_test 
                    PRIVATE
                    random
                    Eigen3::Eigen)

add_test(
	NAME python_ewald3d_energy_test
	COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/ewald3d_energy_test.py --checkpoint ${CMAKE_CURRENT_SOURCE_DIR}/tests/ewald_test.cp
)

add_test(
        NAME python_halfwald_energy_test
        COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/hw_energy_test.py --checkpoint ${CMAKE_CURRENT_SOURCE_DIR}/tests/hw_energy_test.cp
)

add_test(
  NAME ran2_test
  COMMAND $<TARGET_FILE:ran2_test>
)
################################################################################################
