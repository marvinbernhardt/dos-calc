cmake_minimum_required(VERSION 3.8)

project(dos-calc VERSION 0.1 LANGUAGES C)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

find_package(OpenMP REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(FFTW REQUIRED)  # self defined

if(NOT DEFINED ENV{GMXPREFIX})
    message(FATAL_ERROR "You must set GMXPREFIX environment variable")
endif()

add_executable(dos-calc src/dos-calc.c)

target_compile_features(dos-calc PRIVATE c_std_99)

target_include_directories(dos-calc PRIVATE src/)
target_include_directories(dos-calc PRIVATE $ENV{GMXPREFIX}/include)

target_compile_options(dos-calc PRIVATE $<$<C_COMPILER_ID:GNU>:-Wall -O3>)

execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

if(${GIT_BRANCH} STREQUAL "master")
    target_compile_definitions(dos-calc PRIVATE -DVERSION="${PROJECT_VERSION}")
else()
    target_compile_definitions(dos-calc PRIVATE -DVERSION="${PROJECT_VERSION}-${GIT_BRANCH}")
endif()

target_link_libraries(dos-calc PRIVATE OpenMP::OpenMP_C)
target_link_libraries(dos-calc PRIVATE $ENV{GMXPREFIX}/lib64/libgromacs.so)
target_link_libraries(dos-calc PRIVATE ${LAPACK_LIBRARIES})
target_link_libraries(dos-calc PRIVATE ${BLAS_LIBRARIES})
target_link_libraries(dos-calc PRIVATE ${FFTW_LIBRARIES})
target_link_libraries(dos-calc PRIVATE m gfortran)

install(TARGETS dos-calc
        RUNTIME DESTINATION bin
        )
