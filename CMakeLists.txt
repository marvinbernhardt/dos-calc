cmake_minimum_required(VERSION 3.8)

project(dos-calc VERSION 0.2 LANGUAGES C CXX)

option(USE_RPATH "Link to libraries in RPATH" ON)
option(DEBUG "Turn on extensive debugging output" OFF)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

find_package(OpenMP REQUIRED)
find_package(CBLAS REQUIRED)  # self defined
find_package(LAPACKE REQUIRED) # self defined
find_package(FFTW REQUIRED)  # self defined
find_package(xdrfile REQUIRED)

add_executable(dos-calc src/dos-calc.c)
target_include_directories(dos-calc PRIVATE src/)
target_compile_features(dos-calc PUBLIC cxx_std_11)
target_compile_options(dos-calc PRIVATE $<$<C_COMPILER_ID:GNU>:-Wall -O3>)
if(USE_RPATH)
    set_target_properties(dos-calc PROPERTIES
        BUILD_WITH_INSTALL_RPATH TRUE
        INSTALL_RPATH_USE_LINK_PATH TRUE)
endif(USE_RPATH)
if(DEBUG)
    target_compile_definitions(dos-calc PRIVATE -DDEBUG="1")
endif(DEBUG)

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

target_include_directories(dos-calc PUBLIC ${xdrfile_INCLUDE_DIRS})
target_link_libraries(dos-calc PRIVATE OpenMP::OpenMP_C)
target_link_libraries(dos-calc PRIVATE ${CBLAS_LIBRARIES})
target_link_libraries(dos-calc PRIVATE ${LAPACKE_LIBRARIES})
target_link_libraries(dos-calc PRIVATE ${FFTW_LIBRARIES})
target_link_libraries(dos-calc PRIVATE m gfortran)
target_link_libraries(dos-calc PUBLIC ${xdrfile_LIBRARIES})

install(TARGETS dos-calc
        RUNTIME DESTINATION bin)
