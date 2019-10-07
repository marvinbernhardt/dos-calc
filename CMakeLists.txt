cmake_minimum_required(VERSION 3.8)

project(dos-calc VERSION 0.3 LANGUAGES C CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(USE_RPATH "Link to libraries in RPATH" ON)
option(DEBUG "Turn on extensive debugging output at runtime" OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(OpenMP REQUIRED)
find_package(CBLAS REQUIRED)  # self defined
find_package(LAPACKE REQUIRED) # self defined
find_package(FFTW REQUIRED)  # self defined
find_package(xdrfile REQUIRED)
find_package(cJSON REQUIRED)

add_executable(dos-calc src/dos-calc.c)
target_sources(dos-calc PRIVATE src/)
target_include_directories(dos-calc PRIVATE include/)
target_compile_options(dos-calc PUBLIC "$<$<CONFIG:RELEASE>:-O3>")
target_compile_options(dos-calc PRIVATE $<$<C_COMPILER_ID:GNU>:-Wall -Wextra -Wpedantic>)
if(USE_RPATH)
    set_target_properties(dos-calc PROPERTIES
        BUILD_WITH_INSTALL_RPATH TRUE
        INSTALL_RPATH_USE_LINK_PATH TRUE)
endif(USE_RPATH)
if(DEBUG)
    target_compile_definitions(dos-calc PRIVATE -DDEBUG="1")
endif(DEBUG)

target_compile_definitions(dos-calc PRIVATE -DVERSION="${PROJECT_VERSION}")

target_include_directories(dos-calc PUBLIC ${xdrfile_INCLUDE_DIRS})
target_include_directories(dos-calc PUBLIC ${CJSON_INCLUDE_DIRS})
target_link_libraries(dos-calc PRIVATE OpenMP::OpenMP_C)
target_link_libraries(dos-calc PRIVATE ${CBLAS_LIBRARIES})
target_link_libraries(dos-calc PRIVATE ${LAPACKE_LIBRARIES})
target_link_libraries(dos-calc PRIVATE ${FFTW_LIBRARIES})
target_link_libraries(dos-calc PRIVATE m gfortran)
target_link_libraries(dos-calc PUBLIC ${xdrfile_LIBRARIES})
target_link_libraries(dos-calc PUBLIC ${CJSON_LIBRARIES})

install(TARGETS dos-calc
        RUNTIME DESTINATION bin)
