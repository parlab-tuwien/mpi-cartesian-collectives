
cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules/")
   
project(tuwcartesians)

set(SRC_DIR "${CMAKE_SOURCE_DIR}/src")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}")


set(INCLUDE_PLATFORM_CONFIG_FILE "${CMAKE_SOURCE_DIR}/platform_files/default.cmake" 
    CACHE STRING "Configure project to use a specific platform file")
include(${INCLUDE_PLATFORM_CONFIG_FILE})

#find_package(GSL REQUIRED)
#if (GSL_INCLUDE_DIR)
#message (STATUS "GSL INCLUDES: ${GSL_INCLUDE_DIR}")
#else(GSL_INCLUDE_DIR)
#message (FATAL_ERROR "GSL libraries not found.")
#endif(GSL_INCLUDE_DIR)

#INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIR} "src")
#set(COMMON_LIBRARIES ${GSL_LIBRARIES} )


SET(MY_COMPILE_FLAGS "")


set(BENCH_COMMON_SRC_FILES
${SRC_DIR}/tuwcartesians.c
${SRC_DIR}/tuw_bench_setup.c
${SRC_DIR}/tuw_bench_arg_parse.c
)

add_executable(bench_cartesians
${BENCH_COMMON_SRC_FILES}
${SRC_DIR}/tuw_bench_runner_full.c
${SRC_DIR}/bench_cartesians.c
)

add_executable(check_cartesians
${BENCH_COMMON_SRC_FILES}
${SRC_DIR}/check_cartesians.c
)
        
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_SOURCE_GENERATOR "TBZ2")


set(BENCHMARK_NAME ${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH})
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${BENCHMARK_NAME}-src")
set(DIST_ARCHIVE "${BENCHMARK_NAME}-dist")

SET(CPACK_SOURCE_IGNORE_FILES
"/CMakeFiles/"
"/_CPack_Packages/"
"/bin/"
"/include/"
"/lib/"
"/scripts/"
"/create_package.sh"
"/Makefile"
"\\\\.tar\\\\.bz2$"
"/\\\\.git/"
"/\\\\.settings/"
"/\\\\..*"
"CMakeCache.txt"
"cmake_install\\\\.cmake"
"CPackConfig\\\\.cmake"
"CPackSourceConfig\\\\.cmake"
"sources\\\\.org"
"${CPACK_SOURCE_IGNORE_FILES}"
)

INCLUDE(CPack)

add_custom_target(pack COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
