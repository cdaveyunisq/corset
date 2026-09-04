cmake_minimum_required(VERSION 3.16)

# Define your project name along with its modern version numbers
project(CorsetIncremental 
        VERSION 2.0.0 
        LANGUAGES CXX C)

project(corset LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

# Optional: add custom debug flags
if(CMAKE_BUILD_TYPE MATCHES Debug)
    if(MSVC)
        add_compile_options(/Zi /Od)
    else()
    # (-fsanitize=address 
        add_compile_options(-g -O0 -Wall -Wextra -gz=none )
    endif()
endif()




# ---------------------------------------------------------------------------
# Configure the version header file
# ---------------------------------------------------------------------------
# This reads version.h.in and generates version.h inside the build folder
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/src/version.h.in"
    "${CMAKE_CURRENT_BINARY_DIR}/generated/version.h"
)


# Include the External Project module to manage non-CMake dependencies
include(ExternalProject)

# ---------------------------------------------------------------------------
# 1. System Requirements Check (Samtools / HTSlib prerequisites)
# ---------------------------------------------------------------------------
find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(LibLZMA REQUIRED)
find_package(CURL REQUIRED) 
find_package(OpenSSL REQUIRED)

# 2. Path Directives
set(DEPS_DIR "${CMAKE_BINARY_DIR}/thirdparty")
set(INSTALL_PREFIX "${DEPS_DIR}/install")

file(MAKE_DIRECTORY "${INSTALL_PREFIX}/include")
file(MAKE_DIRECTORY "${INSTALL_PREFIX}/lib")

# ---------------------------------------------------------------------------
# 3. External Git Project Declarations (Tag 1.24)
# ---------------------------------------------------------------------------

# A. HTSLIB (The underlying core engine required by both tools)
ExternalProject_Add(htslib_git_target
    GIT_REPOSITORY      "https://github.com/samtools/htslib.git"
    GIT_TAG             "1.24"
    PREFIX              "${DEPS_DIR}/htslib"
   
    # FIX: Native list-wrapping forces the host shell to correctly execute the processing pipeline
    CONFIGURE_COMMAND   /bin/sh -c "autoheader && autoconf && ./configure --prefix=${INSTALL_PREFIX} --disable-plugins"
    
    BUILD_COMMAND       make lib-static
    INSTALL_COMMAND     make install
    
    BUILD_IN_SOURCE     1
    LOG_DOWNLOAD        1
    LOG_CONFIGURE       1 # Enable configuration logs to trace errors easily
    LOG_BUILD           1
)


# 5. Export HTSlib Interface Targets to Corset
set(HTSLIB_STATIC_LIB "${INSTALL_PREFIX}/lib/libhts.a")
set_source_files_properties("${HTSLIB_STATIC_LIB}" PROPERTIES GENERATED TRUE)

add_library(HtslibStatic STATIC IMPORTED GLOBAL)
set_target_properties(HtslibStatic PROPERTIES
    IMPORTED_LOCATION "${HTSLIB_STATIC_LIB}"
    # CRUCIAL FIX: Expose both the base include directory AND the nested htslib directory
    INTERFACE_INCLUDE_DIRECTORIES "${INSTALL_PREFIX}/include;${INSTALL_PREFIX}/include/htslib"
)
add_dependencies(HtslibStatic htslib_git_target)


# ---------------------------------------------------------------------------
# 4. Define Your Corset Executable
# ---------------------------------------------------------------------------
# Explicitly track your new partitioned/incremental source files
add_executable(corset 
    src/corset.cc
    src/Cluster.cc
    src/MakeClusters.cc
    src/Read.cc
    src/Transcript.cc)

add_executable(corset_fasta_ID_changer
    src/corset_fasta_ID_changer.cc)

target_include_directories(corset PRIVATE 
    src
    "${CMAKE_CURRENT_BINARY_DIR}/generated")


# Link your executable against HTSlib and core system threading/compression architectures
target_link_libraries(corset PRIVATE 
    HtslibStatic
    CURL::libcurl
    OpenSSL::Crypto 
    ZLIB::ZLIB
    BZip2::BZip2
    LibLZMA::LibLZMA
    pthread
)

target_include_directories(corset_fasta_ID_changer PRIVATE 
    src
    "${CMAKE_CURRENT_BINARY_DIR}/generated"
)

# Even if this utility doesn't read BAMs, linking it to the same group is safe,
# or you can streamline it to only use standard libraries if it is pure C++
target_link_libraries(corset_fasta_ID_changer PRIVATE 
    pthread
)