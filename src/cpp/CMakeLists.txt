cmake_minimum_required(VERSION 3.16)
project(mm_filters LANGUAGES CXX)

cmake_policy(SET CMP0167 NEW)

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(RDKIT_INCLUDE_DIR "$ENV{RDKIT_INCLUDE_DIR}")
set(RDKIT_LIBRARY_DIR "$ENV{RDKIT_LIBRARY_DIR}")

message(STATUS "RDKit include directory: ${RDKIT_INCLUDE_DIR}")
message(STATUS "RDKit library directory: ${RDKIT_LIBRARY_DIR}")


find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)


# Boost
list(APPEND CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}")

set(BOOST_PY_COMPONENT
    "python${Python_VERSION_MAJOR}${Python_VERSION_MINOR}"
)

find_package(Boost REQUIRED COMPONENTS ${BOOST_PY_COMPONENT})

# Library
add_library(mmces_filters SHARED
    bindings.cpp
    filter0.cpp
    filter1.cpp
    filter2.cpp
    Hungarian.cpp
)

# RDKit
target_include_directories(mmces_filters PRIVATE
    "${RDKIT_INCLUDE_DIR}"
    "${RDKIT_BUILD_INCLUDE_DIR}"
)
target_link_directories(mmces_filters PRIVATE
    "${RDKIT_LIBRARY_DIR}"
)

# Link
target_link_libraries(mmces_filters PRIVATE
    Boost::${BOOST_PY_COMPONENT}
    Python::Module
    RDKitGraphMol
    RDKitDataStructs
    RDKitSmilesParse
)

target_link_options(mmces_filters PRIVATE
    -static-libstdc++
    -static-libgcc
)

set_target_properties(mmces_filters PROPERTIES PREFIX "")

if(WIN32)
    set_target_properties(mmces_filters PROPERTIES SUFFIX ".pyd")
endif()

install(TARGETS mmces_filters DESTINATION .)
