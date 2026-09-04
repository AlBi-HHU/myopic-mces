cmake_minimum_required(VERSION 3.16)
project(mm_filters LANGUAGES CXX)

cmake_policy(SET CMP0167 NEW)

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# paths 
set(RDKIT_LIB_DIR "$ENV{RDKIT_LIB_DIR}")

message(STATUS "RDKit library directory: ${RDKIT_LIB_DIR}")

###################### Python ######################

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

##################### Library searching (and hopefully finding) #####################

file(GLOB BOOST_PYTHON_LIBRARY
    "$ENV/lib/libboost_python${Python_VERSION_MAJOR}${Python_VERSION_MINOR}-*.so.*"
)

file(GLOB RDKIT_GRAPHMOL_LIBRARY
    "${RDKIT_LIB_DIR}/libRDKitGraphMol-*.so.*"
)

file(GLOB RDKIT_DATASTRUCTS_LIBRARY
    "${RDKIT_LIB_DIR}/libRDKitDataStructs-*.so.*"
)

file(GLOB RDKIT_SMILESPARSE_LIBRARY
    "${RDKIT_LIB_DIR}/libRDKitSmilesParse-*.so.*"
)

message(STATUS "Boost.Python: ${BOOST_PYTHON_LIBRARY}")
message(STATUS "RDKit graphmol: ${RDKIT_GRAPHMOL_LIBRARY}")
message(STATUS "RDKit datastructs: ${RDKIT_DATASTRUCTS_LIBRARY}")
message(STATUS "RDKit smilesparse: ${RDKIT_SMILESPARSE_LIBRARY}")

###################### RDKit ######################

set(RDKIT_INCLUDE_DIR "$ENV{RDKIT_INCLUDE_DIR}")
set(RDKIT_BUILD_INCLUDE_DIR "$ENV{RDKIT_BUILD_INCLUDE_DIR}")

message(STATUS "RDKit include directory: ${RDKIT_INCLUDE_DIR}")
message(STATUS "RDKit build include directory: ${RDKIT_BUILD_INCLUDE_DIR}")

###################### Library ######################
add_library(mmces_filters SHARED
    bindings.cpp
    filter0.cpp
    filter1.cpp
    filter2.cpp
    Hungarian.cpp
)

target_compile_definitions(mmces_filters PRIVATE
    _GLIBCXX_USE_CXX11_ABI=0
)

target_include_directories(mmces_filters PRIVATE
    "${RDKIT_INCLUDE_DIR}"
    "${RDKIT_BUILD_INCLUDE_DIR}"
    "${RDKIT_LIB_DIR}"
    "$ENV{CONDA_PREFIX}/include"
)

target_link_directories(mmces_filters PRIVATE
    "${RDKIT_LIB_DIR}"
)

target_link_libraries(mmces_filters PRIVATE
    "${BOOST_PYTHON_LIBRARY}"
    Python::Module
    "${RDKIT_DATASTRUCTS_LIBRARY}"
    "${RDKIT_GRAPHMOL_LIBRARY}"
    "${RDKIT_SMILESPARSE_LIBRARY}"
)

set_target_properties(mmces_filters PROPERTIES PREFIX "")

if(WIN32)
    set_target_properties(mmces_filters PROPERTIES SUFFIX ".pyd")
endif()

###################### Installation ######################

install(TARGETS mmces_filters DESTINATION .)
