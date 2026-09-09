cmake_minimum_required(VERSION 3.16)
project(mm_filters)

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(FAST_BUILD "Testbuild" OFF)

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

# module target
add_library(mmces_filters SHARED 
            bindings.cpp
            filter0.cpp
            filter1.cpp
            filter2.cpp
            Hungarian.cpp)

set_target_properties(mmces_filters PROPERTIES PREFIX "")
if(WIN32)
    set_target_properties(mmces_filters PROPERTIES SUFFIX ".pyd")
endif()

if (FAST_BUILD)

    # folder inclusion paths
    if(DEFINED ENV{CONDA_PREFIX})
        message(STATUS "CONDA detected")
        target_include_directories(mmces_filters PRIVATE 
            $ENV{CONDA_PREFIX}/include
            $ENV{CONDA_PREFIX}/include/rdkit
        )
        target_link_directories(mmces_filters PRIVATE 
            $ENV{CONDA_PREFIX}/lib
        )
    endif()

    # link binaries
    target_link_libraries(mmces_filters PRIVATE 
        boost_python${Python_VERSION_MAJOR}${Python_VERSION_MINOR}
        Python::Module
        RDKitGraphMol
        RDKitDataStructs
        RDKitSmilesParse
    )

else ()
    cmake_policy(SET CMP0167 NEW)
    
    set(RDKIT_LIB_DIR "$ENV{RDKIT_LIB_DIR}")
    set(RDKIT_INCLUDE_DIR "$ENV{RDKIT_INCLUDE_DIR}")
    set(RDKIT_BUILD_INCLUDE_DIR "$ENV{RDKIT_BUILD_INCLUDE_DIR}")

    # longer build for multiple architectures and Python versions
    # includes custom 
    file(GLOB BOOST_PYTHON_LIBRARY
        "$ENV{CONDA_PREFIX}/libboost_python${Python_VERSION_MAJOR}${Python_VERSION_MINOR}-*.so.*"
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

    set(RDKIT_INCLUDE_DIR "$ENV{RDKIT_INCLUDE_DIR}")
    set(RDKIT_BUILD_INCLUDE_DIR "$ENV{RDKIT_BUILD_INCLUDE_DIR}")

    message(STATUS "RDKit include directory: ${RDKIT_INCLUDE_DIR}")
    message(STATUS "RDKit build include directory: ${RDKIT_BUILD_INCLUDE_DIR}")
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

endif()

# makes sure compiled .so file is in correct place for python wheel
install(TARGETS mmces_filters DESTINATION .)
