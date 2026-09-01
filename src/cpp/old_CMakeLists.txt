cmake_minimum_required(VERSION 3.16)
project(mm_filters)

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Python COMPONENTS Interpreter Development REQUIRED)

# module target
add_library(mmces_filters SHARED 
            bindings.cpp
            filter0.cpp
            filter1.cpp
            filter2.cpp
            Hungarian.cpp)

# folder inclusion paths
if(DEFINED ENV{CONDA_PREFIX})
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
    boost_python310
    Python::Python
    RDKitGraphMol
    RDKitDataStructs
    RDKitSmilesParse
)

set_target_properties(mmces_filters PROPERTIES PREFIX "")

# makes sure compiled .so file is in correct place for python wheel
install(TARGETS mmces_filters DESTINATION .)
