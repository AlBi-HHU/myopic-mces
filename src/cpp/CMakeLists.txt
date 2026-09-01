cmake_minimum_required(VERSION 3.16)
project(mm_filters LANGUAGES CXX)
cmake_policy(SET CMP0167 NEW)
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static-libstdc++ -static-libgcc")

include(GNUInstallDirs)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

# folder inclusion paths for conda env on all architectures

if(WIN32)
    set(CONDA_INCLUDE_DIR "$ENV{CONDA_PREFIX}/Library/include")
    set(CONDA_LIB_DIR     "$ENV{CONDA_PREFIX}/Library/lib")
else()
    set(CONDA_INCLUDE_DIR "$ENV{CONDA_PREFIX}/include")
    set(CONDA_LIB_DIR     "$ENV{CONDA_PREFIX}/lib")
endif()

list(APPEND CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}")

set(BOOST_PY_COMPONENT "python${Python_VERSION_MAJOR}${Python_VERSION_MINOR}")
find_package(Boost REQUIRED COMPONENTS ${BOOST_PY_COMPONENT})

add_library(mmces_filters SHARED
    bindings.cpp
    filter0.cpp
    filter1.cpp
    filter2.cpp
    Hungarian.cpp
)

target_include_directories(mmces_filters PRIVATE 
        $ENV{CONDA_PREFIX}/include
        $ENV{CONDA_PREFIX}/include/rdkit
    )
    target_link_directories(mmces_filters PRIVATE 
        $ENV{CONDA_PREFIX}/lib
    )

# link binaries
target_link_libraries(mmces_filters PRIVATE 
    Boost::${BOOST_PY_COMPONENT}
    Python::Module
    RDKitGraphMol
    RDKitDataStructs
    RDKitSmilesParse
)

#
target_link_options(mmces_filters PRIVATE
    -static-libstdc++
    -static-libgcc
)

set_target_properties(mmces_filters PROPERTIES PREFIX "")
if(WIN32)
    set_target_properties(mmces_filters PROPERTIES SUFFIX ".pyd")
endif()

# makes sure compiled .so file is in correct place for python wheel
install(TARGETS mmces_filters DESTINATION .)


