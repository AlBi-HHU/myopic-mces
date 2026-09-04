#include <boost/python.hpp>
#include <iostream>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include "filter0.h"

extern double applyFilter0(RDKit::ROMol* mol1, RDKit::ROMol* mol2);
extern double applyFilter1(RDKit::ROMol* mol1, RDKit::ROMol* mol2);
extern double applyFilter2(RDKit::ROMol* mol1, RDKit::ROMol* mol2);

// boost bindings for all three filters
BOOST_PYTHON_MODULE(mmces_filters) {
    using namespace boost::python;
    
    // map the Python name "filterX" to C++ function "applyFilterX"
    def("filter0", applyFilter0, (boost::python::arg("mol1"), boost::python::arg("mol2")), "Calculates filter 0 for two molecules (atom type difference)");
    def("filter1", applyFilter1, (boost::python::arg("mol1"), boost::python::arg("mol2")), "Calculates filter 1 for two molecules (comparison of bond degree by partition)");
    def("filter2", applyFilter2, (boost::python::arg("mol1"), boost::python::arg("mol2")), "Calculates filter 2 for two molecules (bipartite graph matching via hungarian algorithm)");
}
