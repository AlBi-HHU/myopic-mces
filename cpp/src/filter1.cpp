/**
 * Created by Sarah Steinbach.
 */
#include <vector>
#include <algorithm>
#include <cmath>
#include <GraphMol/ROMol.h>
#include "atom_types.h"
#include "filter1.h"

using namespace std;
using namespace RDKit;

namespace {

struct TypedDegree {
    int type;
    double degree;
};

// (atom type, weighted degree) of all atoms, sorted by type ascending, then degree descending
vector<TypedDegree> typedDegrees(const ROMol* mol, AtomTypeIds& typeIds){
    vector<TypedDegree> res;
    res.reserve(mol->getNumAtoms());
    for (const auto atom : mol->atoms()) {
        double degree = 0.0;
        for (const auto bond : mol->atomBonds(atom)) degree += getWeightedBondDegree(bond);
        res.push_back({typeIds.get(atom), degree});
    }
    sort(res.begin(), res.end(), [](const TypedDegree& a, const TypedDegree& b) {
        return a.type != b.type ? a.type < b.type : a.degree > b.degree;
    });
    return res;
}

} // namespace

// per atom type, atoms are mapped in order of descending weighted degree;
// unmapped atoms cost their weighted degree
double applyFilter1(ROMol* mol1, ROMol* mol2) {
    AtomTypeIds typeIds;
    const vector<TypedDegree> d1 = typedDegrees(mol1, typeIds);
    const vector<TypedDegree> d2 = typedDegrees(mol2, typeIds);

    double difference = 0.0;
    size_t i = 0, j = 0;
    while (i < d1.size() && j < d2.size()) {
        if (d1[i].type < d2[j].type) difference += d1[i++].degree;
        else if (d2[j].type < d1[i].type) difference += d2[j++].degree;
        else difference += abs(d1[i++].degree - d2[j++].degree);
    }
    for (; i < d1.size(); i++) difference += d1[i].degree;
    for (; j < d2.size(); j++) difference += d2[j].degree;
    return difference / 2.0;
}
