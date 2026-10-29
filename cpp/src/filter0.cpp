/**
 * Created by Sarah Steinbach.
 */
#include <vector>
#include <cstdlib>
#include <GraphMol/ROMol.h>
#include "atom_types.h"
#include "filter0.h"

using namespace std;
using namespace RDKit;

// difference of the atom type counts (hydrogens are ignored)
double applyFilter0(ROMol* mol1, ROMol* mol2){
    AtomTypeIds typeIds;
    vector<int> counts;   // count in mol1 minus count in mol2, per atom type
    auto addCounts = [&](const ROMol* mol, int sign) {
        for (const auto atom : mol->atoms()) {
            if (atom->getAtomicNum() == 1) continue;
            const int t = typeIds.get(atom);
            if (t >= (int)counts.size()) counts.resize(t + 1, 0);
            counts[t] += sign;
        }
    };
    addCounts(mol1, 1);
    addCounts(mol2, -1);
    double difference = 0.0;
    for (int c : counts) difference += abs(c);
    return difference;
}
