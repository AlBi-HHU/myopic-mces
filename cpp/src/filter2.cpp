/**
 * Created by Sarah Steinbach.
 */
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <GraphMol/ROMol.h>
#include "atom_types.h"
#include "LAPJV.h"
#include "filter2.h"

using namespace std;
using namespace RDKit;

namespace {

// incident bond of an atom: element of the neighbouring atom and bond weight
struct Neighbor {
    int type;
    double weight;
};

// per-molecule data needed by filter 2, computed once per call
struct MolInfo {
    vector<int> type;                   // element id per atom
    vector<double> weightedDegree;      // sum of incident bond weights per atom
    vector<int> nbStart;                // neighbors of atom i: nb[nbStart[i] .. nbStart[i+1])
    vector<Neighbor> nb;                // sorted by element id ascending, then bond weight descending
    vector<vector<int>> atomsByType;    // atom indices per element id
};

MolInfo buildMolInfo(const ROMol* mol, AtomTypeIds& typeIds){
    MolInfo info;
    const unsigned int n = mol->getNumAtoms();
    info.type.resize(n);
    info.weightedDegree.assign(n, 0.0);
    info.nbStart.resize(n + 1);
    for (const auto atom : mol->atoms()) {
        // element ids are shared between both molecules of a comparison
        info.type[atom->getIdx()] = typeIds.get(atom);
    }
    info.atomsByType.resize(typeIds.size());
    for (unsigned int i = 0; i < n; i++) info.atomsByType[info.type[i]].push_back(i);

    info.nb.reserve(2 * mol->getNumBonds());
    for (unsigned int i = 0; i < n; i++) {
        info.nbStart[i] = info.nb.size();
        const Atom* atom = mol->getAtomWithIdx(i);
        for (const auto bond : mol->atomBonds(atom)) {
            const double w = getWeightedBondDegree(bond);
            info.nb.push_back({info.type[bond->getOtherAtomIdx(i)], w});
            info.weightedDegree[i] += w;
        }
        sort(info.nb.begin() + info.nbStart[i], info.nb.end(), [](const Neighbor& a, const Neighbor& b) {
            return a.type != b.type ? a.type < b.type : a.weight > b.weight;
        });
    }
    info.nbStart[n] = info.nb.size();
    return info;
}

// minimal cost for mapping the incident bonds of atom1 onto those of atom2:
// per neighbor element, bonds are mapped in order of descending weight, unmapped bonds cost their weight
double minEdgeMappingCost(const MolInfo& m1, int atom1, const MolInfo& m2, int atom2){
    const Neighbor* a = &m1.nb[m1.nbStart[atom1]];
    const Neighbor* aEnd = &m1.nb[0] + m1.nbStart[atom1 + 1];
    const Neighbor* b = &m2.nb[m2.nbStart[atom2]];
    const Neighbor* bEnd = &m2.nb[0] + m2.nbStart[atom2 + 1];
    double cost = 0.0;
    while (a != aEnd && b != bEnd) {
        if (a->type < b->type) cost += (a++)->weight;
        else if (b->type < a->type) cost += (b++)->weight;
        else cost += abs((a++)->weight - (b++)->weight);
    }
    for (; a != aEnd; a++) cost += a->weight;
    for (; b != bEnd; b++) cost += b->weight;
    return cost;
}

} // namespace

double applyFilter2(ROMol* mol1, ROMol* mol2){
    AtomTypeIds typeIds;
    MolInfo m1 = buildMolInfo(mol1, typeIds);
    MolInfo m2 = buildMolInfo(mol2, typeIds);
    m1.atomsByType.resize(typeIds.size());

    double difference = 0.0;
    vector<double> costMatrix;
    vector<int> assignment;
    for (int t = 0; t < typeIds.size(); t++) {
        const vector<int>& atoms1 = m1.atomsByType[t];
        const vector<int>& atoms2 = m2.atomsByType[t];
        const int n1 = atoms1.size(), n2 = atoms2.size();

        // atom type only in one molecule: all its bonds have to be deleted
        if (n2 == 0) {
            for (int i : atoms1) difference += m1.weightedDegree[i];
            continue;
        }
        if (n1 == 0) {
            for (int j : atoms2) difference += m2.weightedDegree[j];
            continue;
        }
        if (n1 == 1 && n2 == 1) {
            difference += minEdgeMappingCost(m1, atoms1[0], m2, atoms2[0]);
            continue;
        }

        // costMatrix[i * n + j] := minimal cost for mapping atoms1[i] onto atoms2[j];
        // pseudo atoms without incident bonds fill up the smaller side
        const int n = max(n1, n2);
        costMatrix.assign(n * n, 0.0);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double c;
                if (i < n1 && j < n2) c = minEdgeMappingCost(m1, atoms1[i], m2, atoms2[j]);
                else if (i < n1) c = m1.weightedDegree[atoms1[i]];
                else c = m2.weightedDegree[atoms2[j]];
                costMatrix[i * n + j] = c;
            }
        }
        difference += solveAssignmentLAPJV(costMatrix, n, assignment);
    }
    return difference / 2;
}
