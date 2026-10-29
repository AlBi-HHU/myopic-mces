// Helpers shared by the filters.
#ifndef ATOM_TYPES_H
#define ATOM_TYPES_H

#include <string>
#include <utility>
#include <vector>
#include <GraphMol/ROMol.h>

// Assigns compact ids 0, 1, 2, ... to the atom types (element symbols) occurring in the molecules
// of one comparison. Equivalent to comparing Atom::getSymbol(), but avoids its string construction:
// the symbol is determined by the atomic number, except for dummy atoms (atomic number 0) which can
// carry a label.
class AtomTypeIds {
public:
    AtomTypeIds() { for (int& id : byAtomicNum_) id = -1; }

    int get(const RDKit::Atom* atom) {
        const unsigned int anum = atom->getAtomicNum();
        if (anum != 0 && anum < MAX_ANUM) {
            int& id = byAtomicNum_[anum];
            if (id < 0) id = count_++;
            return id;
        }
        const std::string symbol = atom->getSymbol();
        for (const auto& [s, id] : bySymbol_) if (s == symbol) return id;
        bySymbol_.emplace_back(symbol, count_);
        return count_++;
    }

    int size() const { return count_; }

private:
    static constexpr unsigned int MAX_ANUM = 128;
    int byAtomicNum_[MAX_ANUM];
    std::vector<std::pair<std::string, int>> bySymbol_;
    int count_ = 0;
};

inline double getWeightedBondDegree(const RDKit::Bond* bond){
    switch (bond->getBondType()) {
        case RDKit::Bond::SINGLE:   return 1.0;
        case RDKit::Bond::DOUBLE:   return 2.0;
        case RDKit::Bond::TRIPLE:   return 3.0;
        case RDKit::Bond::AROMATIC: return 1.5;
        default:                    return 0.0;
    }
}

#endif // ATOM_TYPES_H
