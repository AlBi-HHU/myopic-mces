/**
 * Created by Sarah Steinbach.
 */
#include <string>
#include <vector>
#include <GraphMol/ROMol.h>
#include <GraphMol/Bond.h>
#include "filter1.h"

using namespace std;
using namespace RDKit;

// Vorab-Deklarationen der Hilfsfunktionen
double getWeightedDegree(ROMol* mol, Atom* atom);
bool compareByDegree(ROMol* mol, Atom* atom1, Atom* atom2);
map<string, vector<Atom*>> builtPartitions(ROMol* mol);

// Die Hauptfunktion, die wir später in Python aufrufen
double applyFilter1(ROMol* mol1, ROMol* mol2) {
    // Partitionen bauen
    map<string, vector<Atom*>> partition1 = builtPartitions(mol1);
    map<string, vector<Atom*>> partition2 = builtPartitions(mol2);

    double tempDifference = 0;

    for (const auto& atomType : partition1) {
        // Prüfen, ob beide Moleküle den Atomtyp besitzen
        if (partition2.find(atomType.first) != partition2.end()) {
            vector<Atom*> atomsSortedByDegree1 = partition1[atomType.first];
            vector<Atom*> atomsSortedByDegree2 = partition2[atomType.first];

            int minSize = min(atomsSortedByDegree1.size(), atomsSortedByDegree2.size());
            for (int i = 0; i < minSize; i++) {
                double degree1 = getWeightedDegree(mol1, atomsSortedByDegree1[i]);
                double degree2 = getWeightedDegree(mol2, atomsSortedByDegree2[i]);
                tempDifference += abs(degree1 - degree2);
            }

            // Wenn ein Molekül mehr Atome dieses Typs hat, Differenzen aufaddieren
            if (atomsSortedByDegree1.size() < atomsSortedByDegree2.size()) {
                for (size_t i = minSize; i < atomsSortedByDegree2.size(); i++) {
                    tempDifference += getWeightedDegree(mol2, atomsSortedByDegree2[i]);
                }
            } else {
                for (size_t i = minSize; i < atomsSortedByDegree1.size(); i++) {
                    tempDifference += getWeightedDegree(mol1, atomsSortedByDegree1[i]);
                }
            }
        }
        // Nur Molekül 1 hat diesen Atomtyp
        else {
            for (auto atom : partition1[atomType.first]) {
                tempDifference += getWeightedDegree(mol1, atom);
            }
        }
    }

    // Atome aufaddieren, die nur in Molekül 2 vorkommen
    for (const auto& atomType : partition2) {
        if (partition1.find(atomType.first) == partition1.end()) {
            for (auto atom : partition2[atomType.first]) {
                tempDifference += getWeightedDegree(mol2, atom);
            }
        }
    }

    return tempDifference / 2.0;
}

// Erstellt die Partitionen für das Molekül
map<string, vector<Atom*>> builtPartitions(ROMol* mol) {
    map<string, vector<Atom*>> partition;

    for (auto atom : mol->atoms()) {
        partition[atom->getSymbol()].push_back(atom);
    }

    // Sortiert Atome absteigend nach ihrem gewichteten Grad
    for (auto& pair : partition) {
        // Wir übergeben 'mol' per Value an das Lambda, statt 'this' zu nutzen
        sort(pair.second.begin(), pair.second.end(),
            [mol](Atom* atom1, Atom* atom2) {
                return compareByDegree(mol, atom1, atom2);
            });
    }
    return partition;
}

bool compareByDegree(ROMol* mol, Atom* atom1, Atom* atom2) {
    return getWeightedDegree(mol, atom1) > getWeightedDegree(mol, atom2);
}

double getWeightedDegree(ROMol* mol, Atom* atom) {
    double degree = 0.0;

    for (auto& bond : mol->atomBonds(atom)) {
        Bond::BondType bondType = bond->getBondType();

        switch (bondType) {
            case Bond::SINGLE:   degree += 1.0; break;
            case Bond::DOUBLE:   degree += 2.0; break;
            case Bond::TRIPLE:   degree += 3.0; break;
            case Bond::AROMATIC: degree += 1.5; break;
            default:             degree += 0.0; break;
        }
    }
    return degree;
}
