/**
 * Created by Sarah Steinbach.
 */
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <GraphMol/ROMol.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include "Hungarian.h"
#include "filter2.h"

using namespace std;
using namespace RDKit;

double getWeightedBondDegree(ROMol* mol, Bond* bond);
bool compareBondsByDegree(ROMol* mol, Bond* bond1, Bond* bond2);
double getWeightedAtomDegree(ROMol* mol, Atom* atom);
bool compareAtomsByDegree(Atom* atom1, Atom* atom2);
map<string, vector<Atom*>> builtPartitions2(ROMol* mol);
double minEdgeMappingCost(ROMol* mol1, ROMol* mol2, Atom* atom1, Atom* atom2);

unordered_map<string, vector<Bond*>> incidentBondsByElement(Atom* atom, vector<Bond*> incidentBonds);
vector<Bond*> getIncidentBonds(ROMol* mol, Atom* atom);



double applyFilter2(ROMol* mol1, ROMol* mol2){
    //build partitions
    map<string,vector<Atom*>> partition1 = builtPartitions2(mol1); // contains name of partition (key is name of atom) with associated list of degrees
    map<string,vector<Atom*>> partition2 = builtPartitions2(mol2);  

    //calculate lower bound
    double tempDifference = 0.0;

    for (const auto& atomType : partition1){
        //check whether both molecules have atom-type
        //atomType: werte paar, atomType.first: string, 
        
        if (partition2.find(atomType.first)==partition2.end()){ // if molecule 2 does not have that atom type at all 
            for (Atom* atom : atomType.second) { //                 -> just add the weight of the atoms
                double cost;
                atom->getProp("weightedDegree",cost);
                tempDifference += cost; 
            }

        }
        else { // both have at least one of that type
            vector<Atom*> atomsMol1 = atomType.second;
            vector<Atom*> atomsMol2 = partition2[atomType.first];

            const int maxSize = max(atomsMol1.size(),atomsMol2.size());
            // (re)init costmatrix
            vector<vector<double>> costMatrix = vector<vector<double>>(maxSize, vector<double>(maxSize));// costMatrix[i][j] := minimal cost for mapping atomsMol1.get(i) onto atomsMol2.get(j)
            for (int i = 0; i < atomsMol1.size(); i++){
                for(int j = 0; j < atomsMol2.size();j++){
                    costMatrix[i][j] = minEdgeMappingCost(mol1, mol2, atomsMol1[i],atomsMol2[j]);
                }
               
                
            }
            
            // Implicit addition of pseudo atoms which have no incident bonds:
            if (atomsMol1.size() < atomsMol2.size()){
                for (int i = atomsMol1.size(); i<maxSize; i++){
                    for (int j = 0; j<maxSize; j++){
                        double cost;
                        atomsMol2[j]->getProp("weightedDegree",cost);
                        costMatrix[i][j] = cost;
                    }
                    
                }
            } else if (atomsMol1.size() > atomsMol2.size()){
                for (int i = 0; i < maxSize; i++){
                    for (int j = atomsMol2.size(); j<maxSize; j++){
                        double cost;
                        atomsMol1[i]->getProp("weightedDegree",cost);
                        costMatrix[i][j] = cost;
                    }
                }
            }

            

            // Minimum Bipartite Matching:
            HungarianAlgorithm hungarianAlgorithm;
            vector<int> bestAtomMapping;
            hungarianAlgorithm.Solve(costMatrix, bestAtomMapping);

            for(int i = 0; i < maxSize; i++){
                tempDifference += costMatrix[i][bestAtomMapping[i]];
            } 
        }
    }

    //deals with case where one molecule has atom-types which the other one hasn't
    for (const auto atomType : partition2){
        if (partition1.find(atomType.first)==partition1.end()){
            double tempCost;
            for (Atom* atom : atomType.second) {
                atom->getProp("weightedDegree",tempCost);
                tempDifference += tempCost;
            }
        }
        
    }
    return tempDifference/2;
    
}

//builds partition for molecule
//each atom-type gets a vector in which the occurring atoms are listed in descending order of degree
map<string,vector<Atom*>> builtPartitions2(ROMol* mol){
        map<string,vector<Atom*>> partition;

        // if atom isn't part of partition yet -> new pair with type-name and vector with atom is created
        //otherwise atom gets added to vector
        map<string, int> currentIndexMap; 
        for (auto atom : mol->atoms()) {
            partition[atom->getSymbol()].push_back(atom);
            atom->setProp("weightedDegree",0.0);
            atom->setProp("weightedDegree",getWeightedAtomDegree(mol,atom));

            
            // Set unique atom IDs
            string atomType = atom->getSymbol();  // get symbols  (C, H, etc.)
            if (currentIndexMap.count(atomType) > 0) { //symbol has already entry
                currentIndexMap[atomType]++; // -> increment
            } else { // no entry yet
                currentIndexMap[atomType] = 1; // -> initialize
            }
            atom->setProp("ID", atomType + to_string(currentIndexMap[atomType]));
        }

        //sorts atoms in decending order according to their degree
        for (auto &pair : partition){
            string atomTyp = pair.first;     // Atomtyp (key)
            sort(pair.second.begin(), pair.second.end(),
                [mol](Atom* atom1, Atom* atom2) {
                return compareAtomsByDegree(atom1, atom2);
                });
        }
        return partition;
}

bool compareAtomsByDegree(Atom* atom1, Atom* atom2){
    double degree1;
    double degree2;
    atom1->getProp("weightedDegree",degree1);
    atom2->getProp("weightedDegree",degree2);
    return degree1 > degree2;
}
double getWeightedAtomDegree(ROMol*  mol, Atom *atom){
        double degree;
        for (auto &bond : mol -> atomBonds(atom)) {
            atom->getProp("weightedDegree", degree);
            degree += getWeightedBondDegree(mol, bond);
            atom->setProp("weightedDegree", degree);
        }
        return degree;
    }

bool compareBondsByDegree(ROMol*  mol, Bond* bond1, Bond* bond2){
    return getWeightedBondDegree(mol,bond1) > getWeightedBondDegree(mol,bond2);
}

double getWeightedBondDegree(ROMol*  mol, Bond *bond){
    double degree = 0.0;
    Bond::BondType bondType = bond->getBondType();
    switch (bondType) {
        case Bond::SINGLE:
            degree += 1.0; 
            break;
        case Bond::DOUBLE:
            degree += 2.0;  
            break;
        case Bond::TRIPLE:
            degree += 3.0;  
            break;
        case Bond::AROMATIC:
            degree += 1.5;  
            break;
        default:
            degree += 0.0;  
            break;
    }
    return degree;
}

double minEdgeMappingCost(ROMol* mol1, ROMol* mol2, Atom* atom1, Atom* atom2){
    // 1. Partition the incident edges:
    
    string atomID1 = "NotSET", atomID2 = "NotSET";
    if (atom1->hasProp("ID")) {
        atom1->getProp("ID", atomID1);
    }
    if (atom2->hasProp("ID")) {
        atom2->getProp("ID", atomID2);
    }
    unordered_map<string, vector<Bond*>> incidentAtomsWithBonds1 = incidentBondsByElement(atom1, getIncidentBonds(mol1,atom1));
    unordered_map<string, vector<Bond*>> incidentAtomsWithBonds2 = incidentBondsByElement(atom2, getIncidentBonds(mol2,atom2));
    double mappingCost = 0.0;

    // 2. Compute the minimum mapping cost:
    // 2.1. Iterate over all atom types in element2IncBondsMol1 and
    // compute for each atom type the best mapping with minimum cost:
    for(const auto& pair : incidentAtomsWithBonds1){
        vector<Bond*> incBonds1 = pair.second;
        vector<Bond*> incBonds2 = {};
        // check if adjacent atom also exists in incidentAtomsWithBonds2
        auto it = incidentAtomsWithBonds2.find(pair.first);
        if (it != incidentAtomsWithBonds2.end() && !it->second.empty()) incBonds2 = it->second; // then use found bond vector instead of empty vector
        
        sort(incBonds1.begin(), incBonds1.end(),
                [mol1](Bond* bond1, Bond* bond2) {
                return compareBondsByDegree(mol1, bond1, bond2);
                });
        sort(incBonds2.begin(), incBonds2.end(),
                [mol2](Bond* bond1, Bond* bond2) {
                return compareBondsByDegree(mol2, bond1, bond2);
                });
        int minSize = min(incBonds1.size(),incBonds2.size());
        for(int i = 0; i< minSize; i++){ // compare minSize first (and sorted) bonds (that occur in both bond vectors)
            Bond* bond1 = incBonds1[i];
            Bond* bond2 = incBonds2[i];
            double newVal = abs(getWeightedBondDegree(mol1, bond1)-getWeightedBondDegree(mol2,bond2));
            mappingCost += newVal;
        }

        //deals with case where one molecule has more atoms of one atom type than the other
        if(incBonds1.size() < incBonds2.size()){
            for(int i = minSize; i < incBonds2.size(); i++){
                double newVal = getWeightedBondDegree(mol2, incBonds2[i]);
                mappingCost += newVal; 

                string atomID2_1 = "NotSET", atomID2_2 = "NotSET";
                incBonds2[i]->getBeginAtom()->getProp("ID", atomID2_1);
                incBonds2[i]->getEndAtom()->getProp("ID", atomID2_2);
            }
        } else if (incBonds1.size() > incBonds2.size()){
            for(int i = minSize; i < incBonds1.size(); i++){
                double newVal = getWeightedBondDegree(mol1, incBonds1[i]);
                mappingCost += newVal; 
                
                string atomID1_1 = "NotSET", atomID1_2 = "NotSET";
                incBonds1[i]->getBeginAtom()->getProp("ID", atomID1_1);
                incBonds1[i]->getEndAtom()->getProp("ID", atomID1_2);
            }
        }
    }
    
    //deals with case where one molecule has atom-types which the other one hasn't
    for (auto& pair : incidentAtomsWithBonds2){
        auto it = incidentAtomsWithBonds1.find(pair.first);
        if (it == incidentAtomsWithBonds1.end() || it->second.empty()) {
            vector <Bond*> incBonds = pair.second;
            for (auto& incBond : incBonds) {
                double newVal = getWeightedBondDegree(mol2, incBond);
                mappingCost += newVal;

                string atomID2_1 = "NotSET", atomID2_2 = "NotSET";
                incBond->getBeginAtom()->getProp("ID", atomID2_1);
                incBond->getEndAtom()->getProp("ID", atomID2_2);
            } 
        }
    }
    return mappingCost;
}

vector<Bond*> getIncidentBonds(ROMol* mol, Atom* atom){
    vector<Bond*> incidentBonds;
    auto bondIterator  = mol->getAtomBonds(atom);
    
    for (auto bond = bondIterator.first; bond != bondIterator.second; bond++){
        auto bondDescriptor = *bond;
        unsigned int sourceIndex = bondDescriptor.m_source;
        unsigned int targetIndex = bondDescriptor.m_target;
        
        Bond* bondObjekt = mol->getBondBetweenAtoms(sourceIndex,targetIndex);
        incidentBonds.push_back(bondObjekt);
    }
    return incidentBonds;
}

//returns map with adjacent atoms and their bond
unordered_map<string, vector<Bond*>> incidentBondsByElement(Atom* atom, vector<Bond*> incidentBonds){
    unordered_map<string, vector<Bond*>> element2BondList;
    for (auto bond: incidentBonds){
        element2BondList[bond->getOtherAtom(atom)->getSymbol()].push_back(bond);
    }
    return element2BondList;
}



