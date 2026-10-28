/**
 * Created by Sarah Steinbach.
 */
#include <string>
#include <vector>
#include <GraphMol/ROMol.h>
#include <cmath>
#include <map>
#include <fstream>
#include "filter0.h"
#include <GraphMol/SmilesParse/SmilesParse.h>

using namespace std;
using namespace RDKit;

map<string, int> buildAtomTypeList(RDKit::ROMol* mol);

double applyFilter0(ROMol* mol1, ROMol* mol2){
    double difference = 0.0;
    map<string,int> atomTypeList1 = buildAtomTypeList(mol1);
    map<string,int> atomTypeList2 = buildAtomTypeList(mol2);
    for (const auto& atomType : atomTypeList1) {
        //check whether both molecules have atom-type
        if (atomTypeList2.find(atomType.first)!=atomTypeList2.end()) {
            difference += abs(atomType.second-atomTypeList2[atomType.first]);
        //only molecule1 has considered atom-type -> all degrees are added to score    
        }else difference += atomTypeList1[atomType.first];
    }
    //add all degrees to tempDifference from atoms which are only in molecule2
    for (const auto& atomType : atomTypeList2) {
        if (atomTypeList1.find(atomType.first)==atomTypeList1.end()) {
            difference += atomType.second;
        }
    } 
    return difference;
}

map<string, int> buildAtomTypeList(RDKit::ROMol* mol){
    map<string,int> atomTypeList;
    for (auto atom : mol->atoms()) {
        //ignore Hydrogen
        if (atom->getSymbol() !="H") atomTypeList[atom->getSymbol()]+=1;
    }
    return atomTypeList;
}
