/**
 * Created by Sarah Steinbach.
 */
#ifndef FILTER0_H
#define FILTER0_H

#include <string>
#include <map>
#include <cmath>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

using namespace std;
using namespace RDKit;

double applyFilter0(RDKit::ROMol* mol1, RDKit::ROMol* mol2);

#endif // FILTER0_H
