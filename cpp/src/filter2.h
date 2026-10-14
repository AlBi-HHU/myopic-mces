/**
 * Created by Sarah Steinbach.
 */
#ifndef FILTER2_H
#define FILTER2_H

#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <GraphMol/ROMol.h>
#include "Hungarian.h"

using namespace std;
using namespace RDKit;

double applyFilter2(RDKit::ROMol* mol1, RDKit::ROMol* mol2);

#endif // FILTER2_H
