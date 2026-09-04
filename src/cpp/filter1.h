/**
 * Created by Sarah Steinbach.
 */
#ifndef FILTER1_H
#define FILTER1_H

#pragma once

#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <GraphMol/ROMol.h>
#include <GraphMol/Bond.h>

using namespace std;
using namespace RDKit;

double applyFilter1(RDKit::ROMol* mol1, RDKit::ROMol* mol2);

#endif // FILTER1_H
