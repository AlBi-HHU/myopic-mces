// Hungarian.h: minimum cost assignment for square cost matrices.
//
// O(n^3) Hungarian algorithm with dual potentials and shortest augmenting paths
// (cf. https://cp-algorithms.com/graph/hungarian-algorithm.html).

#ifndef HUNGARIAN_H
#define HUNGARIAN_H

#include <vector>

// cost: row-major n x n matrix, cost[i * n + j] := cost of assigning row i to column j
// assignment: filled with assignment[i] = column assigned to row i
// returns the total cost of the optimal assignment
double solveAssignment(const std::vector<double>& cost, int n, std::vector<int>& assignment);

#endif
