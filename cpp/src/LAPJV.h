// LAPJV.h: minimum cost assignment for square cost matrices.
//
// Shortest augmenting path algorithm of R. Jonker and A. Volgenant, "A shortest augmenting path
// algorithm for dense and sparse linear assignment problems", Computing 38, 325-340 (1987).

#ifndef LAPJV_H
#define LAPJV_H

#include <vector>

// cost: row-major n x n matrix, cost[i * n + j] := cost of assigning row i to column j
// assignment: filled with assignment[i] = column assigned to row i
// returns the total cost of the optimal assignment
double solveAssignmentLAPJV(const std::vector<double>& cost, int n, std::vector<int>& assignment);

#endif
