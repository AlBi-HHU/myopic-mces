#include <limits>
#include "LAPJV.h"

double solveAssignmentLAPJV(const std::vector<double>& cost, int n, std::vector<int>& rowsol)
{
    const double INF = std::numeric_limits<double>::infinity();
    const double EPS = 1e-10;
    auto c = [&](int i, int j) { return cost[i * n + j]; };

    rowsol.assign(n, -1);                   // rowsol[i] := column assigned to row i
    std::vector<int> colsol(n, -1);         // colsol[j] := row assigned to column j
    std::vector<double> v(n);               // column duals
    std::vector<int> matches(n, 0), freeRows(n), collist(n), pred(n);
    std::vector<double> d(n);

    // 1. column reduction (in reverse order, as in the original)
    for (int j = n - 1; j >= 0; j--) {
        int imin = 0;
        double min = c(0, j);
        for (int i = 1; i < n; i++) {
            if (c(i, j) < min) { min = c(i, j); imin = i; }
        }
        v[j] = min;
        if (++matches[imin] == 1) {
            rowsol[imin] = j;
            colsol[j] = imin;
        } else if (v[j] < v[rowsol[imin]]) {
            const int j1 = rowsol[imin];
            rowsol[imin] = j;
            colsol[j] = imin;
            colsol[j1] = -1;
        }
    }

    // 2. reduction transfer from assigned to unassigned rows
    int numFree = 0;
    for (int i = 0; i < n; i++) {
        if (matches[i] == 0) {
            freeRows[numFree++] = i;
        } else if (matches[i] == 1) {
            const int j1 = rowsol[i];
            double min = INF;
            for (int j = 0; j < n; j++) {
                if (j != j1 && c(i, j) - v[j] < min) min = c(i, j) - v[j];
            }
            v[j1] -= min;
        }
    }

    // 3. augmenting row reduction, two passes
    for (int loop = 0; loop < 2 && numFree > 0; loop++) {
        int k = 0;
        const int prevNumFree = numFree;
        numFree = 0;
        while (k < prevNumFree) {
            const int i = freeRows[k++];
            // smallest and second smallest reduced cost in row i
            double umin = c(i, 0) - v[0], usubmin = INF;
            int j1 = 0, j2 = -1;
            for (int j = 1; j < n; j++) {
                const double h = c(i, j) - v[j];
                if (h < usubmin) {
                    if (h >= umin) { usubmin = h; j2 = j; }
                    else { usubmin = umin; umin = h; j2 = j1; j1 = j; }
                }
            }
            int i0 = colsol[j1];
            // a difference below EPS counts as a tie: otherwise two rows can steal j1 from each other
            // forever when the price update is too small to change v[j1] in floating point
            const bool decrease = usubmin - umin > EPS;
            if (decrease) {
                // lower the price of j1 so that it is still the best column for i
                v[j1] -= usubmin - umin;
            } else if (i0 >= 0) {
                // tie: take the second best column instead of stealing j1
                j1 = j2;
                i0 = colsol[j2];
            }
            if (i0 >= 0) rowsol[i0] = -1;
            rowsol[i] = j1;
            colsol[j1] = i;
            if (i0 >= 0) {
                // the displaced row is processed again right away if the price decreased, otherwise later
                if (decrease) freeRows[--k] = i0;
                else freeRows[numFree++] = i0;
            }
        }
    }

    // 4. augmentation: shortest augmenting path for each remaining free row
    for (int f = 0; f < numFree; f++) {
        const int freeRow = freeRows[f];
        for (int j = 0; j < n; j++) {
            d[j] = c(freeRow, j) - v[j];
            pred[j] = freeRow;
            collist[j] = j;
        }
        // collist[0 .. low) scanned, [low .. up) to do (distance == min), [up .. n) not yet reached
        int low = 0, up = 0, last = 0, endOfPath = -1;
        double min = 0.0;
        bool unassignedFound = false;
        do {
            if (up == low) {
                // collect the columns with minimal distance
                last = low - 1;
                min = d[collist[up++]];
                for (int k = up; k < n; k++) {
                    const int j = collist[k];
                    const double h = d[j];
                    if (h <= min) {
                        if (h < min) { up = low; min = h; }
                        collist[k] = collist[up];
                        collist[up++] = j;
                    }
                }
                for (int k = low; k < up; k++) {
                    if (colsol[collist[k]] < 0) {
                        endOfPath = collist[k];
                        unassignedFound = true;
                        break;
                    }
                }
            }
            if (!unassignedFound) {
                // scan a column with minimal distance: relax via its assigned row
                const int j1 = collist[low++];
                const int i = colsol[j1];
                const double h = c(i, j1) - v[j1] - min;
                for (int k = up; k < n; k++) {
                    const int j = collist[k];
                    const double v2 = c(i, j) - v[j] - h;
                    if (v2 < d[j]) {
                        pred[j] = i;
                        if (v2 == min) {
                            if (colsol[j] < 0) {
                                endOfPath = j;
                                unassignedFound = true;
                                break;
                            }
                            collist[k] = collist[up];
                            collist[up++] = j;
                        }
                        d[j] = v2;
                    }
                }
            }
        } while (!unassignedFound);

        // update the prices of the scanned columns
        for (int k = 0; k <= last; k++) {
            const int j1 = collist[k];
            v[j1] += d[j1] - min;
        }
        // augment along the path
        int i;
        do {
            i = pred[endOfPath];
            colsol[endOfPath] = i;
            const int j1 = endOfPath;
            endOfPath = rowsol[i];
            rowsol[i] = j1;
        } while (i != freeRow);
    }

    double total = 0.0;
    for (int i = 0; i < n; i++) total += c(i, rowsol[i]);
    return total;
}
