#include <algorithm>
#include <limits>
#include "Hungarian.h"

double solveAssignment(const std::vector<double>& cost, int n, std::vector<int>& assignment)
{
    const double INF = std::numeric_limits<double>::infinity();
    // 1-based indices, index 0 is a virtual column used as the root of each augmenting path
    std::vector<double> u(n + 1, 0.0), v(n + 1, 0.0), minv(n + 1);
    std::vector<int> p(n + 1, 0), way(n + 1, 0);   // p[j] := row matched to column j
    std::vector<char> used(n + 1);

    // Initial feasible potentials (row and column reduction) and greedy matching on tight entries
    // (reduced cost 0). The cost matrices of filter 2 contain many zeros (atoms with identical
    // neighborhoods), so most rows are matched here and only few augmenting paths remain.
    std::vector<int> rowMatched(n + 1, 0);
    for (int i = 1; i <= n; i++) {
        const double* row = &cost[(i - 1) * n];
        u[i] = *std::min_element(row, row + n);
    }
    for (int j = 1; j <= n; j++) {
        double m = INF;
        for (int i = 1; i <= n; i++) m = std::min(m, cost[(i - 1) * n + (j - 1)] - u[i]);
        v[j] = m;
    }
    for (int i = 1; i <= n; i++) {
        const double* row = &cost[(i - 1) * n];
        for (int j = 1; j <= n; j++) {
            if (p[j] == 0 && row[j - 1] - u[i] - v[j] == 0.0) {
                p[j] = i;
                rowMatched[i] = 1;
                break;
            }
        }
    }

    for (int i = 1; i <= n; i++) {
        if (rowMatched[i]) continue;
        p[0] = i;
        int j0 = 0;
        std::fill(minv.begin(), minv.end(), INF);
        std::fill(used.begin(), used.end(), 0);
        do {
            used[j0] = 1;
            const int i0 = p[j0];
            const double* row = &cost[(i0 - 1) * n];
            double delta = INF;
            int j1 = 0;
            for (int j = 1; j <= n; j++) {
                if (used[j]) continue;
                const double cur = row[j - 1] - u[i0] - v[j];
                if (cur < minv[j]) { minv[j] = cur; way[j] = j0; }
                if (minv[j] < delta) { delta = minv[j]; j1 = j; }
            }
            for (int j = 0; j <= n; j++) {
                if (used[j]) { u[p[j]] += delta; v[j] -= delta; }
                else minv[j] -= delta;
            }
            j0 = j1;
        } while (p[j0] != 0);
        // augment along the found path
        do {
            const int j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        } while (j0 != 0);
    }

    assignment.assign(n, -1);
    double total = 0.0;
    for (int j = 1; j <= n; j++) {
        assignment[p[j] - 1] = j - 1;
        total += cost[(p[j] - 1) * n + (j - 1)];
    }
    return total;
}
