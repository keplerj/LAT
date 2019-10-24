// "parallel_full_search.hh" -- implements template function parallel_full_search as part of the L(inear) A(rrangement) T(oolbox) library.
//
// Copyright (C) 2019 Georgios N Printezis
//
// This file is part of the LAT library. This library is free
// software; you can distribute it and/or modify it under the
// the terms of the GNU General Public License as published by the
// Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.
//
// Contact me on johakepl@gmail.com.

#ifndef PARALLEL_FULL_SEARCH_HH
#define PARALLEL_FULL_SEARCH_HH

#include "Graph.hh"
#include <omp.h>
#include <cmath>
#include <vector>

namespace lat {
    
template<typename R, typename Z>
std::vector<Z> parallel_full_search(const Graph<R, Z> & G, std::vector<Z> sequence) {
    R min_cost = la(G(sequence));

    Z z =0; Z cnt = 0;

    while (z < 1) {
        z++;

        std::cout << cnt << ' '  << min_cost << '\n';

        R new_cost;

        parallel_select_best_neighbor(G, sequence, new_cost);

        if (new_cost < min_cost) {
            cnt++;

            z = 0;

            min_cost = new_cost;
        }
    }

    return sequence;
}

template<typename R, typename Z, typename ZIter, typename RIter>
void parallel_select_best_local_neighbor(const Graph<R, Z> & G, std::vector<Z> sequence, 
                                         const Z n, const Z m, const Z proc, const Z range, const Z mid, 
                                         ZIter I_begin, ZIter J_begin, RIter min_costs_begin, R min_cost) {
    const Z start = proc * range;

    const Z finish = proc != (m - 1) ? start + range : mid;

    const Z s_start = std::max(n - 1 - finish, finish);

    const Z s_finish = n - 1 - start;

    Z ii = 0, jj = 0;

    for (Z i = start; i < finish; i++) {
        for (Z j = i + 1; j < n; j++) {
            std::swap(sequence[i], sequence[j]);

            const R new_cost = la(G(sequence));

            if (new_cost < min_cost) {
                min_cost = new_cost;

                ii = i; jj = j;
            }

            std::swap(sequence[i], sequence[j]);
        }
    }

    for (Z i = s_start; i < s_finish; i++) {
        for (Z j = i + 1; j < n; j++) {
            std::swap(sequence[i], sequence[j]);

            const R new_cost = la(G(sequence));

            if (new_cost < min_cost) {
                min_cost = new_cost;

                ii = i; jj = j;
            }

            std::swap(sequence[i], sequence[j]);
        }
    }

    *(I_begin + proc) = ii; *(J_begin + proc) = jj; *(min_costs_begin + proc) = min_cost;
}

template<typename Iter>
const auto argmin(Iter start, Iter finish) {
    return std::distance(start, std::min_element(start, finish));
}

template<typename R, typename Z>
void parallel_select_best_neighbor(const Graph<R, Z> & G, std::vector<Z> & sequence, R & min_cost) {
    const Z n = numnodes(G);

    const Z m = omp_get_num_procs();

    min_cost = la(G(sequence));

    std::vector<Z> I(m, 0);

    std::vector<Z> J(m, 0);

    std::vector<R> min_costs(m, min_cost);

    const Z range = std::lround(n / (2. * m));

    const Z mid = n / 2;

#   pragma omp parallel for
    for (Z proc = 0; proc < m; proc++) {
        parallel_select_best_local_neighbor(G, sequence, n, m, proc, range, mid, 
                                            I.begin(), J.begin(), min_costs.begin(), min_cost);
    }

    const Z min_idx = argmin(min_costs.begin(), min_costs.end());

    std::swap(sequence[I[min_idx]], sequence[J[min_idx]]);

    min_cost = min_costs[min_idx];
}

}

#endif
