// "full_search.hh" -- implements template function full_search for use in the L(inear) A(rrangement) T(oolbox) library.
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
 
#ifndef FULL_SEARCH_HH
#define FULL_SEARCH_HH

#include "Graph.hh"

namespace lat {
    
template<typename R, typename Z>
std::vector<Z> full_search(const Graph<R, Z> & G, std::vector<Z> sequence) {
    R min_cost = la(G(sequence));

    Z z =0, cnt = 0;

    while (z < 1) {
        z++;

        std::cout << cnt << ' '  << min_cost << '\n';

        select_best_neighbor(G, sequence);

        R new_cost = la(G(sequence));

        if (new_cost < min_cost) {
            cnt++;

            z = 0;

            min_cost = new_cost;
        }
    }

    return sequence;
}

template<typename R, typename Z>
void select_best_neighbor(const Graph<R, Z> & G, std::vector<Z> & sequence) {
    Z n = numnodes(G);

    R min_cost = la(G(sequence));

    Z ii =0, jj = 0;

    for (auto i = 0; i < n - 1; i++) {
        for (auto j = i + 1; j < n; j++) {
            std::swap(sequence[i], sequence[j]);

            R new_cost = la(G(sequence));

            if (new_cost < min_cost) {
                min_cost = new_cost;

                ii = i; jj = j;
            }
            
            std::swap(sequence[i], sequence[j]);
        }
    }

    std::swap(sequence[ii], sequence[jj]);
}

}

#endif

