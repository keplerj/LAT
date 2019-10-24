 // "successive_augmentation.hh" -- implements template function successive_augmentation as part of the L(inear) A(rrangement) T(oolbox) library.
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
 
#ifndef SUCCESSIVE_AUGMENTATION
#define SUCCESSIVE_AUGMENTATION

#include "Graph.hh"

namespace lat {

template<typename R, typename Z>
std::vector<Z> successive_augmentation(const Graph<R, Z> & G, const std::vector<Z> & initial_sequence) {
    Z n = numnodes(G);

    std::vector<Z> sequence;

    sequence.reserve(n);

    Z cend, mid1, mid2, pos;

    R min_cost;

    if (n % 2) {
       cend = 1;

       mid1 = mid2 = n / 2 + 1;

       sequence.push_back(initial_sequence[mid1]);
    }
    else {
        cend = 2;

        mid1 = n / 2 - 1;

        mid2 = n / 2;

        sequence.push_back(initial_sequence[mid1]);

        sequence.push_back(initial_sequence[mid2]);
    }

    for (int i = 0; i < mid1; i++) {
        cend++;

        //std::cout << cend << '\n';

        sequence.push_back(initial_sequence[mid1 - 1 - i]);

        min_cost = la(G(sequence));

        pos = cend;

        for (int j = cend - 1; j > 0; j--) {
            std::swap(sequence[j], sequence[j - 1]);

            R cost = la(G(sequence));

            if (cost < min_cost) {
                min_cost = cost;

                pos = j;
            }
        }

        std::rotate(sequence.begin(), sequence.begin() + 1, sequence.begin() + pos);

        cend++;

        //std::cout << cend << '\n';

        sequence.push_back(initial_sequence[mid2 + 1 + i]);

        min_cost = la(G(sequence));

        pos = cend;

        for (int j = cend - 1; j > 0; j--) {
            std::swap(sequence[j], sequence[j - 1]);

            R cost = la(G(sequence));

            if (cost < min_cost) {
                min_cost = cost;

                pos = j;
            }
        }

        std::rotate(sequence.begin(), sequence.begin() + 1, sequence.begin() + pos);
    }

    return sequence;
}

}

#endif
