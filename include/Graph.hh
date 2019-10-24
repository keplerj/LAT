// "Graph.hh" -- implements template class Graph as part of the L(inear) A(rrangement) T(oolbox) library.
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


#ifndef GRAPH_HH
#define GRAPH_HH

#define NDEBUG

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cassert>
#include <omp.h>

namespace lat {

template<typename real, typename integer>
class Graph final {
public:
    Graph(const std::vector<real> & _A, 
          const std::vector<integer> & _IA, 
          const std::vector<integer> & _JA,
          const integer _rows, const integer _cols) : 
    A{_A}, IA{_IA}, JA{_JA}, rows{_rows}, cols{_cols} { ; }

    Graph(const Graph<real, integer> & _G) : 
    A{_G.A}, IA{_G.IA}, JA{_G.JA}, rows{_G.rows}, cols{_G.cols} { ; }

    Graph<real, integer> & operator=(const Graph<real, integer> & _G) {
        A = _G.A; IA = _G.IA; JA = _G.JA;

        rows = _G.rows; cols = _G.cols;

        return (*this);
    }

    const Graph<real, integer> operator()(const std::vector<integer> & p) const {
        Graph<real, integer> res{*this};

        return res.perm(p);
    }

    template<typename R, typename Z>
    friend const Z nnz(const Graph<R, Z> & G);

    template<typename R, typename Z>
    friend const Z numnodes(const Graph<R, Z> & G);

    template<typename R, typename Z>
    friend const R la(const Graph<R, Z> & G);

    template<typename R, typename Z>
    friend const R stable_la(const Graph<R, Z> & G);

    void print() const;
        
    ~Graph() { ; }

private:
    std::vector<real> A;

    std::vector<integer> IA, JA;

    integer rows, cols;

    Graph<real, integer> & c_perm(const std::vector<integer> & p); 

    Graph<real, integer> & r_perm(const std::vector<integer> & p); 

    Graph<real, integer> & perm(const std::vector<integer> & p) {
        return (*this).c_perm(p).r_perm(p);
    }

    Graph<real, integer> & perm(const std::vector<integer> & rp, const std::vector<integer> & cp) {
        return (*this).c_perm(cp).r_perm(rp);
    }
};


template<typename real, typename integer>
Graph<real, integer> & Graph<real, integer>::c_perm(const std::vector<integer> & p) {
    const integer pcols = p.size();

    std::vector<integer> resJA(pcols + 1);
                
    resJA[0] = 0;

    for (integer j = 0; j < pcols; j++) {
        const integer pj = p[j];

        const integer ub = JA[pj + 1], lb = JA[pj];

        resJA[j + 1] = (ub - lb);
    }

    std::partial_sum(resJA.begin(), resJA.end(), resJA.begin());
                        
    const integer nonzeros = resJA[pcols];

    std::vector<real> resA(nonzeros);

    std::vector<integer> resIA(nonzeros);

    for (integer j = 0; j < pcols; j++) {
        const integer pj = p[j];

        const integer ub = JA[pj + 1], lb1 = JA[pj], lb2 = resJA[j];

        std::copy(A.begin() + lb1, A.begin() + ub, resA.begin() + lb2);

        std::copy(IA.begin() + lb1, IA.begin() + ub, resIA.begin() + lb2);
    }
    
    A = std::move(resA); IA = std::move(resIA); JA = std::move(resJA); cols = pcols;

    return (*this);
}

template<typename real, typename integer>
Graph<real, integer> & Graph<real, integer>::r_perm(const std::vector<integer> & p) {
    const integer prows = p.size();

    if (prows < rows) {

    std::vector<integer> h(rows, - 1);

    for (integer i = 0; i < prows; i++) {
        h[p[i]] = i;
    }
    
    const integer nonzeros = IA.size();

    std::vector<integer> resIA(nonzeros);

    std::vector<real> resA(nonzeros, - 1);

    for (integer i = 0; i < nonzeros; i++) {
        auto hIAi = h[IA[i]];

        resIA[i] = hIAi;

        if (hIAi > - 1) {
            resA[i] = A[i];
        }
    }

    std::vector<integer> resJA(cols + 1);

    resJA[0] = 0;

    for (integer j = 1; j <= cols; j++) {
        const integer ub = JA[j], lb = JA[j - 1];

        const integer elements = std::count_if(resIA.begin() + lb, resIA.begin() + ub, [] (integer i) {
                                                                                           return i > - 1;
                                                                                       });
        resJA[j] = elements;
    }

    std::partial_sum(resJA.begin(), resJA.end(), resJA.begin());

    resIA.erase(std::remove(resIA.begin(), resIA.end(), - 1), resIA.end());

    resA.erase(std::remove(resA.begin(), resA.end(), - 1), resA.end());

    A = std::move(resA); IA = std::move(resIA); JA = std::move(resJA); rows = prows;

    }
    else {

    std::vector<integer> h(rows);

    for (integer i = 0; i < rows; i++) {
        h[p[i]] = i;
    }
    
    const integer nonzeros = nnz(*this);

    std::vector<integer> resIA(nonzeros);

    for (integer i = 0; i < nonzeros; i++) {
        resIA[i] = h[IA[i]];
    }

    IA = std::move(resIA);

    }

    return (*this);
}

template<typename R, typename Z>
const Z nnz(const Graph<R, Z> & G) {
    return G.A.size();
}

template<typename R, typename Z>
const Z numnodes(const Graph<R, Z> & G) {
    return G.cols;
}

template<typename R, typename Z>
const R la(const Graph<R, Z> & G) {
    //std::vector<R> costs(nnz(G));

    const std::vector<R> & A = G.A;

    const std::vector<Z> & IA = G.IA;

    const std::vector<Z> & JA = G.JA;

    const Z cols = G.cols;

    R total_cost = .0;

    for (Z j = 0; j < cols; j++) {
        const Z ub = JA[j + 1], lb = JA[j];

        for (Z i = lb; i < ub; i++) {
            total_cost += A[i] * std::fabs(j - IA[i]);
        }
    }

    //const R total_cost = std::accumulate(costs.begin(), costs.end(), 0.);

    return total_cost;
}

template<typename R, typename Z>
const R stable_la(const Graph<R, Z> & G) {
    std::vector<R> costs(nnz(G));

    const std::vector<R> & A = G.A;

    const std::vector<Z> & IA = G.IA;

    const std::vector<Z> & JA = G.JA;

    const Z cols = G.cols;

    for (Z j = 0; j < cols; j++) {
        const Z ub = JA[j + 1], lb = JA[j];

        for (Z i = lb; i < ub; i++) {
            costs[i] = A[i] * std::fabs(j - IA[i]);
        }
    }

    std::sort(costs.begin(), costs.end());

    const R total_cost = std::accumulate(costs.begin(), costs.end(), 0.);

    return total_cost;
}

template<typename real, typename integer>
void Graph<real, integer>::print() const {
    for (integer j = 0; j < cols; j++) {
        const integer ub = JA[j + 1], lb = JA[j];

        for (integer i = lb; i < ub; i++) {
            std::cout << "(" << IA[i] << ", " << j << ")\t" << A[i] << "\n";                
        }
    }
}

}

#endif
