// "load_mtx.hh" -- implements template functions for reading and writing .mtx files as part of the L(inear) A(rrangement) T(oolbox) library.
//
// Copyright (C) 2018 - 2019 Georgios N Printezis
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

#ifndef LOAD_MTX_HH
#define LOAD_MTX_HH

#include <fstream>
#include <algorithm>
#include <cstdlib>
#include "Graph.hh"

namespace lat {

template<typename real, typename integer>
const Graph<real, integer> load_mtx(std::string file_name) {
    std::ifstream fin(file_name);

    if (!fin.is_open()) {
        std::cerr << "unable to read file : " << file_name << '\n';

        std::exit(EXIT_FAILURE);
    }
    else {
        std::clog << "successfully opened file " << file_name << '\n';
    }

    integer rows, cols, nnz;

    while (fin.peek() == '%') {
        fin.ignore(2048, '\n');
    }

    fin >> rows >> cols >> nnz;

    std::vector<real> A(nnz);

    std::vector<integer> I(nnz);

    std::vector<integer> J(nnz);

    for (integer i = 0; i < nnz; i++) {
        fin >> I[i] >> J[i] >> A[i];
    }

    fin.close();

    std::clog << "matrix loaded to memory";

    std::transform(I.begin(), I.end(), I.begin(), [] (integer i) { 
                                                    return i -= 1; 
                                                  });
                                               
    std::transform(J.begin(), J.end(), J.begin(), [] (integer i) { 
                                                    return i -= 1; 
                                                  });

    std::vector<integer> JA(cols + 1);

    for (integer i = 0; i < nnz; i++) {
        JA[J[i]] += 1;
    }

    std::partial_sum(JA.begin(), JA.end(), JA.begin());

    return Graph<real, integer>(A, I, JA, rows, cols);
}

template<typename real, typename integer>
const std::vector<integer> load_mtx_sequence(std::string &file_name) {
    std::ifstream fin(file_name);

    if (!fin.is_open()) {
        std::cerr << "unable to read file : " << file_name << '\n';

        std::exit(EXIT_FAILURE);
    }
    else {
        std::clog << "successfully opened file " << file_name << '\n';
    }

    integer rows;

    real cost;

    while (fin.peek() == '%') {
        fin.ignore(2048, '\n');
    }

    fin >> rows >> cost;

    std::vector<integer> s(rows);

    for (integer i = 0; i < rows; i++) {
        fin >> s[i]; 
    }

    fin.close();

    std::clog << "sequence loaded to memory\n";

    std::clog << "best cost : " << std::fixed << cost << '\n';

    std::transform(s.begin(), s.end(), s.begin(), [] (integer i) { 
                                                    return i -= 1; 
                                                  });

    return s;
}

template<typename real, typename integer>
const Graph<real, integer> load_mtx_bin(std::string & file_name) {
    std::ifstream fin(file_name);

    if (!fin.is_open()) {
        std::cerr << "unable to read file : " << file_name << '\n';

        std::exit(EXIT_FAILURE);
    }
    else {
        std::clog << "successfully opened file " << file_name << '\n';
    }

    integer rows, cols, nnz;

    while (fin.peek() == '%') {
        fin.ignore(2048, '\n');
    }

    fin >> rows >> cols >> nnz;

    std::vector<integer> I(nnz);

    std::vector<integer> J(nnz);

    for (integer i = 0; i < nnz; i++) {
        fin >> I[i] >> J[i];
    }

    fin.close();

    std::clog << "matrix loaded to memory\n";

    std::transform(I.begin(), I.end(), I.begin(), [] (integer i) { 
                                                    return i -= 1; 
                                                  });

    std::transform(J.begin(), J.end(), J.begin(), [] (integer i) { 
                                                    return i -= 1; 
                                                  });

    std::vector<real> A(nnz, 1.);    

    std::vector<integer> JA(cols + 1);

    for (integer i = 0; i < nnz; i++) {
        JA[J[i]] += 1;
    }

    std::partial_sum(JA.begin(), JA.end(), JA.begin());

    return Graph<real, integer>(A, I, JA, rows, cols);
}

template<typename real, typename integer>
void write_mtx_sequence(const std::string & file_name, const real cost, const std::vector<integer> & s) {
    std::ofstream file;

    file.open(file_name);

    if (!file.is_open()) {
        std::cerr << "unable to open file for output : " << file_name << '\n';

        std::exit(EXIT_FAILURE);
    }
    else {
        std::clog << "successfully opened file " << file_name << '\n';
    }

    file << s.size() << ' ' << std::fixed << cost << '\n';

    for (auto n : s) {
        file << n + 1 << '\n'; 
    }

    file.close();

    std::clog << "successful write\n";
}

}
#endif
