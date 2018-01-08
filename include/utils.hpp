#pragma once

#include <array>
#include <functional>
#include <vector>

#include "pcg_random.hpp"
#include "pcg_extras.hpp"

/*****************************************************************************/
/**                             FREE FUNCTIONS                              **/
/*****************************************************************************/

typedef typename std::vector<double> doublev;

int is_odd(int a) {
    return a % 2;
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

long mod(long a, long b) {
    int r = a % b;
    return r < 0 ? r + b : r;
}


template <typename T, int Ndims>
std::array<T, 2 * Ndims> concat_arrays (const std::array<T, Ndims> &a1, const std::array<T, Ndims> &a2) {
    std::array<T, 2 * Ndims> r;
    
    for(int i = 0; i < Ndims; i++) {
        r[i]         = a1[i];
        r[i + Ndims] = a2[i];
    }
    
    return r;
}

/*****************************************************************************/
/**                                 CLASSES                                 **/
/*****************************************************************************/


class seeded_rand {
    private:
        // Seed with a real random value, if available
        pcg_extras::seed_seq_from<std::random_device> seed_source;
        
        // Using 64 bit permuted congruential generator
        pcg64 rng;

    public:
        seeded_rand():rng(seed_source) {}
        
        pcg64 & operator()() {
            return rng;
        }
};


template <typename T, int Ndims>
class ndarray {
    typedef typename std::array<int, Ndims> idxtype;
    
    private:
        std::vector<T> elems;
    
        std::array<int, Ndims> dims;
        int nsites;
    public:
        ndarray()
            :nsites(1), dims{}
        {}

        ndarray(const idxtype & dims_)
            :nsites(1), dims{}
        {
            resize(dims_);
        }

        int resize(const idxtype & dims_) {
            if(dims != dims_) {                
                nsites = 1;
                
                for(int i = 0; i < Ndims; i++) {
                    dims[i] = dims_[i];
                    nsites *= dims[i];
                }
                elems.resize(nsites);
            }
            
            return nsites;
        }

        ndarray(const ndarray<T, Ndims> & nd)
            :elems(nd.elems),
             dims(nd.dims),
             nsites(nd.nsites)
        {}

        T & operator[](const idxtype & idx) {
            return elems[address_of(idx)];
        }
        const T & operator[](const idxtype & idx) const {
            return elems[address_of(idx)];
        }

        const std::array<int, Ndims> & get_dims() const {
            return dims;
        }

        int address_of(const idxtype & idx) const {
            int block_size = 1;        
            int pos = 0;
            
            for(int i = Ndims - 1; i >= 0; i--) {
                pos += block_size * idx[i];

                block_size *= dims[i];
            }
            
            return pos;
        }
        
        /* The lattice is partitioned into two sublattices, such that no two members of each are
           ever nearest neighbors. This function returns the id of the sublattice a node belongs to. */
        static int nodegroup(const idxtype & idx) {
            int parity = 0;
            
            for(int i = 0; i < Ndims; i++) {
                parity += is_odd(idx[i]);
            }

            return is_odd(parity);
        }

        /* useful for traversing */
        idxtype index_of(int k) const {
            idxtype ret;
            int oldk=k;
            
            for(int i = Ndims - 1; i >= 0; i--) {
                ret[i] = k % dims[i];
                k /= dims[i];
            }

            return ret;
        }
        idxtype index_of(typename std::vector<T>::const_iterator it) const {
            return index_of(it - elems.begin());
        }

        idxtype neighbor(const idxtype & idx, int direction, int delta, bool wrap = false) const {
            idxtype ret(idx);
            
            ret[direction] += delta;
            
            if(wrap)
                ret[direction] = mod(ret[direction], dims[direction]);

            return ret;
        }

        typename std::vector<T>::reference front() const {
            return elems.front();
        }
        typename std::vector<T>::reference back() const {
            return elems.front();
        }

        typename std::vector<T>::const_iterator begin() const {
            return elems.begin();
        }
        typename std::vector<T>::const_iterator end() const {
            return elems.end();
        }

        typename std::vector<T>::size_type size() const {
            return elems.size();
        }

        bool empty() const {
            return elems.empty();
        }

        void map(std::function<void (T &)> operation) {
            for(T a : elems)
                operation(a);
        }

        T set_all(T a) {
            for(T & i : elems)
                i = a;

            return a;
        }
};