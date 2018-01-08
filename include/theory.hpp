#pragma once

#include <cstdio>
#include <array>
#include "lattice.hpp"
#include "constants.hpp"

template <int Ndims, int Nfields>
class lattice;

template <int Ndims, int Nfields>
class theory {
    protected: 
        typedef typename std::array<int, Ndims> idxtype;
        typedef ::lattice<Ndims, Nfields> lattice;
    public:

        theory() {}

        virtual ~theory() {}
        
        unsigned int bc[MAXDIMS];
        
        /* theory dependent constant parameters */
        virtual std::array<double, Nfields> recommended_rw_step_size() const = 0;
        virtual int get_ndims() const = 0;
        virtual int get_nfields() const = 0;
        
        virtual double potential(const lattice & lat, const idxtype & idx) const = 0;
        virtual double lagrangian(const lattice & lat, const idxtype & idx) const = 0;
        virtual void print_params(FILE *f = stdout) const = 0;
};