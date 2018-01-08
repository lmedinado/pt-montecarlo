#pragma once

#include <random>
#include "theory.hpp"
#include "utils.hpp"


template <int Ndims, int Nfields>
class theory;

template <int Nfields>
struct phis {
    double phi[Nfields];
};

template <int Ndims, int Nfields>
class lattice {
    typedef typename std::array<int, Ndims> idxtype;
    typedef ::theory<Ndims, Nfields> theory;
    typedef ::phis<Nfields> phis;
    
    private:
        ndarray<phis, Ndims> fields;

        theory & th;

        idxtype dims;
        
        std::uniform_int_distribution<int> intdist;

    public:
        lattice(const lattice & lat)
            :fields(ndarray<phis, Ndims>(lat.fields)),
             th(lat.th), dims(lat.dims),

             intdist(0, lat.fields.size() - 1)
        {}

        lattice(const idxtype & dims_, theory & th_)
            :fields(ndarray<phis, Ndims>(dims_)),
             th(th_), dims(dims_), 
             
             intdist(0, fields.size() - 1)
        {}

        int address_of(const idxtype & idx) const {
            return fields.address_of(idx);
        }
        
        idxtype index_of(int k) const {
            return fields.index_of(k);
        }
        
        idxtype first_site() const {
            return fields.index_of(0);
        }
        idxtype last_site() const {
            return fields.index_of(get_n_sites() - 1);
        }

        idxtype next_site(const idxtype & idx) const {
            return fields.index_of(mod(fields.address_of(idx) + 1, get_n_sites()));
        }

        idxtype neighbor(const idxtype & idx, int direction, int delta, bool wrap = false) const {
            return fields.neighbor(idx, direction, delta, wrap);
        }

        const std::array<int, Ndims> & get_dims() const {
            return dims;
        }

        int get_n_sites() const {            
            return fields.size();
        }

        int max_dist() const {
            return *std::max_element(dims.begin(), dims.end());
        }

        int random_site(std::mt19937 & rng) {
            return intdist(rng);
        }

        double set_phi(int field, const idxtype & idx, double val) {
            return fields[idx].phi[field] = val;
        }

        /* phi(idx) takes the boundary conditions into consideration. */
        double phi(int field, const idxtype & idx) const {
            idxtype idxw;
            int sgn = 1;
            
            for(int i = 0; i < Ndims; i++) {
                /* the index wraps around */
                idxw[i] = mod(idx[i], dims[i]);
            
                switch(th.bc[i]) {
                    case PBC:
                        break;
                    case APBC:
                        sgn *= (idx[i] < 0? -1 : 1) * (is_odd(idx[i] / dims[i]) ? -1 : 1);
                        break;
                    case DIRC:
                        sgn *= ((idx[i] >= 0 and idx[i] < dims[i]) ? 1 : 0);
                        break;
                }
            }
            
            return sgn * fields[idxw].phi[field];
        }

        phis phi(const idxtype & idx) const {
            phis ret;

            for(int field = 0; field < th.nfields(); field++)
                ret.phi[field] = phi(field, idx);

            return ret;
        }

        double deltaphi(int field, const idxtype & idx, int dir = 0) const {
            return phi(field, neighbor(idx, dir, 1)) - phi(field, idx);
        }

        double lagrangian(const idxtype & idx) const {
            return th.lagrangian(*this, idx);
        }

        double action() const {
            double s = 0.0;

            for(auto it = fields.begin(); it < fields.end(); it++)
                s += lagrangian(fields.index_of(it - fields.begin()));

            return s;
        }

        /* measurement routines */
        double probe_onepoint_function(int field, const idxtype & idx) const {
            return phi(field, idx);
        }

        double probe_twopoint_function(int field1, int field2, const idxtype & idx1, const idxtype & idx2) const {
            return phi(field1, idx1) * phi(field2, idx2);
        }

        void fprint_config(FILE *f = stdout) const {          
           
            fprintf(f, "\t{");
            for(int field = 0; field < Nfields; field++) {
                fprintf(f, "(* begin field #%d *)", field);

                idxtype idx = first_site();
                fprintf(f, "{");
                do {
                    fprintf(f, "{");
                    for(int i = 0; i < Ndims; i++)
                        fprintf(f, "%d, ", idx[i]);

                    fprintf(f, "%f", phi(field, idx));
                    
                    fprintf(f, "}");
                    if(next_site(idx) != first_site())
                        fprintf(f, ", ");
                    
                } while((idx = next_site(idx)) != first_site());
                fprintf(f, "}");

                fprintf(f, "(* end field #%d *)", field);
                
                if(field != Nfields - 1)
                    fprintf(f, "\n,\n ");
            }
            fprintf(f, "}\n");

        }
};