#include "theory.hpp"

template <int Ndims, int Nfields>
class free_field_Nd : public theory<Ndims, Nfields> {
    typedef typename std::array<int, Ndims> idxtype;
    typedef ::lattice<Ndims, Nfields> lattice;
    
    public:
        double m0 = 1.0;
        double spacing = 1.0;

        double mu2_1 = 0.035;
        double mu2_2 = 0.055;
        double gg = 0.0;
        
        //unsigned int bc[MAXDIMS] = {PBC, PBC, PBC, PBC};
        unsigned int bc[MAXDIMS] {PBC, PBC, PBC, PBC};
        
        static constexpr int ndims = Ndims;
        static constexpr int nfields = Nfields;
        
//        std::array<double, nfields> rec_rw_step = {{4.0, 4.0}};
        std::array<double, nfields> rec_rw_step {{4.0, 4.0}};
        
        std::array<double, nfields> recommended_rw_step_size() const { return rec_rw_step; }
        
        int get_ndims() const { return ndims; }
        int get_nfields() const { return nfields; }
        
        double kinetic_term(const lattice & lat, int field, const idxtype & idx) const {
            double r = 0.0;
            
            for(int direction = 0; direction < Ndims; direction++)
                r += 0.5 * m0 * pow(lat.deltaphi(field, idx, direction), 2) / (spacing * spacing);
            
            return r;
        }
        
        double lagrangian(const lattice & lat, const idxtype & idx) const {
            return spacing * (kinetic_term(lat, 0, idx) + kinetic_term(lat, 1, idx) + potential(lat, idx));
        }

        double potential(const lattice & lat, const idxtype & idx) const {
            return 0.5 * mu2_1 * pow(lat.phi(0, idx), 2) + 0.5 * mu2_2 * pow(lat.phi(1, idx), 2);
        }
        
        void print_params(FILE *f = stdout) const {
            fprintf(f, "mu2_1: %f, mu2_2: %f, g: %f", mu2_1, mu2_2, gg);
        } 
};

template <int Ndims>
class basic_ig_coupling : public theory<Ndims, 1 + Ndims> {  
    protected:
        typedef typename std::array<int, Ndims> idxtype;
        typedef ::lattice<Ndims, 1 + Ndims> lattice;
    public:
        double m0 = 1.0;
        double spacing = 1.0;
        
        //unsigned int bc[MAXDIMS] = {{PBC, PBC, PBC, PBC}};
        unsigned int bc[MAXDIMS] {PBC, PBC, PBC, PBC};
        
        static constexpr int ndims = Ndims;
        static constexpr int nfields = 1 + ndims;
        
        int get_ndims() const { return ndims; }
        int get_nfields() const { return nfields; }

        double kinetic_term(const lattice & lat, const idxtype & idx) const {
            double r = 0.0;
            
            for(int direction = 0; direction < Ndims; direction++)
                r += 0.5 * m0 * pow(lat.deltaphi(0, idx, direction), 2) / (spacing * spacing);
            
            return r;
        }
        
        double pi_kinetic_term(const lattice & lat, const idxtype & idx) const {
            double r = 0.0;
            
            for(int direction = 0; direction < Ndims; direction++)
                r += 0.5 * m0 * pow(lat.phi(1 + direction, idx), 2) / (spacing * spacing);
            
            return r;
        }
        
        double lagrangian(const lattice & lat, const idxtype & idx) const {
            return spacing * (kinetic_term(lat, idx) + pi_kinetic_term(lat, idx) + potential(lat, idx));
        }
        
        virtual double potential(const lattice & lat, const idxtype & idx) const = 0;
        
};

class free_ig_coupling : public basic_ig_coupling<1> {
    public:
        double mu2_1 = 0.001;
        //double mu2_2 = 0.25;
        double mu2_2 = 0.002;
        double gg = 0.1;
        
        /* recommended step sizes for these parameters      */
        /*                                                  */
        /* const double mu2_1 = 0.0015;                     */
        /* const double mu2_2 = 0.0025;                     */
        /* const double gg = 0.10;                          */
        /*                                                  */
        
        std::array<double, nfields> rec_rw_step = {{2.0, 0.2}};
//        std::array<double, nfields> rec_rw_step {2.0, 0.2};
        //std::array<double, nfields> rec_rw_step = {4.0, 2.0};
        
        std::array<double, nfields> recommended_rw_step_size() const { return rec_rw_step; }
        
        double potential(const lattice & lat, const idxtype & idx) const {
            return 0.5 * mu2_1 * lat.phi(0, idx) * lat.phi(0, idx) +
               + (0.5 / mu2_2) * pow(lat.deltaphi(1, idx) + gg * lat.phi(0, idx), 2);
        }
        
        void print_params(FILE *f = stdout) const {
            fprintf(f, "mu2_1: %f, mu2_2: %f, g: %f", mu2_1, mu2_2, gg);
        } 
};

class icy_1d : public basic_ig_coupling<1> {
    public:
        double mu2_1 = 0.080;
        double mu2_2 = 0.070;
        double gg = 0.1;

        
        /* recommended step sizes for these parameters      */
        /*                                                  */
        /* const double mu2_1 = 0.0015;                     */
        /* const double mu2_2 = 0.0025;                     */
        /* const double gg = 2.8;                           */
        /*                                                  */
        
        std::array<double, nfields> rec_rw_step = {{3.7, 1.0}};
        //std::array<double, nfields> rec_rw_step {3.7, 1.0};
        
        std::array<double, nfields> recommended_rw_step_size() const { return rec_rw_step; }
        
        double potential(const lattice &lat, const idxtype & idx) const {
            return 0.5 * mu2_1 * lat.phi(0, idx) * lat.phi(0, idx) +
               + (0.5 / mu2_2) * pow(lat.deltaphi(1, idx) + pow(gg * lat.phi(0, idx), 2), 2);
        }
 
        void print_params(FILE *f = stdout) const {
            fprintf(f, "mu2_1: %f, mu2_2: %f, g: %f", mu2_1, mu2_2, gg);
        } 
};

class icy_2d : public basic_ig_coupling<2> {
    public:
        double mu2_1 = 0.0015;
        double mu2_2 = 0.0015;
        double gg = 0.02;
       
        std::array<double, nfields> rec_rw_step = {{2.2, 0.2, 0.2}};
        //std::array<double, nfields> rec_rw_step {2.2, 0.2, 0.2};
        
        std::array<double, nfields> recommended_rw_step_size() const { return rec_rw_step; }
        
        double potential(const lattice & lat, const idxtype & idx) const {
            double del_mu_pi = 0.0;
            
            for(int direction = 0; direction < 2; direction++) {
                del_mu_pi += lat.deltaphi(1 + direction, idx, direction);
            }
            
            return 0.5 * mu2_1 * pow(lat.phi(0, idx), 2)
                + (0.5 / mu2_2) * pow(del_mu_pi + gg * pow(lat.phi(0, idx), 2), 2);
        }
 
        void print_params(FILE *f = stdout) const {
            fprintf(f, "mu2_2: %f mu2_2: %f, g: %f", mu2_1, mu2_2, gg);
        } 
};

class ig_double_well : public basic_ig_coupling<1> {
    public:
        double mu2_2 = 0.5;
        double lambda = 1.0;
        double aa = 1.0;
        double gg = 1.0;

        std::array<double, nfields> rec_rw_step = {{0.6, 0.8}};
        //std::array<double, nfields> rec_rw_step {0.6, 0.8};
        
        std::array<double, nfields> recommended_rw_step_size() const { return rec_rw_step; }
        
        double potential(const lattice &lat, const idxtype & idx) const {
            return lambda * pow((pow(lat.phi(0, idx), 2) - aa*aa), 2)
                + (0.5 / mu2_2) * pow(lat.deltaphi(1, idx) +  gg * lat.phi(0, idx), 2);
        }
 
        void print_params(FILE *f = stdout) const {
            fprintf(f, "aa: %f, mu2_2: %f, lambda: %f, g: %f", aa, mu2_2, lambda, gg);
        } 
};


class ig_double_well_2d : public basic_ig_coupling<2> {
    public:
        /* these are for T = 1.0 */
        double mu2_2 = 0.5;
        double lambda = 0.1;
        double vv = 3.0;
        double gg = 1.0;
        double hh = 0.0;

       
        std::array<double, nfields> rec_rw_step = {{2.3, 2.6, 2.6}};
        //std::array<double, nfields> rec_rw_step {2.3, 2.6, 2.6};
        //std::array<double, nfields> rec_rw_step = { {2.4, 4.5, 4.5} };
        
        std::array<double, nfields> recommended_rw_step_size() const { return rec_rw_step; }
        
        double potential(const lattice & lat, const idxtype & idx) const {
            double del_mu_pi = 0.0;
            
            for(int direction = 0; direction < 2; direction++) {
                del_mu_pi += lat.deltaphi(1 + direction, idx, direction);
            }
            
            return lambda * pow((pow(lat.phi(0, idx), 2) - vv*vv), 2) + hh * lat.phi(0, idx)
                + (0.5 / mu2_2) * pow(del_mu_pi + gg * lat.phi(0, idx), 2);
        }
 
        void print_params(FILE *f = stdout) const {
            fprintf(f, "vv: %f, mu2_2: %f, lambda: %f, g: %f, h: %f", vv, mu2_2, lambda, gg, hh);
        } 
};

class ig_double_well_3d : public basic_ig_coupling<3> {
    public:
        double mu2_2 = 0.5;
        double lambda = 0.1;
        double vv = 3.0;
        double gg = 0.9;
       
        std::array<double, nfields> rec_rw_step = {{2.4, 4.5, 4.5, 4.5}};
        //std::array<double, nfields> rec_rw_step {2.4, 4.5, 4.5, 4.5};
        
        std::array<double, nfields> recommended_rw_step_size() const { return rec_rw_step; }
        
        double potential(const lattice & lat, const idxtype & idx) const {
            double del_mu_pi = 0.0;
            
            for(int direction = 0; direction < 3; direction++) {
                del_mu_pi += lat.deltaphi(1 + direction, idx, direction);
            }
            
            return lambda * pow((pow(lat.phi(0, idx), 2) - vv*vv), 2)
                + (0.5 / mu2_2) * pow(del_mu_pi + gg * lat.phi(0, idx), 2);
        }
 
        void print_params(FILE *f = stdout) const {
            fprintf(f, "vv: %f, mu2_2: %f, lambda: %f, g: %f", vv, mu2_2, lambda, gg);
        } 
};


