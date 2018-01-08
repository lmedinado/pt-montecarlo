#pragma once 
#include <stack>

#include "lattice.hpp"
#include "monte.hpp"
#include "utils.hpp"


struct observable {
    bool measured;
    double mean;
    double stddev;
    double stderror;
    doublev measurements;
};

/*****************************************************************************/
/**                             FREE FUNCTIONS                              **/
/*****************************************************************************/


double vector_mean(doublev v) {
    double ret = 0.0;

    for(int i = 0; i < v.size(); i++)
        ret += v[i];

    return ret / v.size();
}

double vector_stddev(doublev v, double mean) {
    double ret = 0.0;

    for(int i = 0; i < v.size(); i++)
        ret += (v[i] - mean)*(v[i] - mean);

    ret /= (v.size() - 1);

    return sqrt(ret);
}
double vector_stddev(doublev v) {
    return vector_stddev(v, vector_mean(v));
}

double vector_stderr(doublev v) {
    return vector_stddev(v) / sqrt(v.size());
}

/*****************************************************************************/
/**                                 CLASSES                                 **/
/*****************************************************************************/

template <int Ndims, int Nfields>
class stats {
    typedef typename std::array<int, Ndims> idxtype;
    typedef ::lattice<Ndims, Nfields> lattice;
    typedef ::monte<Ndims, Nfields> monte;
    
    private:
        lattice & lat;
        int n_sites;

        monte & mc;

        std::vector<lattice> configs;

        observable action, zero_crossings, connected_components;

        ndarray<observable, Ndims> onep_function[Nfields];        
        ndarray<observable, 2 * Ndims> twop_function[Nfields][Nfields];
        
        observable sp_averaged_onep_function[Nfields];

        ndarray<observable, 2> sp_averaged_twop_function[Nfields][Nfields];
        ndarray<observable, 2> wall_to_wall_function[Nfields][Nfields];

        doublev tmeasurements;
        
        /* bootstrap */
        bool bootstrap = false;
        int resamples = 0;
        
        doublev resampled;

    public:
        stats(lattice & lat_, monte & mc_, unsigned int flags = DEFAULT)
            :lat(lat_),
             mc(mc_),
             n_sites(lat_.get_n_sites())
        {
            clear_stats();
        }

        void clear_stats() {
            action.measured = false;
            action.measurements.clear();
            
            zero_crossings.measured = false;
            zero_crossings.measurements.clear();
            
            connected_components.measured = false;
            connected_components.measurements.clear();
            
            tmeasurements.clear();

            auto dims = lat.get_dims();
            auto clear_observable = [] (observable & o) {
                o.measured = false;
                o.measurements.clear();
            };
            
            for(int field1 = 0; field1 < Nfields; field1++) {
                
                sp_averaged_onep_function[field1].measured = false;
                
                onep_function[field1].resize(dims);
                onep_function[field1].map(clear_observable);

                for(int field2 = 0; field2 < Nfields; field2++) {
                    
                    twop_function[field1][field2].resize(concat_arrays<int, Ndims>(dims, dims));
                    twop_function[field1][field2].map(clear_observable);

                    sp_averaged_twop_function[field1][field2].resize({Ndims, lat.max_dist()});
                    sp_averaged_twop_function[field1][field2].map(clear_observable);
                  
                    wall_to_wall_function[field1][field2].resize({Ndims, lat.max_dist()});
                    wall_to_wall_function[field1][field2].map(clear_observable);
                    
                }
            }
        }

        void generate_configs(int nconfigs, int autocorrel) {

            char st[] = "Generating configurations";
            int window_w = 80 - std::strlen(st);
            int oldfrac = 0, newfrac = 0;

            printf("%s", st);

            configs.reserve(nconfigs);
            
            for(int i = 0; i < nconfigs; i++) {
                
                mc.metropolis(autocorrel);

                auto l(lat);
                configs.push_back(l);

                /* progress bar */
                newfrac = (int)((float) window_w  * i / nconfigs);
                if(newfrac > oldfrac) {
                    oldfrac = newfrac;
                    printf(".");
                }
            }
            printf("\nGeneration complete.\n");
            
            /* old measurements are no longer valid */
            clear_stats();
        }


        /* bootstrap estimation */
        
        void enable_bootstrap() { bootstrap = true; }
        void disable_bootstrap() { bootstrap = false; }
        
        void enable_bootstrap(int resamples_) {
            resamples = resamples_;
            enable_bootstrap();
        }

        const observable & bootstrap_estimate(observable & o, doublev measurements) {
            auto dist = std::uniform_int_distribution<int>(0, measurements.size() - 1);
            
            doublev b_mean;

            for(int i = 0; i < resamples; i++) {
                double sum = 0.0; 
                for(int j = 0; j < measurements.size(); j++) {
                    sum += measurements[dist(mc.rand())];
                }
                
                b_mean.push_back(sum / measurements.size());
            }
            
            o.mean = vector_mean(b_mean);
            o.stddev = vector_stddev(b_mean, o.mean);
            o.stderror = o.stddev; /* bootstrap stddev is already std error of the mean */

            return o;
        }


        /* generic measurement routine */

        const observable & measure_observable(observable & o, std::function<double (const lattice &)> probe, bool store_measurements = false) {

            if(!o.measured) {
                doublev & measurements(store_measurements ? o.measurements : tmeasurements);
                
                double sum = 0.0;
                
                if(measurements.size() != configs.size())
                    measurements.resize(configs.size());
                
                int i = 0;
                for(const auto & c : configs) {
                    
                    double r = probe(c);
                    
                    measurements[i++] = r;
                    sum += r;
                }
                
                if(bootstrap) {
                    bootstrap_estimate(o, measurements);
                } else {                
                    o.mean = sum / configs.size();
                    o.stddev = vector_stddev(measurements, o.mean);
                    o.stderror = o.stddev / sqrt(configs.size());
                }
                
                o.measured = true;
            }
            
            return o;
        }

        /* measurement routines for specific observables */

        const observable & measure_action(bool store_measurements = false) {

            return measure_observable(action,
                
                [&] (const lattice & c) -> double {
                    return c.action();
                },
            store_measurements);
        }

        const observable & measure_onepoint_function(int field, const idxtype & idx, bool store_measurements = false) {
            
            //printf("measuring 1pfunction at {%d, %d}", idx[0], idx[1]);
            
            return measure_observable(
                onep_function[field][idx],
                
                [&] (const lattice & c) -> double {
                    return c.probe_onepoint_function(field, idx);
                },
            store_measurements);
        }

        const observable & measure_twopoint_function(int field1, int field2, const idxtype & idx1, const idxtype & idx2, bool store_measurements = false) {
            
            /* measure ahead of time for cache coherence */
            idxtype idx = configs[0].first_site();
            do {
                measure_onepoint_function(field1, idx);
                measure_onepoint_function(field2, idx);
            } while((idx = configs[0].next_site(idx)) != configs[0].first_site());

            return measure_observable(
                twop_function[field1][field2][concat_arrays<int, Ndims>(idx1, idx2)],
                
                [&] (const lattice & c) -> double {
                    return c.probe_twopoint_function(field1, field2, idx1, idx2) - measure_onepoint_function(field1, idx1).mean*measure_onepoint_function(field2, idx2).mean;
                },
            store_measurements);

        }

        const observable & measure_onepoint_function(int field, bool store_measurements = false) {
            
            measure_observable(
                sp_averaged_onep_function[field],
                
                [&] (const lattice & c) -> double {
                    double r = 0.0;
                    
                    idxtype idx = c.first_site();
                    do {
                        r += c.probe_onepoint_function(field, idx);
                    } while((idx = c.next_site(idx)) != c.first_site()); 
                    
                    return r / n_sites;
                },
            store_measurements);
            
            // wrong but ok for now
            // needs covariance matrix
            if(!bootstrap) {
                sp_averaged_onep_function[field].stddev *= sqrt(n_sites);
                sp_averaged_onep_function[field].stderror *= sqrt(n_sites);
            }
            
            return sp_averaged_onep_function[field];
        }

        const observable & measure_twopoint_function(int field1, int field2, int delta_idx, int direction = 0, bool store_measurements = false) {

            /* measure ahead of time for cache coherence */
            idxtype idx = configs[0].first_site();
            do {
                measure_onepoint_function(field1, idx);
                measure_onepoint_function(field2, idx);
            } while((idx = configs[0].next_site(idx)) != configs[0].first_site());

            measure_observable(
                sp_averaged_twop_function[field1][field2][{direction, delta_idx}],
                
                [&] (const lattice & c) -> double {
                    double r = 0.0;
                    
                    idxtype idx = c.first_site();
                    do {
                        r += c.probe_twopoint_function(field1, field2, idx, c.neighbor(idx, direction, delta_idx, true))
                           - measure_onepoint_function(field1, idx).mean*measure_onepoint_function(field2, c.neighbor(idx, direction, delta_idx, true)).mean;
                    } while((idx = c.next_site(idx)) != c.first_site()); 
                    
                    return r / n_sites;
                },
            store_measurements);

            // wrong but ok for now
            // needs covariance matrix
            if(!bootstrap) {
                sp_averaged_twop_function[field1][field2][{direction, delta_idx}].stddev *= sqrt(n_sites);
                sp_averaged_twop_function[field1][field2][{direction, delta_idx}].stderror *= sqrt(n_sites);
            }
            //sp_averaged_twop_function[field1][field2][n_sites - {direction, delta_idx}] = sp_averaged_twop_function[field1][field2][{direction, delta_idx}];
                       
            return sp_averaged_twop_function[field1][field2][{direction, delta_idx}];
        }

        const observable & measure_wall_to_wall_function(int field1, int field2, int delta_idx, int direction = 0, bool store_measurements = false) {
            
            measure_observable(
                wall_to_wall_function[field1][field2][{direction, delta_idx}],
                
                [&] (const lattice & c) -> double {
                    double r = 0.0;
                    
                    std::vector<double> v1(c.get_dims()[direction]), v2(c.get_dims()[direction]);
                    
                    idxtype idx = c.first_site();
                    do {
                        v1[idx[direction]] += c.probe_onepoint_function(field1, idx);
                        v2[idx[direction]] += c.probe_onepoint_function(field2, idx);
                    } while((idx = c.next_site(idx)) != c.first_site()); 
                    
                    for(int i = 0; i < c.get_dims()[direction]; i++)
                        r += v1[i] * v2[mod(i + delta_idx, c.get_dims()[direction])];
                    
                    return r / n_sites;                   
                },
            store_measurements);
            
            return wall_to_wall_function[field1][field2][{direction, delta_idx}];
        }

        const observable & measure_zero_crossings(int field, bool store_measurements = false) {
            measure_observable(
                zero_crossings,
                
                [&] (const lattice & c) -> double {
                    double r = 0.0;
                    
                    idxtype idx = c.first_site();
                    do {
                        idxtype nextidx = c.next_site(idx);
                        
                        /* if sign changes, increment r */
                        if(sgn(c.phi(field, idx)) != sgn(c.phi(field, nextidx)))
                            r += 1.0;
                    } while((idx = c.next_site(idx)) != c.first_site()); 
                    
                    return r;
                },
            store_measurements);
            
            return zero_crossings;
        }

        const observable & measure_connected_components(int field, bool store_measurements = false) {
            
            measure_observable(connected_components,
            
            [&] (const lattice & c) -> double {
                double r = 0.0;
                
                ndarray<int, Ndims> colors(c.get_dims());
                colors.set_all(-1);
                
                int current_color = -1;
                
                idxtype idx = c.first_site();    
                std::stack<idxtype> to_visit;
                
                /* use a DFS to paint each lattice site with a color which
                   indexes which connected component it is a part of.      */
                do {
                    if(colors[idx] == -1) {

                        /* new connected component, new color */
                        current_color++;
                        
                        to_visit.push(idx);
                        while(!to_visit.empty()) {
                            const idxtype pidx = to_visit.top();
                            to_visit.pop();
                            
                            /* if not visited, visit */
                            if(colors[pidx] == -1)  {
                                colors[pidx] = current_color;
                                
                                /* for each neighbor of pidx with phi(field, neighbor) having the same sign as phi(field, pidx) */
                                for(int direction = 0; direction < Ndims; direction++) {
                                    for(int delta_idx : {-1, 1}) {
                                        const idxtype & nidx = c.neighbor(pidx, direction, delta_idx, true);
                                        
                                        if(sgn(c.phi(field, pidx)) == sgn(c.phi(field, nidx))) {
                                            
                                            /* if not yet visited, push onto stack */
                                            if(colors[nidx] == -1)  
                                                to_visit.push(nidx);
                                        }
                                    }
                                }
                            }
                        }
                    }
                } while((idx = c.next_site(idx)) != c.first_site()); 
                
                return (double) current_color + 1.0;
            },
            store_measurements);
            
            return connected_components;
        }

        void fprint_configs(FILE *f = stdout) {
            fprintf(f, "{ ");

            for(int i = 0; i < configs.size(); i++) {
                configs[i].fprint_config(f);

                if(i != configs.size() - 1)
                    fprintf(f, ", ");
            }

            fprintf(f, " }\n");
        }
              

};
