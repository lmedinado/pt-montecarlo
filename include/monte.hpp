#pragma once
#include <array>
#include <omp.h>

#include "utils.hpp"

extern seeded_rand my_rand;
#pragma omp threadprivate(my_rand)
seeded_rand my_rand;


template <int Ndims, int Nfields>
class monte {
    typedef typename std::array<int, Ndims> idxtype;
    typedef ::lattice<Ndims, Nfields> lattice;
    
    private:
        double temperature;

        lattice & mainlat, newlat;
        
        /* statistics */
        int naccepts[Nfields], nrejects[Nfields];
        
        /* parameters describing the type of local search to be performed */
        unsigned int search_algo;
        std::array<double, Nfields> rw_step_size;
        
        //std::mt19937 mt_rand;
        std::uniform_real_distribution<double> realdist01;
        //static thread_local seeded_rand my_rand;

    public:
        monte(lattice & lat, std::array<double, Nfields> && rw_step_size_, unsigned int flags = DEFAULT)
            :mainlat(lat),
             newlat(lat),
             
             search_algo(flags & SEARCH_MASK), rw_step_size(rw_step_size_),

             realdist01(0, 1)
             
        {
            //std::random_device seed;
            //std::array<int, std::mt19937::state_size> seed_data;
            //std::generate_n(seed_data.data(), seed_data.size(), std::ref(seed));
            //std::seed_seq seq(seed_data.begin(), seed_data.end());
            //mt_rand.seed(seq);
            
            switch (flags & START_MASK) {
                case COLD_START:
                    init(COLD_START);
                    break;
                case HOT_START:
                    init(HOT_START);
                    break;
            }
        }

        //std::mt19937 & rand() {
        //    return mt_rand;
        //}
        
        auto & rand() {
            return my_rand();
        }

        void init(int type){
            idxtype idx = mainlat.first_site();

            switch(type){
                case COLD_START:
                    printf("Cold start\n");
                    
                    do {
                        for(int field = 0; field < Nfields; field++) {
                            //mainlat.set_phi(field, idx, newlat.set_phi(field, idx, -3.0));
                            mainlat.set_phi(field, idx, newlat.set_phi(field, idx, 0.0));
                        }

                    } while((idx = mainlat.next_site(idx)) != mainlat.first_site());

                    temperature = 0.0;
                    break;

                case HOT_START:
                    printf("Hot start\n");

                    do {
                        for(int field = 0; field < Nfields; field++)
                            mainlat.set_phi(field, idx, newlat.set_phi(field, idx, (realdist01(my_rand()) - 0.5) * 2 * INFTY));

                    } while((idx = mainlat.next_site(idx)) != mainlat.first_site());

                    temperature = std::numeric_limits<double>::infinity();
                    break;
            }

            for(int i = 0; i < Nfields; i++)
                naccepts[i] = nrejects[i] = 0;
        }


        /* sets the temperature and runs thermalization steps */
        void thermalize(double temp, int steps) {
            temperature = temp;

            printf("Thermalizing at T = %f. ", temperature);

            metropolis(steps);
        }

        void local_update(double & x, double rw_step_size) {
            switch(search_algo) {
                case RANDOM_WALK:
                    x += rw_step_size * (realdist01(my_rand()) - 0.5);
                    break;
                case UNIFORM_SEARCH:
                    x = (realdist01(my_rand()) - 0.5) * 2 * INFTY;
                    break;
            }
        }

        double local_update(lattice &lat, int field, std::array<int, Ndims> idx, double rw_step_size) {
            double x;
            
            switch(search_algo) {
                case RANDOM_WALK:
                    x = lat.phi(field, idx) + rw_step_size * (realdist01(my_rand()) - 0.5);
                    break;
                case UNIFORM_SEARCH:
                    x = (realdist01(my_rand()) - 0.5) * 2 * INFTY;
                    break;
            }
            
            return lat.set_phi(field, idx, x);
        }

        auto select_site(const idxtype & idx) {
            /* grabs the subsequent lattice site */
            return mainlat.next_site(idx);
        }

        void metropolis(int steps) {

            //idxtype idx = mainlat.first_site();
            

            
            for(int i = 0; i < steps; i++) {
                for(int ngroup = 0; ngroup < 2; ngroup++) {
                    //do {
                    #pragma omp parallel for schedule(static)// schedule(dynamic,1024)
                    //for(int addr = ngroup; addr < mainlat.get_n_sites(); addr += 2) {
                    for(int addr = ngroup; addr < mainlat.get_n_sites(); addr++) {
//printf("numthreads:%d\n",omp_get_num_threads());
                        idxtype idx = mainlat.index_of(addr);
                        if(ndarray<int, Ndims>::nodegroup(idx) == ngroup)
                            continue;
                        //idxtype idx = mainlat.index_of(addr);
                        //if(ndarray<int, Ndims>::nodegroup(idx) == ngroup) {
                          //  addr++;
                          //  idx = mainlat.index_of(addr);
                       // }
                        
                        for(int field = 0; field < Nfields; field++) {
                            /* updates temp lattice */
                            local_update(newlat, field, idx, rw_step_size[field]);

                            if(exp(-delta_s(idx) / temperature) > realdist01(my_rand())) {    /* accept */
                                mainlat.set_phi(field, idx, newlat.phi(field, idx));
                                naccepts[field]++;
                            } else {    /* reject */
                                newlat.set_phi(field, idx, mainlat.phi(field, idx));
                                nrejects[field]++;
                            }
                        }
                    }
                    //} while((idx = select_site(idx)) != mainlat.first_site());
                    #pragma omp barrier
                }
            }

        }

        double delta_s(idxtype idx) {
            double r = newlat.lagrangian(idx) - mainlat.lagrangian(idx);
            
            for (int direction = 0; direction < Ndims; direction++) {
                r += newlat.lagrangian(newlat.neighbor(idx, direction, -1)) - mainlat.lagrangian(mainlat.neighbor(idx, direction, -1));
            }
            
            return r;
        }


        /* ATTN */
        void report_acceptance() {
            printf("\n\nPhi acceptance rate: %.1f%% (%d accepts, %d rejects)\n", (double) 100.0 * naccepts[0]/( naccepts[0] + nrejects[0]), naccepts[0], nrejects[0]);
            printf("Pi acceptance rate: %.1f%% (%d accepts, %d rejects)\n", (double) 100.0 * naccepts[1]/( naccepts[1] + nrejects[1]), naccepts[1], nrejects[1]);
            if(Nfields > 2)printf("Pi2 acceptance rate: %.1f%% (%d accepts, %d rejects)\n", (double) 100.0 * naccepts[2]/( naccepts[2] + nrejects[2]), naccepts[2], nrejects[2]);
            printf("\n");
        }
};
