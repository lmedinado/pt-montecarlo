#define _CRT_SECURE_NO_WARNINGS

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iterator>
#include <limits>
#include <numeric>
#include <random>
#include "include/constants.hpp"
#include "include/utils.hpp"
#include "include/theory.hpp"
#include "include/models.hpp"
#include "include/monte.hpp"
#include "include/stats.hpp"


const bool PRINT_ERRORS = false;
const bool PRINT_CONFIG_TO_FILE = false;

const int NMEAS = 15000;
const int NSWEEPS = 2;


int main(int argc, char *argv[]){

    ig_double_well_2d th;
    //ig_double_well_2d th;
    //icy_2d th;
    //icy_1d th;
    //free_ig_coupling th;
    
    //free_field_Nd<2> th;

    
    //{64,64}, {256}
    lattice<th.ndims, th.nfields> lat({64,64}, th);
    monte<th.ndims, th.nfields> mc(lat, th.recommended_rw_step_size(), HOT_START);
    stats<th.ndims, th.nfields> st(lat, mc);

    double temperature = 1.0;

    int nmeas = NMEAS;
    int nsweeps = NSWEEPS;

    mc.thermalize(temperature, 1024);

    mc.report_acceptance();

    printf("\n");

    st.generate_configs(20, 200);


    printf("\n\naction mean: %f\n", st.measure_action().mean / lat.get_n_sites() );
    printf("action stderr: %f\n", st.measure_action().stderror / sqrt(lat.get_n_sites()));

    printf("S/T = %f\n", st.measure_action().mean / lat.get_n_sites() / temperature);


    printf("\n1p function(phi1) mean: %f\n", st.measure_onepoint_function(0).mean);

    printf("1p function(phi1) stderr: %f\n", st.measure_onepoint_function(0).stderror);
    
    for(int field1 = 0; field1 < 1; field1++) {
        for(int field2 = 0; field2 < 1; field2++) {
            
            for(int direction = 0; direction < 1; direction++) { //th.ndims
                printf("\n\t(* <%s,%s>_%d  ", (field1? "pi" : "phi"), (field2? "pi" : "phi"), direction);
                th.print_params();
                printf(" *)");
                
                printf("\n\t{");
                for(int i = 0; i < 32; i++) {
                    if(PRINT_ERRORS)
                        printf("{%f, %f}", st.measure_twopoint_function(field1, field2, i, direction).mean, st.measure_twopoint_function(field1, field2, i, direction).stderror);
                        //{std::array<int, 1> idx1{0}, idx2{i}; printf("{%f, %f}, ", st.measure_twopoint_function(field1, field2, idx1, idx2).mean, st.measure_twopoint_function(field1, field2, idx1, idx2).stderror);}
                    else
                        printf("%f", st.measure_twopoint_function(field1, field2, i, direction).mean);
                    //printf("%f, ", st.measure_wall_to_wall_function(field1, field2, i).mean);
                    if(i != 31)
                        printf(", ");
                }
                printf("}\n");

                printf("\n");
            }
        }
    }

    printf("\n");
    printf("{g, zero crossings}: {%f, %f},\n", th.gg, st.measure_zero_crossings(0).mean);
    printf("zero crossings stddev: %f", st.measure_zero_crossings(0).stddev);

    printf("\n");
    printf("{g, # of domains}: {%f, %f},\n", th.gg, st.measure_connected_components(0).mean);
    printf("# of domains stddev: %f", st.measure_connected_components(0).stddev);

    
    
    // FILE *f = fopen("out.txt", "w");
    // lat.fprint_config(f);
    // fclose(f);


    if(PRINT_CONFIG_TO_FILE or true) {
        FILE *f;

        f = fopen("out.txt", "w");

        printf("\nSaving configurations to file out.txt\n");
        st.fprint_configs(f);

        fclose(f);
    }

    return 0;
}
