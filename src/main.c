#include "include.h"


bool compare_rates(int err, struct parameters params, double Tmin, double Tmax, int N) {
    bool all_equal = true;

    double f1[NNUCREAC+1], f2[NNUCREAC+1];
    double r1[NNUCREAC+1], r2[NNUCREAC+1];

    for ( int i = 0; i <= NNUCREAC; i++ ) {
        f1[i] = 0.;
        f2[i] = 0.;
        r1[i] = 0.;
        r2[i] = 0.;
    }

    for ( int j = 0; j < N; ++j ) {
        double logT = log(Tmin) + log(Tmax/Tmin)*j/(N-1);

        double T = exp(logT);

        rate_pn(err, f1, r1, T, 0.71*T, params.life_neutron);
        rate_pn(err, f2, r2, T, 0.71*T, params.life_neutron);

        rate_weak(err, f1);
        rate_weak_test(err, f2);

        rate_all(err, f1, T);
        rate_all_test(err, f2, T);

        for ( int i = 1; i <= 100; ++i ) {
            if ( f1[i] != f2[i] || r1[i] != r2[i] ) {
                printf("%d %.6e %.6e %.6e %.6e\n", i, f1[i], f2[i], r1[i], r2[i]);

                all_equal = false;
            }
        }
    }

    return all_equal;
}


int main(int argc, char **argv) {
    // Default parameters
    double eta            = 6.137;
    char  *cosmo_file     = "data/sm_cosmo_file.dat";
    int    clines         = 6174;
    char  *abundance_file = "";

    // Parse the command-line arguments
    if ( argc >= 5 ) {
        abundance_file = argv[4];
    }

    if ( argc >= 4 ) {
        sscanf(argv[3], "%d", &clines);
        cosmo_file = argv[2];
    }

    if ( argc >= 2 ) {
        sscanf(argv[1], "%lf", &eta);
    }


    // Set the relevant parameters for BBN
    struct parameters params;
    params.eta0           = eta*1e-10;
    params.life_neutron   = 880.2;
    params.method         = "RK2";
    params.decay_neutrons = true;

    // Read the cosmo-file
    load_cosmo_data(cosmo_file, clines);

    // Testing
    if ( false ) for ( int err = 0; err <= 2; err++ ) compare_rates(err, params, 1e-5, 1e10, 100000);

    // Run the calculation
    double Y0[NNUC+1], Y0_high[NNUC+1],Y0_low[NNUC+1];
    // -->
    final_abundances(0, params, Y0     );
    final_abundances(1, params, Y0_high);
    final_abundances(2, params, Y0_low );

    // Print the results to the screen
    for ( int i = 1; i < 10; i++ ) {
        printf("%.6e %.6e %.6e\n", Y0[i], Y0_high[i], Y0_low[i]);
    }

    // If a filename was provided, also save the results
    // to the specified abundance-file
    if ( strcmp(abundance_file, "") != 0 ) {
        FILE* file = fopen(abundance_file, "w");

        if ( file == NULL ) {
            perror("Could not open the provided abundance-file:");

            exit(1);
        }

        for ( int i = 1; i < 10; i++ ) {
            fprintf(file, "%.6e %.6e %.6e\n", Y0[i], Y0_high[i], Y0_low[i]);
        }
    
        fclose(file);
    }

    free_cosmo_data();

    return 0;
}
