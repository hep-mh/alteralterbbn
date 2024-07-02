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
            if ( f1[i] != f2[i] ) {
                printf("f: %d %.6e\n", i, fabs(f1[i]-f2[i])/f1[i]);

                all_equal = false;
            }

            if ( r1[i] != r2[i] ) {
                printf("r: %d %.6e\n", i, fabs(r1[i]-r2[i])/r1[i]);

                all_equal = false;
            }
        }
    }

    return all_equal;
}


int main(int argc, char **argv) {
    // Default parameters
    double  eta          = 6.137;
    char   *io_directory = "io/sm";

    // Parse the command-line arguments
    if ( argc >= 3 ) {
        io_directory = argv[2];
    }

    if ( argc >= 2 ) {
        sscanf(argv[1], "%lf", &eta);
    }

    // -->
    char *cosmo_file_name     = join_strings(io_directory, "/cosmo_file.dat");
    char *abundance_file_name = join_strings(io_directory, "/abundance_file.dat");
    char *param_file_name     = join_strings(io_directory, "/param_file.dat");

    // Set the relevant parameters for BBN
    struct parameters params;
    params.eta0           = eta*1e-10;
    params.life_neutron   = 879.4;
    params.method         = "RK2";
    params.decay_neutrons = true;

    // Load the cosmological data
    load_cosmo_data(cosmo_file_name);

    // Testing
    if ( false ) for ( int err = 0; err <= 2; err++ ) compare_rates(err, params, 1e-5, 1e10, 100000);

    // Run the calculation
    double Y0m[NNUC+1], Y0h[NNUC+1],Y0l[NNUC+1];
    // -->
    final_abundances(0, params, Y0m);
    final_abundances(1, params, Y0h);
    final_abundances(2, params, Y0l);

    // Free the cosmological data
    free_cosmo_data();

    // Print the results to the screen
    for ( int i = 1; i < 10; i++ ) {
        printf("%.6e %.6e %.6e\n", Y0m[i], Y0h[i], Y0l[i]);
    }


    // WRITE ABUNDANCE_FILE ///////////////////////////////////////////////////
    FILE *abundance_file = fopen(abundance_file_name, "w");

    if ( abundance_file == NULL ) {
        fprintf(stderr, "ERROR: Could not open the file '%s': %s\n", abundance_file_name, strerror(errno));

        exit(EXIT_FAILURE);
    }

    for ( int i = 1; i < 10; i++ ) {
        fprintf(abundance_file, "%.6e %.6e %.6e\n", Y0m[i], Y0h[i], Y0l[i]);
    }

    fclose(abundance_file);


    // WRITE PARAM_FILE ///////////////////////////////////////////////////////
    FILE *param_file = fopen(param_file_name, "w");

    if ( param_file == NULL ) {
        fprintf(stderr, "ERROR: Could not open the file '%s': %s\n", param_file_name, strerror(errno));

        exit(EXIT_FAILURE);
    }

    fprintf(param_file, "eta=%lfe-10\n", eta);

    fclose(param_file);


    return 0;
}
