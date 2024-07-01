#include "include.h"

// An arrray to store info about the cosmological evolution
double *cosmo_array;
int COSMO_ROWS;

// The logarithmic spacing between two neighboring points in T
double delta_logx_cosmo;

// A flag to determine if a cosmo_file was loaded
bool cosmo_file_loaded = false;


void read_cosmo_file(char *filename, int nrows) {
    FILE* f = fopen(filename, "r");

    if ( f == NULL ) {
        perror("Could not open the provided cosmo-file. Exit!");
        
        exit(1);
    }

    COSMO_ROWS = nrows;
    cosmo_array = (double *) malloc(COSMO_ROWS * COSMO_COLS * sizeof(double));
    //arr[i * rows + j] = ...; // logically equivalent to arr[i][j]

    for( int row = 0; row < COSMO_ROWS; row++ ) {
        fscanf(f, "%lf %lf %lf %lf %lf %lf %lf", &cosmo_array[row*COSMO_COLS+0], &cosmo_array[row*COSMO_COLS+1], &cosmo_array[row*COSMO_COLS+2], &cosmo_array[row*COSMO_COLS+3], &cosmo_array[row*COSMO_COLS+4], &cosmo_array[row*COSMO_COLS+5], &cosmo_array[row*COSMO_COLS+6]);
    }

    delta_logx_cosmo = (log(cosmo_array[(COSMO_ROWS-1)*COSMO_COLS]) - log(cosmo_array[0]))/(COSMO_ROWS - 1);

    fclose(f);

    cosmo_file_loaded = true;
}


double interp_cosmo_array(int i_col, double x) {
    if ( !cosmo_file_loaded ) {
        perror("Cannot interpolate data, since no cosmo-file has been loaded. Exit!");

        exit(1);
    }

    double r = (log(x) - log(cosmo_array[0]))/delta_logx_cosmo;
    int r_div = (int)floor(r);
    double r_mod = r - r_div;

    if ( r_div < 0 || r_div > COSMO_ROWS - 2 ) {
        perror("Index out of range in interp_cosmo_array. Exit!\n");

        exit(1);
    }

    return cosmo_array[r_div*COSMO_COLS + i_col] * (1 - r_mod) +  cosmo_array[(r_div+1)*COSMO_COLS + i_col]*r_mod;
}


double cosmo_t_T(double T) {
    if ( !cosmo_file_loaded ) {
        perror("Cannot calculate t(T) since no cosmo-file has been loaded. Exit!");

        exit(1);
    }

    double logt_min = log(cosmo_array[COSMO_COL_t]);
    double logt_max = log(cosmo_array[(COSMO_ROWS-1)*COSMO_COLS + COSMO_COL_t]);

    double eps = 1e-5;
    int N_itermax = 100;
    double logt_tmp, T_tmp;

    for( int i = 0; i < N_itermax; i++ ) {
        logt_tmp = 0.5*(logt_min + logt_max);
        T_tmp = interp_cosmo_array(COSMO_COL_T, exp(logt_tmp));
        if ( fabs((T - T_tmp)/T) < eps ) {
            return exp(logt_tmp);
        }
        if ( T_tmp > T ) {
            logt_min = logt_tmp;
        } else {
            logt_max = logt_tmp;
        }
    }

    perror("The bisection method to determine t(T) does not converge. Exit!\n");
    
    exit(1);

    return 0.;
}
