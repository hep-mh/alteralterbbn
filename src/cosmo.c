#include "include.h"

// An arrray to store info about the cosmological evolution
double **cosmo_data;
int COSMO_ROWS;

// A flag to determine if a cosmo_file was loaded
bool cosmo_data_loaded = false;


SortOrder determine_sort_order(double *arr, int size) {
    if (size < 2) return ASCENDING;

    return arr[0] < arr[1] ? ASCENDING : DESCENDING;
}


int find_index(double *arr, int size, double x) {
    if (size < 2) return -1; // Cannot determine position for arrays smaller than 2

    SortOrder order = determine_sort_order(arr, size);
    
    int left = 0;
    int right = size - 1;

    while (left <= right) {
        int mid = left + (right - left) / 2;

        if ((order == ASCENDING && arr[mid] <= x && (mid == size - 1 || arr[mid + 1] > x)) ||
            (order == DESCENDING && arr[mid] >= x && (mid == size - 1 || arr[mid + 1] < x))) {
            return mid;
        }
        
        if ((order == ASCENDING && arr[mid] < x) || (order == DESCENDING && arr[mid] > x)) {
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }

    return -1;
}


void load_cosmo_data(char *filename, int nrows) {
    FILE* f = fopen(filename, "r");

    if ( f == NULL ) {
        perror("Could not open the provided cosmo-file:");
        
        exit(1);
    }

    COSMO_ROWS = nrows;

    cosmo_data = malloc(COSMO_COLS * sizeof(double*));
    for (int col = 0; col < COSMO_COLS; col++ ) {
        cosmo_data[col] = malloc(COSMO_ROWS * sizeof(double));
    }

    for( int row = 0; row < COSMO_ROWS; row++ ) {
        fscanf(f, "%lf %lf %lf %lf %lf %lf %lf", &cosmo_data[0][row], &cosmo_data[1][row], &cosmo_data[2][row], &cosmo_data[3][row], &cosmo_data[4][row], &cosmo_data[5][row], &cosmo_data[6][row]);
    }

    fclose(f);

    cosmo_data_loaded = true;
}


void free_cosmo_data() {
    for (int col = 0; col < COSMO_COLS; col++ ) {
        free(cosmo_data[col]);
    }

    free(cosmo_data);

    cosmo_data_loaded = false;
}


double interp_cosmo_data(double x, int xc, int yc) {
    int ix = find_index(cosmo_data[xc], COSMO_ROWS, x);

    if ( ix == -1 || ix == COSMO_ROWS - 1 ) {
        fprintf(stderr, "Cannot interpolate data: index out of range\n");

        exit(1);
    }

    double x1 = cosmo_data[xc][ix], x2 = cosmo_data[xc][ix+1];
    double y1 = cosmo_data[yc][ix], y2 = cosmo_data[yc][ix+1];

    int sign = ( y1 < 0 && y2 < 0 ) ? -1 : 1;
    // -->
    y1 *= sign;
    y2 *= sign;

    // -->
    double m = log(y2/y1)/log(x2/x1);

    return sign * y2 * pow(x/x2, m);
}


double temperature(double t) {
    return interp_cosmo_data(t, COSMO_COL_t, COSMO_COL_T);
}


double neutrino_temperature(double t) {
    return interp_cosmo_data(t, COSMO_COL_t, COSMO_COL_Tnu);
}


double dTdt(double t) {
    return interp_cosmo_data(t, COSMO_COL_t, COSMO_COL_dTdt);
}


double nb_eta_final_ratio(double t) {
    return interp_cosmo_data(t, COSMO_COL_t, COSMO_COL_nBEtaFinalRatio);
}


double cosmo_t_T(double T) {
    if ( !cosmo_data_loaded ) {
        fprintf(stderr, "ERROR: Cannot calculate t(T) since no cosmo-file has been loaded\n");

        exit(1);
    }

    double logt_min = log(cosmo_data[COSMO_COL_t][0]);
    double logt_max = log(cosmo_data[COSMO_COL_t][COSMO_ROWS-1]);

    double eps = 1e-5;
    int N_itermax = 100;
    double logt_tmp, T_tmp;

    for( int i = 0; i < N_itermax; i++ ) {
        logt_tmp = 0.5*(logt_min + logt_max);
        T_tmp = temperature(exp(logt_tmp));
        if ( fabs((T - T_tmp)/T) < eps ) {
            return exp(logt_tmp);
        }
        if ( T_tmp > T ) {
            logt_min = logt_tmp;
        } else {
            logt_max = logt_tmp;
        }
    }

    fprintf(stderr, "ERROR: The bisection method to determine t(T) does not converge\n");
    
    exit(1);

    return 0.;
}
