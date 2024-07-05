#include "include.h"

// The data that is stored in the provided cosmo-file
double **cosmo_data;
int      COSMO_ROWS;

// A flag to determine if a cosmo_file was loaded
bool cosmo_data_loaded = false;


SortOrder determine_sort_order(const double *arr, int size) {
    if (size < 2) return ASCENDING;

    return arr[0] < arr[1] ? ASCENDING : DESCENDING;
}


int find_index(const double *arr, int size, double x) {
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


int extract_line_number(const char *filename) {
    FILE *file = fopen(filename, "r");
    if ( file == NULL ) {
        fprintf(stderr, "%s Could not determine number of lines in file '%s': %s.\n", ERROR, filename, strerror(errno));
        
        exit(EXIT_FAILURE);
    }

    long line_count = 0;

    int c;
    while ( (c = fgetc(file)) != EOF ) {
        if (c == '\n') {
            line_count++;
        }
    }

    if ( ferror(file) ) {
        fprintf(stderr, "%s Could not determine number of lines in file '%s': %s.\n", ERROR, filename, strerror(errno));

        exit(EXIT_FAILURE);
    }

    fclose(file);

    return line_count;
}


void load_cosmo_data(const char *filename) {
    FILE *f = fopen(filename, "r");

    if ( f == NULL ) {
        fprintf(stderr, "%s Could not open file '%s': %s.\n", ERROR, filename, strerror(errno));
        
        exit(EXIT_FAILURE);
    }

    // Read the number of lines in the cosmo-file
    COSMO_ROWS = extract_line_number(filename);

    // Allocate space for the cosmo_data array
    cosmo_data = malloc( COSMO_COLS * sizeof(double*) );
    for ( int col = 0; col < COSMO_COLS; col++ ) {
        cosmo_data[col] = malloc( COSMO_ROWS * sizeof(double) );
    }

    // Read the data
    for ( int row = 0; row < COSMO_ROWS; row++ ) {
        int read_columns = fscanf(f, "%lf %lf %lf %lf %lf %lf", &cosmo_data[0][row], &cosmo_data[1][row], &cosmo_data[2][row], &cosmo_data[3][row], &cosmo_data[4][row], &cosmo_data[5][row]);
    
        if ( read_columns != COSMO_COLS ) {
            fprintf(stderr, "%s Wrong number of columns/row in file '%s'.\n", ERROR, filename);

            exit(80);
        }

        // Adapt the units (hbar in GeV s)
        cosmo_data[0][row] *= 1;         // s     --> s
        cosmo_data[1][row] *= 1e-3;      // MeV   --> GeV
        cosmo_data[2][row] *= 1e-6/hbar; // MeV^2 --> GeV/s
        cosmo_data[3][row] *= 1e-3;      // MeV   --> GeV
        cosmo_data[4][row] *= 1e-3/hbar; // MeV   --> 1/s
        cosmo_data[5][row] *= 1e-9;      // MeV^3 ---> GeV^3

        // Consume the rest of the line
        fscanf(f, "%*[^\n]\n");
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
    if ( !cosmo_data_loaded ) {
        fprintf(stderr, "%s Cannot perform interpolation: data not initialized.\n", ERROR);

        exit(83);
    }

    int ix = find_index(cosmo_data[xc], COSMO_ROWS, x);

    if ( ix == -1 ) {
        fprintf(stderr, "%s Cannot perform interpolation: index out of range for x = %.3e.\n", ERROR, x);

        exit(84);
    }

    bool is_between = ( cosmo_data[xc][ix] <= x && x <= cosmo_data[xc][ix+1] ) ||
                      ( cosmo_data[xc][ix] >= x && x >= cosmo_data[xc][ix+1] );
    // Check if ix has the desired properties
    if ( ix > COSMO_ROWS - 2 || !is_between ) {
        fprintf(stderr, "%s Could not find i such that x is between x[i] and x[i+1]\n", ERROR);

        exit(85);
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


double time(double T) {
    return interp_cosmo_data(T, COSMO_COL_T, COSMO_COL_t);
}


double time_from_bisection(double T) {
    if ( !cosmo_data_loaded ) {
        fprintf(stderr, "%s Cannot calculate t(T): data not initialized.\n", ERROR);

        exit(85);
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

    fprintf(stderr, "%s Cannot calculate t(T): bisection method does not converge.\n", ERROR);
    
    exit(86);

    return 0.;
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
