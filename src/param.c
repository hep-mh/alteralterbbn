#include "include.h"

#define MAX_LINE_LENGTH 256


double load_eta(const char *filename) {
    char line[MAX_LINE_LENGTH];

    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "%s Could not open file '%s': %s\n", ERROR, filename, strerror(errno));
        
        exit(EXIT_FAILURE);
    }

    double eta = -1;
    while ( fgets(line, sizeof(line), file) ) {
        if ( sscanf(line, "eta=%lf", &eta) == 1 ) {
            break;
        }
    }

    fclose(file);
    
    return eta;
}
