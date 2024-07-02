#include "include.h"

#define MAX_LINE_LENGTH 256


double load_eta(const char *filename) {
    char line[MAX_LINE_LENGTH];

    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        return EXIT_FAILURE;
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
