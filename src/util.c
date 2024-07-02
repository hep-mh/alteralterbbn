#include "include.h"


double max(double x, double y) {
    if( x < y ) return y; else return x;
}


double min(double x, double y) {
    if( x < y ) return x; else return y;
}


double factorial(int n) {
    if ( n == 0 ) return 1.;
    if ( n == 1 ) return 1.;
    if ( n == 2 ) return 2.;
    if ( n == 3 ) return 6.;
    if ( n == 4 ) return 24.;
    
    return n*factorial(n-1); 
}


char* join_strings(const char* str1, const char* str2) {
    // Calculate the length of the new string
    size_t len1 = strlen(str1);
    size_t len2 = strlen(str2);
    size_t total_len = len1 + len2 + 1; // +1 for the null terminator

    // Allocate memory for the new string
    char *result = (char*) malloc( total_len * sizeof(char) );

    // Copy the first string into the result
    strcpy(result, str1);

    // Concatenate the second string to the result
    strcat(result, str2);

    return result;
}
