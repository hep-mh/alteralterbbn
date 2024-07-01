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