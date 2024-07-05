#include "include.h"


int linearize(double T, double reacparam[NNUCREAC+1][10], double f[NNUCREAC+1], double r[NNUCREAC+1], int loop, int inc, int ip, double dt, double Y0[NNUC+1], double Y[NNUC+1], double dY_dt[NNUC+1], double rhob) {    
    int i, j, g, h, k, l, n, ind, rn1, rn2, rn3, rn4, rn5, rn6;
    double cn1, cn2, cn3, cn4, cn5, cn6;
    double yY[NNUC+1];
    // -->
    cn1 = cn2 = cn3 = cn4 = cn5 = cn6 = 0.;

    int fail;

    int type[NNUCREAC+1], n1[NNUCREAC+1], n2[NNUCREAC+1], n3[NNUCREAC+1], n4[NNUCREAC+1], n5[NNUCREAC+1], n6[NNUCREAC+1];
    double rev[NNUCREAC+1], q9[NNUCREAC+1];
    double a[NNUC+1][NNUC+1], b[NNUC+1], yx[NNUC+1];
    double x[NNUC+1], a0[NNUC+1][NNUC+1];
    
    double cx, sum, xdy, t;
    int icnvm, nord, test;

    for ( i = 1; i <= NNUCREAC; i++ ) {
        type[i] = (int)reacparam[i][1];
        n1[i]   = (int)reacparam[i][2];
        n2[i]   = (int)reacparam[i][3];
        n3[i]   = (int)reacparam[i][4];
        n4[i]   = (int)reacparam[i][5];
        n5[i]   = (int)reacparam[i][6];
        n6[i]   = (int)reacparam[i][7];
        rev[i]  = reacparam[i][8];
        q9[i]   = reacparam[i][9];
    }
    
    memset(a, 0., sizeof(double) * (NNUC+1) * (NNUC+1));

    for ( n = 1; n <= NNUCREAC; n++ ) {
        ind = type[n];
        i = n1[n];
        j = n2[n];
        g = n3[n];
        h = n4[n];
        k = n5[n];
        l = n6[n];
        if ( i <= NNUC && l <= NNUC ) {
            rn6 = ind%10;
            rn5 = (ind%100-rn6)/10;
            rn4 = (ind%1000-10*rn5-rn6)/100;
            rn3 = (ind%10000-100*rn4-10*rn5-rn6)/1000;
            rn2 = (ind%100000-1000*rn3-100*rn4-10*rn5-rn6)/10000;
            rn1 = (ind-10000*rn2-1000*rn3-100*rn4-10*rn5-rn6)/100000;
            
            if ( ind != 100001 ) r[n]=rev[n]*pow(rhob,rn4+rn5+rn6-1.)*pow(0.987e10*pow(T,1.5),rn1+rn2+rn3-rn4-rn5-rn6)*exp(-q9[n]/T)*f[n];
            f[n] = pow(rhob,rn1+rn2+rn3-1.)*f[n];
            cn1 = (rn1*pow(Y[i],rn1-1.)*pow(Y[j],rn2)*pow(Y[g],rn3))/((rn1+rn2+rn3)*factorial(rn1)*factorial(rn2)*factorial(rn3))*f[n]*dt;
            if ( rn2== 0 ) cn2 = 0.; else cn2 = (rn2*pow(Y[j],rn2-1.)*pow(Y[i],rn1)*pow(Y[g],rn3))/((rn1+rn2+rn3)*factorial(rn1)*factorial(rn2)*factorial(rn3))*f[n]*dt;
            if ( rn3== 0 ) cn3 = 0.; else cn3 = (rn3*pow(Y[g],rn1-1.)*pow(Y[j],rn2)*pow(Y[i],rn1))/((rn1+rn2+rn3)*factorial(rn1)*factorial(rn2)*factorial(rn3))*f[n]*dt;
            if ( rn4== 0 ) cn4 = 0.; else cn4 = (rn4*pow(Y[h],rn4-1.)*pow(Y[k],rn5)*pow(Y[l],rn6))/((rn4+rn5+rn6)*factorial(rn4)*factorial(rn5)*factorial(rn6))*r[n]*dt;
            if ( rn5== 0 ) cn5 = 0.; else cn5 = (rn5*pow(Y[k],rn5-1.)*pow(Y[h],rn4)*pow(Y[l],rn6))/((rn4+rn5+rn6)*factorial(rn4)*factorial(rn5)*factorial(rn6))*r[n]*dt;
            cn6=(rn6*pow(Y[l],rn6-1.)*pow(Y[k],rn5)*pow(Y[h],rn4))/((rn4+rn5+rn6)*factorial(rn4)*factorial(rn5)*factorial(rn6))*r[n]*dt;

            // Invert the indexes
            i = NNUC+1-i;
            j = NNUC+1-j;
            g = NNUC+1-g;
            h = NNUC+1-h;
            k = NNUC+1-k;
            l = NNUC+1-l;

            // Fill i (n1) nuclide column
            a[i][i] += rn1*cn1;
            if ( j <= NNUC ) a[j][i] += rn2*cn1;
            if ( g <= NNUC ) a[g][i] += rn3*cn1;
            if ( h <= NNUC ) a[h][i] -= rn4*cn1;
            if ( k <= NNUC ) a[k][i] -= rn5*cn1;
            a[l][i] -= rn6*cn1;
            
            // Fill j (n2) nuclide column
            if ( j <= NNUC ) {
                a[i][j] += rn1*cn2;
                a[j][j] += rn2*cn2;
                if ( g <= NNUC ) a[g][j] += rn3*cn2;
                if ( h <= NNUC ) a[h][j] -= rn4*cn2;
                if ( k <= NNUC ) a[k][j] -= rn5*cn2;
                a[l][j] -= rn6*cn2;
            }
            
            // Fill g (n3) nuclide column
            if ( g <= NNUC ){
                a[i][g] += rn1*cn3;
                if ( j <= NNUC ) a[j][g] += rn2*cn3;
                a[g][g] += rn3*cn3;
                if ( h <= NNUC ) a[h][g] -= rn4*cn3;
                if ( k <= NNUC ) a[k][g] -= rn5*cn3;
                a[l][g] -= rn6*cn3;
            }

            // Fill h (n4) nuclide column
            if ( h <= NNUC ) {
                a[i][h] -= rn1*cn4;
                if ( j <= NNUC ) a[j][h] -= rn2*cn4;
                if ( g <= NNUC ) a[g][h] -= rn3*cn4;
                a[h][h] += rn4*cn4;
                if ( k <= NNUC ) a[k][h] += rn5*cn4;
                a[l][h] += rn6*cn4;
            }

            // Fill k (n5) nuclide column
            if ( k <= NNUC ) {
                a[i][k] -= rn1*cn5;
                if ( j <= NNUC ) a[j][k] -= rn2*cn5;
                if ( g <= NNUC ) a[g][k] -= rn3*cn5;
                if ( h <= NNUC ) a[h][k] += rn4*cn5;
                a[k][k] += rn5*cn5;
                a[l][k] += rn6*cn5;
            }

            // Fill l (n6) nuclide column
            a[i][l] -= rn1*cn6;
            if ( j <= NNUC ) a[j][l] -= rn2*cn6;
            if ( g <= NNUC ) a[g][l] -= rn3*cn6;
            if ( h <= NNUC ) a[h][l] += rn4*cn6;
            if ( k <= NNUC ) a[k][l] += rn5*cn6;
            a[l][l] += rn6*cn6;
        }
    }
   
    // Finish the A matrix
    for ( i = 1; i <= NNUC; i++ ) {   
        // Add identity matrix
        if ( (a[i][i] += 1.) == 0. )  {
            // If zeros at pivot points, terminate matrix evaluation
            return i;
        }; 

        // Initial abundances 
        b[NNUC+1-i]=Y0[i];
    }

    if ( loop == 1 ) icnvm = ip; else icnvm = 0;
    
    nord = 0;
    fail = 0;
    // Set RH and solution vectors to initial values
    for ( i = 1; i <= NNUC; i++) {
        x[i] = b[i];
        yx[i] = 0.;
    }
    // Save the matrix
    if ( icnvm == inc ) memcpy(a0, a, sizeof(a));

    // Triangularize the matrix
    for ( i = 1; i <= NNUC; i++ ) {
        for ( j = i+1; j <= NNUC; j++ ) {
            if ( a[j][i] != 0. ){
                cx = a[j][i]/a[i][i];
                for ( k = i+1; k <= NNUC; k++ ) a[j][k] -= cx*a[i][k];
                a[j][i] = cx;
                x[j] -= cx*x[i];
            }
        }
    }

    // Back substitution
    do {   
        x[NNUC] /= a[NNUC][NNUC];
        yx[NNUC] += x[NNUC];
        
        for( i = NNUC-1; i >= 1; i-- ) {
            sum = 0.;
            for ( j = i+1; j <= NNUC; j++ ) sum += a[i][j]*x[j];

            x[i] = (x[i]-sum)/a[i][i];
            yx[i] += x[i];
        }

        test = 1;
    
        if ( icnvm == inc ) {
            for ( i = 1; i <= NNUC; i++ ) {
                if ( yx[i] != 0. ) {
                    xdy = fabs(x[i]/yx[i]);
                    
                    if ( xdy > 2.e-4 ) {
                        if ( nord < 1 ) {
                            nord++;
                            
                            for ( j = 1; j <= NNUC; j++ ) {
                                t = 0.;
                                for ( k = 1; k <= NNUC; k++ ) t += a0[j][k]*yx[k];
                                x[j] = b[j]-t;
                            }

                            for ( j = 2; j <= NNUC; j++ ) for ( k = j+1; k <= NNUC; k++ ) x[k] -= a[k][j]*x[j];
                            break;
                        } else return -1;
                    } else test=0;
                } else test=0;
            }
        } else test=0;
    } while(test);

    // Calculate the derivatives of of the abundances
    for ( i = 1; i <= NNUC; i++ ) {
        yY[i] = yx[NNUC+1-i];

        dY_dt[i] = (yY[i]-Y0[i])/dt;
    }

    return fail;
}


int run_nucleosynthesis(int err, Parameters params, double ratioH[NNUC+1]) {
    int i;

    double dY_dt0[NNUC+1], dY_dt[NNUC+1], Y0[NNUC+1], Y[NNUC+1];
    for( i = 0 ; i <= NNUC; i++ ) {
        ratioH[i] = 0.;

        dY_dt0[i] = 0.;
        dY_dt[i]  = 0.;
        Y0[i]     = 0.;
        Y[i]      = 0.;
    }

    double f[NNUCREAC+1],r[NNUCREAC+1];
    for( i = 0; i <= NNUCREAC; i++) {
        f[i] = 0.;
        r[i] = 0.;
    }

    double dtl, dtmin;

    // Nuclides: 1=n, 2=p, 3=H2, 4=H3, 5=He3, 6=He4, 7=Li6, 8=Li7, 9=Be7, 10=Li8, 11=B8, 12=Be9, 13=B10, 14=B11, 15=C11, 16=B12, 17=C12, 18=N12, 19=C13, 20=N13, 21=C14, 22=N14, 23=O14, 24=N15, 25=O15, 26=O16
    // Atomic number A
    double Am[NNUC+1] = {0., 1., 1., 2., 3., 3., 4., 6., 7., 7., 8., 8., 9., 10., 11., 11., 12., 12., 12., 13., 13., 14., 14., 14., 15., 15., 16.};

    double reacparam[NNUCREAC+1][10] = {
        // type: #n1#n2#n3#n4#n5#n6
        // n1: incoming nuclide number
        // n2: incoming light nuclide number
        // n3: incoming lightest nuclide number
        // n4: outgoing lightest nuclide number
        // n5: outgoing light nuclide number
        // n6: outgoing nuclide number
        // rev: reverse reaction coefficient
        // q: energy release in reaction in 10**9 Kelvin (K = ev/k_B)
        {0,0,0,0,0,0,0,0,0.,0.},                   // none
        {1,100001,1,0,0,0,0,2,0.,0.},              // n -> p
        {2,100001,4,0,0,0,0,5,0.,0.},              // H3 -> e- + v + He3
        {3,100002,10,0,0,0,0,6,0.,0.},             // Li8 -> e- + v + 2 He4
        {4,100001,16,0,0,0,0,17,0.,0.},            // B12 -> e- + v + C12
        {5,100001,21,0,0,0,0,22,0.,0.},            // C14 -> e- + v + N14
        {6,100002,11,0,0,0,0,6,0.,0.},             // B8 -> e+ + v + 2 He4
        {7,100001,15,0,0,0,0,14,0.,0.},            // C11 -> e+ + v + B11
        {8,100001,18,0,0,0,0,17,0.,0.},            // N12 -> e+ + v + C12
        {9,100001,20,0,0,0,0,19,0.,0.},            // N13 -> e+ + v + C13
        {10,100001,23,0,0,0,0,22,0.,0.},           // O14 -> e+ + v + N14
        {11,100001,25,0,0,0,0,24,0.,0.},           // O15 -> e+ + v + N15
        {12,110001,2,1,0,0,0,3,0.477,25.815},      // H + n -> g + H2
        {13,110001,3,1,0,0,0,4,1.65,72.612},       // H2 + n -> g + H3
        {14,110001,5,1,0,0,0,6,2.63,238.794},      // He3 + n -> g + He4
        {15,110001,7,1,0,0,0,8,1.20,84.132},       // Li6 + n -> g + Li7
        {16,110011,5,1,0,0,2,4,1.001,8.863},       // He3 + n -> p + H3
        {17,110011,9,1,0,0,2,8,1.001,19.080},      // Be7 + n -> p + Li7
        {18,110011,7,1,0,0,4,6,1.068,55.503},      // Li6 + n -> He4 + H3
        {19,110002,9,1,0,0,0,6,4.68,220.382},      // Be7 + n -> He4 + He4
        {20,110001,3,2,0,0,0,5,1.65,63.749},       // H2 + p -> g + He3
        {21,110001,4,2,0,0,0,6,2.63,229.931},      // H3 + p -> g + He4
        {22,110001,7,2,0,0,0,9,1.20,65.053},       // Li6 + p -> g + Be7
        {23,110011,7,2,0,0,5,6,1.067,46.640},      // Li6 + p -> He4 + He3
        {24,110002,8,2,0,0,0,6,4.68,201.302},      // Li7 + p -> He4 + He4
        {25,110001,6,3,0,0,0,7,1.55,17.109},       // H2 + He4 -> g + Li6
        {26,110001,6,4,0,0,0,8,1.13,28.629},       // H3 + He4 -> g + Li7
        {27,110001,6,5,0,0,0,9,1.13,18.412},       // He3 + He4 -> g + Be7
        {28,200011,3,0,0,0,1,5,1.73,37.934},       // 2 H2 -> n + He3
        {29,200011,3,0,0,0,2,4,1.73,46.798},       // 2 H2 -> p + H3
        {30,110011,4,3,0,0,1,6,5.51,204.116},      // H3 + H2 -> n + He4
        {31,110011,5,3,0,0,2,6,5.51,212.979},      // He3 + H2 -> p + He4
        {32,200021,5,0,0,0,2,6,3.35,149.229},      // 2 He3 -> 2 p + He4
        {33,110012,8,3,0,0,1,6,9.81,175.487},      // Li7 + H2 -> n + He4 + He4
        {34,110012,9,3,0,0,2,6,9.83,194.566},      // Be7 + H2 -> p + He4 + He4
        {35,110001,5,4,0,0,0,7,2.47,183.290},      // He3 + H3 -> g + Li6
        {36,110011,7,3,0,0,1,9,2.52,39.237},       // Li6 + H2 -> n + Be7
        {37,110011,7,3,0,0,2,8,2.52,58.317},       // Li6 + H2 -> p + Li7
        {38,110011,5,4,0,0,3,6,1.59,166.181},      // He3 + H3 -> H2 + He4
        {39,200021,4,0,0,0,1,6,3.34,131.503},      // 2 H3 -> 2n + He4
        {40,110111,5,4,0,1,2,6,3.34,140.366},      // He3 + H3 -> n + p + He4
        {41,110011,8,4,0,0,1,12,3.55,121.136},     // Li7 + H3 -> n + Be9
        {42,110011,9,4,0,0,2,12,3.55,140.215},     // Be7 + H3 -> p + Be9
        {43,110011,8,5,0,0,2,12,3.55,129.999},     // Li7 + He3 -> p + Be9
        {44,110001,8,1,0,0,0,10,1.33,23.589},      // Li7 + n -> g + Li8
        {45,110001,13,1,0,0,0,14,3.07,132.920},    // B10 + n -> g + B11
        {46,110001,14,1,0,0,0,16,2.37,39.111},     // B11 + n -> g + B12
        {47,110011,15,1,0,0,2,14,1.001,32.086},    // C11 + n -> p + B11
        {48,110011,13,1,0,0,6,8,0.755,32.371},     // B10 + n -> He4 + Li7
        {49,110001,9,2,0,0,0,11,1.32,1.595},       // Be7 + p -> g + B8
        {50,110001,12,2,0,0,0,13,0.986,76.424},    // Be9 + p -> g + B10
        {51,110001,13,2,0,0,0,15,3.07,100.834},    // B10 + p -> g + C11
        {52,110001,14,2,0,0,0,17,7.10,185.173},    // B11 + p -> g + C12
        {53,110001,15,2,0,0,0,18,2.37,6.979},      // C11 + p -> g + N12
        {54,110011,16,2,0,0,1,17,3.00,146.061},    // B12 + p -> n + C12
        {55,110011,12,2,0,0,6,7,0.618,24.663},     // Be9 + p -> He4 + Li6
        {56,110011,13,2,0,0,6,9,0.754,13.291},     // B10 + p -> He4 + Be7
        {57,110011,16,2,0,0,6,12,0.291,79.903},    // B12 + p -> He4 + Be9
        {58,110001,7,6,0,0,0,13,1.60,51.761},      // Li6 + He4 -> g + B10
        {59,110001,8,6,0,0,0,14,4.07,100.549},     // Li7 + He4 -> g + B11
        {60,110001,9,6,0,0,0,15,4.07,87.543},      // Be7 + He4 -> g + C11
        {61,110011,11,6,0,0,2,15,3.07,85.948},     // B8 + He4 -> p + C11
        {62,110011,10,6,0,0,1,14,3.07,76.960},     // Li8 + He4 -> n + B11
        {63,110011,12,6,0,0,1,17,10.28,66.158},    // Be9 + He4 -> n + C12
        {64,110011,12,3,0,0,1,13,2.06,50.609},     // Be9 + H2 -> n + B10
        {65,110011,13,3,0,0,2,14,6.42,107.105},    // B10 + H2 -> p + B11
        {66,110011,14,3,0,0,1,17,14.85,159.357},   // B11 + H2 -> n + C12
        {67,210001,6,1,0,0,0,12,0.600,18.262},     // 2 He4 + n -> g + Be9
        {68,300001,6,0,0,0,0,17,2.06,84.420},      // 3 He4 -> g + C12
        {69,110012,10,2,0,0,1,6,3.54,177.713},     // Li8 + p -> n + He4 + He4
        {70,110012,11,1,0,0,2,6,3.55,218.787},     // B8 + n -> p + He4 + He4
        {71,110012,12,2,0,0,3,6,0.796,7.554},      // Be9 + p -> H2 + He4 + He4
        {72,110003,14,2,0,0,0,6,3.45,100.753},     // B11 + p -> 2 He4 + He4
        {73,110003,15,1,0,0,0,6,3.46,132.838},     // C11 + n -> 2 He4 + He4
        {74,110001,17,1,0,0,0,19,0.898,57.400},    // C12 + n -> g + C13
        {75,110001,19,1,0,0,0,21,3.62,94.884},     // C13 + n -> g + C14
        {76,110001,22,1,0,0,0,24,2.74,125.715},    // N14 + n -> g + N15
        {77,110011,20,1,0,0,2,19,1.001,34.846},    // N13 + n -> p + C13
        {78,110011,22,1,0,0,2,21,3.00,7.263},      // N14 + n -> p + C14
        {79,110011,25,1,0,0,2,24,1.001,41.037},    // O15 + n -> p + N15
        {80,110011,25,1,0,0,6,17,0.707,98.659},    // O15 + n -> He4 + C12
        {81,110001,17,2,0,0,0,20,0.896,22.554},    // C12 + p -> g + N13
        {82,110001,19,2,0,0,0,22,1.21,87.621},     // C13 + p -> g + N14
        {83,110001,21,2,0,0,0,24,0.912,118.452},   // C14 + p -> g + N15
        {84,110001,20,2,0,0,0,23,3.62,53.705},     // N13 + p -> g + O14
        {85,110001,22,2,0,0,0,25,2.73,84.678},     // N14 + p -> g + O15
        {86,110011,24,2,0,0,0,26,3.67,140.733},    // N15 + p -> g + O16
        {87,110011,24,2,0,0,6,17,0.706,57.622},	   // N15 + p -> He4 + C12
        {88,110001,17,6,0,0,0,26,5.20,83.111},     // C12 + He4 -> g + O16
        {89,110011,13,6,0,0,2,19,9.35,47.134},     // B10 + He4 -> p + C13
        {90,110011,14,6,0,0,2,21,11.03,9.098},     // B11 + He4 -> p + C14
        {91,110011,15,6,0,0,2,22,3.68,33.921},     // C11 + He4 -> p + N14
        {92,110011,18,6,0,0,2,25,4.25,111.620},    // N12 + He4 -> p + O15
        {93,110011,20,6,0,0,2,26,5.80,60.557},     // N13 + He4 -> p + O16
        {94,110011,13,6,0,0,1,20,9.34,12.288},     // B10 + He4 -> n + N13
        {95,110011,14,6,0,0,1,22,3.67,1.835},      // B11 + He4 -> n + N14
        {96,110011,16,6,0,0,1,24,4.25,88.439},     // B12 + He4 -> n + N15
        {97,110011,19,6,0,0,1,26,5.79,25.711},     // C13 + He4 -> n + O16
        {98,110011,14,3,0,0,2,16,4.96,13.296},     // B11 + H2 -> p + B12
        {99,110011,17,3,0,0,2,19,1.88,31.585},     // C12 + H2 -> p + C13
        {100,110011,19,3,0,0,2,21,7.58,69.069}     // C13 + H2 -> p + C14
    
    };

    int inc   = 50; // CHANGED: 100
    int ltime = 0;
    int is    = 1;
    int ip    = inc;
    int it    = 0;

    // Values corresponding to 'failsafe = 2'
    double cy  = 0.1;     // Limiting value of dY/dt
    double ct  = 0.005;   // Limiting value of dT/dt
    double dt0 = 1.e-4;   // CHANGED: 1e-10
    int nitmax = 50000;   // CHANGED: 100
     
    int fail = 0;

    // Different values require also changing nitmax
    double T9i = 100.; // Initial temperature, default: 27, ~ 8.6 MeV
    double T9f = 0.01; // Final temperature                 ~ 0.9 keV
    // -->
    // Initialize the temperatures for the loop
    double T9   = T9i;
    double Tnu  = T9;
    double TGeV = T9/GeV_to_GK;
    // -->
    double t  = time(TGeV); // s
    double dt = dt0;

    double Ytmin = 1.e-30;

    if ( dMpn / T9 > 58. ) {
        Y[1] = Ytmin;
        Y[2] = 1.;
    } else if ( dMpn / T9 < -58. ) {
        Y[1] = 1.;
        Y[2] = Ytmin;
    } else {
        Y[1] = 1./( exp( dMpn/T9) + 1. );
        Y[2] = 1./( exp(-dMpn/T9) + 1. );
    }
    // -->
    Y0[1] = Y[1];
    Y0[2] = Y[2];

    // Calculate the initial baryon density in g/cm^3.
    double rhob0 = GeV4_to_gcm3 * params.eta0 * M_amu * nb_eta_final_ratio(t);
    // -->
    Y[3]  = Y[1]*Y[2]*rhob0*exp(reacparam[12][9]/T9)/pow(T9,1.5)/(reacparam[12][8]*0.987e10);
    Y0[3] = Y[3];
    // Initial statistical equilibrium of deuterium pre-BBN

    for ( i = 4; i <= NNUC; ++i ) {
        Y[i]  = Ytmin;
        Y0[i] = Y[i];
    }

    // Calculate the weak rates (independent of temperature)
    rate_weak(err, f);

    double rho_baryons, dT9_dt, dlnT9_dt;
    while( ltime == 0 ) {
        for( int loop = 1; loop <= 2; loop++ ) {
            // Update the cosmological quantities
            Tnu    = GeV_to_GK * neutrino_temperature(t); // GK
            TGeV   = T9/GeV_to_GK;                        // GeV
            dT9_dt = GeV_to_GK * dTdt(t);                 // GK/s
            // -->
            dlnT9_dt = dT9_dt/T9; // = dlnT_dt, units of s

            rho_baryons = GeV4_to_gcm3 * params.eta0 * M_amu * nb_eta_final_ratio(t); // g/cm^3.

            // Calculate the n->p conversion rate
            rate_pn(err, f, r, T9, Tnu, params.life_neutron);
            // Calculate the nuclear rates
            rate_all(err, f, T9);

            // Linearize the nonlinear ODE system
            fail = linearize(T9, reacparam, f, r, loop, inc, ip, dt, Y0, Y, dY_dt, rho_baryons);
            if ( fail > 0 ) {
                return fail;
            }

            // Condition to break the loop
            // Here dt and 1/dlnT9_dt are in units of s
            // --> 1e-16 is unitless
            if ( T9 <= T9f || dt < fabs(1e-16/dlnT9_dt) || ip == inc ) {
                it++;

                if ( it == nitmax || ip < inc ) {
                    ltime = 1;
                }

                // Calculate the baryon-to-photon ratio
                ratioH[0] = params.eta0 * nb_eta_final_ratio(t)/(0.243587656435 * TGeV*TGeV*TGeV);

                // Calculate the ratios 'element abundance / hydrogen abundance'
                for( int j = 1; j <= NNUC; j++ ) {
                    ratioH[j] = Y[j]/Y[2];
                }

                // Handle the special cases of hydrogen and helium-4
                ratioH[2] = Y[2]*Am[2]; // Am[2] = 1
                ratioH[6] = Y[6]*Am[6]; // Am[6] = 4
            }

            if ( loop == 1 ) {
                // Determine the next time step dt
                if ( ip == inc ) ip = 0;

                ip++;
                is++;

                if ( is > 3 ) {
                    dtmin = fabs(1./dlnT9_dt)*ct; // Units check out
                    for( i = 1; i <= NNUC; i++ ) {
                        if ( dY_dt[i] != 0. && Y[i] > Ytmin ) {
                            dtl = ( fabs(Y[i]/dY_dt[i]) ) * cy * ( pow(log10(Y[i])/log10(Ytmin), 2.) + 1. );

                            if ( dtl < dtmin ) {
                                dtmin = dtl;
                            }
                        }
                    }

                    if ( dtmin > dt*1.5 ) {
                        dtmin = dt*1.5;
                    }

                    // -->
                    dt = dtmin;
                }

                // CHANGED:
                // Enforce a minimal time step
                if ( dt/t < 5.0e-4 ) {
                    dt = 5.0e-4*t;
                }

                // Adapt the time
                t += dt;

                // Calculate the new temperature for the next iteration step
                T9 = GeV_to_GK * temperature(t);

                for( i = 1; i <= NNUC; i++ ) {
                    Y0[i]     = Y[i];
                    dY_dt0[i] = dY_dt[i];
                    Y[i]      = Y0[i] + dY_dt0[i]*dt;

                    if ( Y[i] < Ytmin ) {
                        Y[i] = Ytmin;
                    }
                }
            } else { // if ( loop==2)
                for( i = 1; i <= NNUC; i++ ) {
                    Y[i] = Y0[i] + ( dY_dt[i] + dY_dt0[i] )*0.5*dt;

                    if ( Y[i] < Ytmin ) {
                        Y[i] = Ytmin;
                    }
                }
            }
        } // end of loop over 'loop = 1,2'


    } // while( ltime == 0 )

    // Post processing
    for ( i = 0; i <= NNUC; i++ ) {
        ratioH[i] = fabs(ratioH[i]);
    }

    return fail;
}


void get_final_abundances(int err, Parameters params, double Y0[NNUC+1]) {
    // Define an array to store the ratios
    double ratioH[NNUC+1];

    // Run the BBN code
    run_nucleosynthesis(err, params, ratioH);

    // Extract the baryon-to-photon ratio
    Y0[0] = ratioH[0];

    // Handle the special case of hydrogen
    Y0[2] = ratioH[2];

    for ( int i = 1; i <= NNUC; i++ ) {
        if ( i == 2 ) continue;

        Y0[i] = ratioH[i]*Y0[2];
    }
        
    // Handle the special case of helium-4
    Y0[6] = ratioH[6]/4;

    // Perform the neutron decay
    if ( params.decay_neutrons ) {
        Y0[2] += Y0[1];
        Y0[1]  = 0.;
    }

    /* Result of this function:
        Y0[0]  final baryon-to-photon ratio
        Y0[1]  n_n   / n_b
        Y0[2]  n_H   / n_b
        Y0[3]  n_D   / n_b
        Y0[4]  n_3H  / n_b
        Y0[5]  n_3He / n_b
        Y0[6]  n_4He / n_b = Yp/4
        Y0[7]  n_6Li / n_b
        Y0[8]  n_7Li / n_b
        Y0[9]  n_7Be / n_b
    */

}
