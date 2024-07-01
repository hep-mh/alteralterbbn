#pragma once

#include <math.h>    // sqrt, exp, pow, ...
#include <stdbool.h> // bool
#include <stdio.h>   // printf
#include <stdlib.h>  // malloc, exit
#include <string.h>  // memcpy


// Physical parameters
#define pi      3.1415926535897932385  // --
#define hbar    6.582119514e-25        // GeV s
#define M_amu   0.931494               // GeV, atomic mass unit
#define m_e     510.9989461e-6         // GeV
#define k_B     8.617330e-5            // GeV/GK
#define alphaem 0.007297353            // --
#define dMpn    15.009288310184841     // GeV, mass difference of proton and neutron

// Parameters for unit conversions
#define GeV_to_GK    11604.52165624311
#define GeV4_to_gcm3 2.3201150650611885e+17

// Algorithm specific parameters
#define NNUCREAC 100
#define NNUC      26

// Input specific parameters
#define COSMO_COLS                7

#define COSMO_COL_t               0 // time                 in s
#define COSMO_COL_T               1 // temperature          in GeV
#define COSMO_COL_dTdt            2 // dT/dt                in GeV/s
#define COSMO_COL_Tnu             3 // neutrino temperature in GeV
#define COSMO_COL_nphi            4 // n_phi                in GeV^3
#define COSMO_COL_H               5 // Hubble rate          in 1/s
#define COSMO_COL_nBEtaFinalRatio 6 // nBaryon/etaFinal     in GeV^3


struct parameters {
    double eta0;                // baryon-to-photon ratio
    double life_neutron;        // neutron lifetime
    char *method;               // integration method
    bool decay_neutrons;        // wether to decay neutrons after BBN
};


// cosmo.c //////////////////////////////////////////////////////////////////////////////
extern double **cosmo_array;
extern int COSMO_ROWS;
extern double delta_logx_cosmo;
extern bool cosmo_file_loaded;

typedef enum { ASCENDING, DESCENDING } SortOrder;

SortOrder determine_sort_order(double *arr, int size);
int       find_index(double *arr, int size, double x);
void      load_cosmo_data(char *filename, int ncols);
void      free_cosmo_data();
double    interp_cosmo_data(double x, int xc, int yc);

double    interp_cosmo_array(int i_col, double x);
double    cosmo_t_T(double T);


// util.c ///////////////////////////////////////////////////////////////////////////////
double max(double x, double y);
double min(double x, double y);
double factorial(int n);


// rates.c //////////////////////////////////////////////////////////////////////////////
void   rate_weak(int err, double f[NNUCREAC+1]);
double rate_pn_enu(int type, double T9, double Tnu, int beta_samples, double fierz);
void   rate_pn_noerr(double f[NNUCREAC+1], double r[NNUCREAC+1], double T9, double Tnu, double taun, int beta_samples, double fierz);
void   rate_pn(int err, double f[NNUCREAC+1], double r[NNUCREAC+1], double T9, double Tnu, double taun);
void   rate_all(int err, double f[NNUCREAC+1], double T9);


// test.c TEMP
void rate_weak_test(int err, double f[]);
void rate_all_test(int err, double f[], double T9);


// bbn.c ////////////////////////////////////////////////////////////////////////////////
int  linearize(double T9, double reacparam[NNUCREAC+1][10], double f[NNUCREAC+1], double r[NNUCREAC+1], int loop, int inc, int ip, double dt, double y0[NNUC+1], double y[NNUC+1], double dydt[NNUC+1], double rhob);
int  nucl(int err, struct parameters params, double ratioH[NNUC+1]);
void final_abundances(int err, struct parameters params, double Y0[NNUC+1]);

