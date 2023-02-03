#ifndef NKBB_H_
#define NKBB_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <float.h>

const double Pi = M_PI;
double xobs, yobs;
double epsi3, a13, a22, a52;
double spin, isco;

#include <fitsio.h>

#define NUM_PARAM_NKBB 11
#define PARAM_DEFAULT 0.0

#define TR_TABLE_PATH "/opt/tools/nkbb"

#define FITSTABLE_NDEFPAR 30 							// nonkerr modifications
#define NK_FILENAME_0 "Trf_Kerr_Thermal_Mdl.fits"
#define NK_FILENAME_1 "Trf_Johannsen_a13_thermal.fits"
#define NK_FILENAME_2 "Trf_Johannsen_a22_thermal.fits"
#define NK_FILENAME_3 "Trf_Johannsen_e3_thermal.fits"
#define NK_FILENAME_4 "Trf_conformal_gravity.fits"
#define NK_FILENAME_5 "Trf_KRZ_d1.fits"
#define NK_FILENAME_6 "Trf_KRZ_d2.fits"

#define F_FILENAME_0 "dimlessF_Kerr_Mdl.fits"
#define F_FILENAME_1 "dimlessF_a13.fits"
#define F_FILENAME_2 "dimlessF_a22.fits"
#define F_FILENAME_3 "dimlessF_epsi3.fits"

#define FITSTABLE_NA 30
#define FITSTABLE_NMU0 22
#define FITSTABLE_NR 100
#define FITSTABLE_NG 40

double robs[FITSTABLE_NR], gmin[FITSTABLE_NR], gmax[FITSTABLE_NR], trff1[FITSTABLE_NR][FITSTABLE_NG], trff2[FITSTABLE_NR][FITSTABLE_NG], cosne1[FITSTABLE_NR][FITSTABLE_NG], cosne2[FITSTABLE_NR][FITSTABLE_NG], F[FITSTABLE_NR];

#define TR_ERROR(msg,status) (tr_error(__func__, msg,status))

#define CHECK_TR_ERROR(msg,status) (check_tr_error(__func__, msg,status))

#define CHECK_STATUS_RET(status, retval) \
  if (EXIT_SUCCESS!=status) return(retval);

#define CHECK_STATUS_VOID(status) \
  if (EXIT_SUCCESS!=status) return;

#define CHECK_STATUS_BREAK(status) \
  if (EXIT_SUCCESS!=status) break;

#define CHECK_MALLOC_VOID_STATUS(a,status) \
	if (NULL==a) { \
		TR_ERROR("memory allocation failed",status); \
		return;\
	}

#define CHECK_MALLOC_RET_STATUS(a,status,retval) \
	if (NULL==a) { \
		TR_ERROR("memory allocation failed",status); \
		return retval;\
	}

#define CHECK_STATUS_RET(status, retval) \
  if (EXIT_SUCCESS!=status) return(retval);




typedef struct{

// "eta"
// "a"
// "i"
// "Mbh"
// "Mdd"
// "Dbh"
// "h"
// "rflag"
// "lflag"

	double eta;
	double a;
	double incl;
	double Mbh;
	double Mdd;
	double Dbh;
	double h;
	double rflag;
	double lflag;
	int defpartype;
	double defpar;

} kerrbbParam;

typedef struct{
	double* r;
	double* gmin;
	double* gmax;
	double** trff1;
	double** trff2;
	double** cosne1;
	double** cosne2;
}fitsDat;

// reltable.h
typedef struct{
	float* a; // spin
	int n_a;

	double** defpar; // deformation parameter  // need to be modified yet 
	int n_defpar; 

	float* mu0; // inclination
	int n_mu0;

	fitsDat**** arr; // relline data array

	// dimensions of relline array
	int n_r;
	int n_g;

}fitsTable;


/** the F single data structure */
typedef struct{
	double* F;
}FDat;

/** the F table structure */
typedef struct{

	float* a; // spin
	int n_a;

	double** defpar; // deformation parameter  // need to be modified yet 
	int n_defpar; 

	FDat*** arr; // relline data array

}FTable;



double cache_limit = 1E-8;

void tr_error(const char* const func, const char* const msg, int* status);
void check_tr_error(const char* const func, const char* const msg, int* status);
char* get_tr_table_path( void );
void set_defpar(double defpar);
void set_defpartype(int defpartype);
double get_defpar();
int get_defpartype();
void set_defparval(double defparval);
double get_defparval();
void set_cached_defparval();
char * get_filename();
char * get_filenameF();
char * get_defparcolname();
void set_cached_defpar();
void set_cached_defpartype();
int check_caching_defpar();
void free_fitsTable(fitsTable* tab);
void read_tr_table(char* filename, fitsTable** inp_tab, int* status);
fitsTable* new_fitsTable(int n_a, int n_mu0, int n_r, int n_g, int n_defpar, int* status);


void get_interpolated_fitsTable(double a, double defpar, double mu0, int* status);

double interp_lin_3d_nk_float(double v1, double v2, double v3, double v4, double v5, double v6, double v7, double v8,
  double d, double d1, double d2, double d3, double d4, double d5, double d6, double d7, double d8);

/* binary search */
int binary_search_float(float* arr,int n,float val);




void init_par_nkbb(kerrbbParam** krb_param, const double* inp_par, const int n_parameter, int* status);
void tdnkbb(const double* ener0, const int n_ener0, double* photar, const double* parameter, const int n_parameter, int* status);
void lmodnkbb(const double* ener0, const int n_ener0, const double* parameter, int ifl, double* photar, double* photer, const char* init);

void nkbb_kernel(double* ener0, double* photar, int n_ener0, kerrbbParam* krb_param, int* status );
kerrbbParam* new_kerrbbParam(int* status);
void free_kerrbbParam(kerrbbParam* param);

#endif
