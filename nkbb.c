#include "nkbb.h"



double DEFPAR;
double DEFPARVAL;
int DEFPARTYPE;
double d[9];
double cached_defpar;
double cached_defparval;
int cached_defpartype;

fitsTable* tr_table=NULL;
FTable* F_table=NULL;




void tr_error(const char* const func, const char* const msg, int* status){
	*status = EXIT_FAILURE;
	printf(" *** error in tr (%s): %s!\n", func, msg);
}


void check_tr_error(const char* const func, const char* const msg, int* status){
	if (*status!=EXIT_SUCCESS){
		*status = EXIT_FAILURE;
		printf(" *** error in tr (%s): %s!\n", func, msg);
	}
}


/** get the relxill table path (dynamically from env variable)  **/
char* get_tr_table_path( void ){
	char* path;
	path = getenv("TR_TABLE_PATH");
	if (path!=NULL){
		return path;
	} else {
		return TR_TABLE_PATH;
	}
}




void set_defpar(double defpar) {
	DEFPAR = defpar;
}

void set_defpartype(int defpartype) {
	DEFPARTYPE = defpartype;
}


double get_defpar(){
	return DEFPAR;
}

int get_defpartype(){
	return DEFPARTYPE;
}

void set_defparval(double defparval) {
	DEFPARVAL = defparval;
}

double get_defparval(){
	return DEFPARVAL;
}

void set_cached_defparval() {
	cached_defparval = get_defparval();
}



char * get_filename(){
	if (get_defpartype() == 0) {
		return NK_FILENAME_0;
	}
	if (get_defpartype() == 1) {
		return NK_FILENAME_1;
	}
	if (get_defpartype() == 2) {
		return NK_FILENAME_2;
	}
    if (get_defpartype() == 3) {
        return NK_FILENAME_3;
    }
    if (get_defpartype() == 4) {
        return NK_FILENAME_4;
    }
    if (get_defpartype() == 5) {
        return NK_FILENAME_5;
    }
    if (get_defpartype() == 6) {
        return NK_FILENAME_6;
    }
	return "";
}

char * get_filenameF(){
	if (get_defpartype() == 0) {
		return F_FILENAME_0;
	}
	if (get_defpartype() == 1) {
		return F_FILENAME_1;
	}
	if (get_defpartype() == 2) {
		return F_FILENAME_2;
	}
    if (get_defpartype() == 3) {
        return F_FILENAME_3;
    }
	return "";
}

char * get_defparcolname(){
	if (get_defpartype() == 0) {
		return "Mdl";
	}
	if (get_defpartype() == 1) {
		return "a13";
	}
	if (get_defpartype() == 2) {
		return "a22";
	}
    if (get_defpartype() == 3) {
        return "epsi3";
    }
    if (get_defpartype() == 4) {
        return "L";
    }
    if (get_defpartype() == 5) {
        return "d1";
    }
    if (get_defpartype() == 6) {
        return "d2";
    }
    return "";
}

void set_cached_defpar() {
	cached_defpar = get_defpar();
}

void set_cached_defpartype() {
	cached_defpartype = get_defpartype();
}

int check_caching_defpar() {
	if (fabs(cached_defpar - get_defpar()) > cache_limit || cached_defpartype != get_defpartype())
		return 1;
	else 
		return 0;
}

int check_caching_defpartype() {
	if (cached_defpartype != get_defpartype())
		return 1;
	else 
		return 0;
}



static void cpy(double a[], float b[], int n)
{
	int i;
	for(i = 0; i<n; i++)
	{
		a[i] = (double)b[i];
	}
}



// geometry.c
double POLYHEDRON_VOLUME_3D(double coord[][3], int order_max, int face_num, int node[][6], int node_num, int order[])
{

	// int dim_num = 3;
	int face, n1, n2, n3, v;
	double volume = 0.0;

	for (face = 0; face < face_num; face++)
	{
		n3 = node[order[face] - 1][face];

		for (v = 0; v < order[face] - 2; v++)
		{
			n1 = node[v][face];
			n2 = node[v+1][face];
			volume = volume 
				+ coord[n1][0] * (coord[n2][1] * coord[n3][2] - coord[n3][1] * coord[n2][2])
				+ coord[n2][0] * (coord[n3][1] * coord[n1][2] - coord[n1][1] * coord[n3][2])
				+ coord[n3][0] * (coord[n1][1] * coord[n2][2] - coord[n2][1] * coord[n1][2]);
		}
	}
	volume = volume / 6.0;
	return volume;

}

//	nonkerr modifications
/** get a new and empty nk fits table (structure will be allocated)  */
FTable* new_FTable(int n_a, int n_defpar, int* status){
	FTable* tab = (FTable*) malloc (sizeof(FTable));
	CHECK_MALLOC_RET_STATUS(tab,status,tab);

	// we know which dimensions the table should have
	tab->n_a   =  n_a;
	tab->n_defpar = n_defpar;

	tab->a = NULL;
	tab->defpar = NULL;

	tab->arr=NULL;

	tab->arr = (FDat***) malloc (sizeof(FDat**)*tab->n_a);
	CHECK_MALLOC_RET_STATUS(tab->arr,status,tab);

	int ii; int jj;
	for (ii=0; ii<tab->n_a; ii++){

		tab->arr[ii] = (FDat**) malloc (sizeof(FDat*)*tab->n_defpar);
		CHECK_MALLOC_RET_STATUS(tab->arr[ii],status,tab);

		for (jj=0; jj<tab->n_defpar; jj++){
			tab->arr[ii][jj] = NULL;
			CHECK_STATUS_RET(*status,tab);
		}

	}

	// printf(" new NK object is initialized\n");
	return tab;
}


static FDat* new_FDat(int nr, int* status){
	FDat* dat = (FDat*) malloc (sizeof(FDat));
	CHECK_MALLOC_RET_STATUS(dat,status,dat);

	dat->F = (double*) malloc( sizeof(double) * nr);
	CHECK_MALLOC_RET_STATUS(dat->F,status,dat);

	return dat;
}


// static void free_FDat(FDat* dat, int nr){
// 	if (dat!=NULL){
// 		free(dat->F);
// 	}
// }


void free_FTable(FTable* tab){
	if(tab!=NULL){
		if (tab->arr!=NULL){
			int ii;
			for (ii=0; ii<tab->n_a; ii++){
				if (tab->arr[ii] !=NULL){
					int jj;
					for (jj=0; jj<tab->n_defpar; jj++){
						free(tab->arr[ii][jj]->F);
						free(tab->arr[ii][jj]);
					}
					free(tab->arr[ii]);
				}
			}
			free(tab->arr);
		}
		free(tab->a);
		free(tab->defpar);
		free(tab);
	}
}




/** load one single data extension from the relline table   */
static FDat* load_single_FDat(fitsfile* fptr, char* extname, int nhdu, int* status){

	// int extver = 0;
	// fits_movnam_hdu(fptr, BINARY_TBL, extname, extver ,status);
	int exttype;
	fits_movabs_hdu(fptr,nhdu,&exttype, status);
	if (*status!=EXIT_SUCCESS){
		printf(" *** error moving to extension %s\n",extname);
		return NULL;
	}

	// get the column id-number
	int colnum_F;
	if(fits_get_colnum(fptr, CASEINSEN, extname, &colnum_F, status)) return NULL;
	

	// check the number of rows (need to coincide with FITSTABLE_NR
	// long n;
	// if (fits_get_num_rows(fptr, &n, status)) return NULL;
	// if (n != FITSTABLE_NR){
	// 	TR_ERROR("inconsistent number of rows in tr table",status);
	// 	printf("    -> expecting %i, but found %ld in extensions %s",FITSTABLE_NR,n,extname);
	// 	return NULL;
	// }

	// allocate the memory for the table
	FDat* dat = new_FDat(FITSTABLE_NR, status);
	CHECK_STATUS_RET(*status,NULL);

	// now load the table column by column

    // (1) start with the 1D columns
    int anynul=0;
    double nullval=0.0;
    LONGLONG nelem = (LONGLONG) FITSTABLE_NR;
    float vall[FITSTABLE_NR];
    fits_read_col(fptr, TFLOAT, colnum_F, 1, 1, nelem ,&nullval,vall, &anynul, status);
    cpy(dat->F, vall, FITSTABLE_NR);
    CHECK_STATUS_RET(*status,dat);
    // fits_read_col(fptr, TFLOAT, colnum_gmin, 1, 1, nelem ,&nullval,vall, &anynul, status);
    // cpy(dat->gmin, vall, FITSTABLE_NR);
    // CHECK_STATUS_RET(*status,dat);
    // fits_read_col(fptr, TFLOAT, colnum_gmax, 1, 1, nelem ,&nullval,vall, &anynul, status);
    // cpy(dat->gmax, vall, FITSTABLE_NR);
    // CHECK_STATUS_RET(*status,dat);

    // (2) and finally the 2D columns

    // load_single_fitsDat_2dcol(fptr,dat->trff2,FITSTABLE_NR,FITSTABLE_NG,colnum_trff2,status);
    // CHECK_STATUS_RET(*status,dat);
    // load_single_fitsDat_2dcol(fptr,dat->trff1,FITSTABLE_NR,FITSTABLE_NG,colnum_trff1,status);
    // CHECK_STATUS_RET(*status,dat);
    // load_single_fitsDat_2dcol(fptr,dat->cosne1,FITSTABLE_NR,FITSTABLE_NG,colnum_cosne1,status);
    // CHECK_STATUS_RET(*status,dat);
    // load_single_fitsDat_2dcol(fptr,dat->cosne2,FITSTABLE_NR,FITSTABLE_NG,colnum_cosne2,status);
    // CHECK_STATUS_RET(*status,dat);

    return dat;

}




void read_Ftable(char* filename, FTable** inp_tab, int* status){

	FTable* tab = (*inp_tab);
	fitsfile *fptr=NULL;

	char* fullfilename=NULL;
	char* extname=NULL;

	// printf(" point 4\n");

	do{ // Errot handling loop
		if (tab != NULL && !check_caching_defpartype()){
			TR_ERROR("tr table already loaded",status);
			break;
		}

		if (tab != NULL){
			free_FTable(tab);
			tab = NULL;
		}


		// tab = new_fitsTable(FITSTABLE_NA,FITSTABLE_NMU0,FITSTABLE_NR,FITSTABLE_NG,FITSTABLE_NDEFPAR,status);
		tab = new_FTable(FITSTABLE_NA,FITSTABLE_NDEFPAR,status);
		CHECK_STATUS_BREAK(*status);

		// should be set by previous routine
		assert(tab!=NULL);
		assert(tab->arr!=NULL);

		// get the full filename
		if (asprintf(&fullfilename, "%s/%s", get_tr_table_path() , get_filenameF()) == -1){
			TR_ERROR("failed to construct full path of the rel table",status);
			break;
		}

		// open the file
		if (fits_open_table(&fptr, fullfilename, READONLY, status)) {
			CHECK_TR_ERROR("opening of the F table failed",status);
			printf("    either the full path given (%s) is wrong \n",fullfilename);
			// printf("    or you need to download the table ** %s **  from \n",filename);
			// printf("    http://www.physics.fudan.edu.cn/tps/people/bambi/Site/RELXILL_NK.html \n");
			break;
		}




		// nonkerr modifications  *** yet to be modified ***

		// printf("ready to read\n");


		int _nhdu_ = 2;
		int exttype;
		fits_movabs_hdu(fptr, _nhdu_, &exttype, status);
		if (*status!=EXIT_SUCCESS){
			printf(" *** error moving to spin and defpar table %s\n",extname);
			return;
		}
		int colnum_a, colnum_defpar;
		long _n_;
		if(fits_get_colnum(fptr, CASEINSEN, "a", &colnum_a, status)) { printf(" error at fits get colnum a"); break; }
		// printf(" colnum_a -- %d\n", colnum_a);

		if(fits_get_colnum(fptr, CASEINSEN, get_defparcolname(), &colnum_defpar, status)) { printf(" error at fits get colnum defpar"); break; }				// *** need to be checked more
		// printf(" colnum_defpar -- %d\n", colnum_defpar);

		if (fits_get_num_rows(fptr, &_n_, status)) return;

		if (_n_ != FITSTABLE_NA){
			// printf(" %d %d\n", _n_, FITSTABLE_NA);
			TR_ERROR("inconsistent number of rows in spin table",status);
			return;
		}

		// printf(" length of a, defpar tables -- %d\n", (int)_n_);

		tab->a = (double*) malloc(tab->n_a * sizeof(double));

		int anynul=0;
		double nullval=0.0;
		LONGLONG nelem = (LONGLONG) _n_;
		float tmp[FITSTABLE_NA];
		fits_read_col(fptr, TFLOAT, colnum_a, 1, 1, nelem ,&nullval,tmp, &anynul, status);
		cpy(tab->a, tmp, FITSTABLE_NA);

		if (*status!=EXIT_SUCCESS){
			printf(" *** error reading spin tables %s\n",extname);
			return;
		}

		tab->defpar = (double**) malloc( sizeof(double*) * tab->n_a);
		int iii;
		for (iii = 0; iii < tab->n_a; iii++){
			tab->defpar[iii] = (double*) malloc( sizeof(double) * tab->n_defpar);
		}

		anynul=0;
		nullval=0.0;
   		nelem = (LONGLONG) FITSTABLE_NDEFPAR;
   		float vall[FITSTABLE_NDEFPAR];

		for (iii=0; iii<FITSTABLE_NA;iii++){
	        	if(fits_read_col(fptr, TFLOAT, colnum_defpar, iii+1, 1, nelem ,&nullval,vall, &anynul, status)) return;
	        	cpy(tab->defpar[iii], vall, FITSTABLE_NDEFPAR);
		}

		if (*status!=EXIT_SUCCESS){
			printf(" *** error reading defpar table %s\n",extname);
			return;
		}

		_nhdu_++;
		int ii; int jj;
		for (ii=0; ii<tab->n_a; ii++){
			for (jj=0; jj<tab->n_defpar; jj++){

				if (asprintf(&extname, "F_%i_%i", ii+1,jj+1) == -1){
					TR_ERROR("failed to construct full path the F table",status);
					break;
				}

				assert(tab->arr[ii][jj]==NULL);
				// int nhdu = (ii)*tab->n_mu0+jj+4;
				tab->arr[ii][jj] = new_FDat(FITSTABLE_NR, status);
				tab->arr[ii][jj] = load_single_FDat(fptr, extname, _nhdu_, status);
				free(extname);
				if (*status!=EXIT_SUCCESS){
					TR_ERROR("failed to load data from the rel table into memory",status);
					break;
				}
			}
		}



	} while(0);

	// printf(" point 5\n");

	if (*status==EXIT_SUCCESS){
		// assigne the value
		(*inp_tab) = tab;
	} else {
		free_FTable(tab);
	}
	free(fullfilename);

	if (fptr!=NULL) {fits_close_file(fptr,status);}
	// printf(" point 6\n");

	return;
}





static void load_single_fitsDat_2dcol(fitsfile* fptr, double** val,int n1, int n2, int colnum, int* status){

    int anynul=0;
    double nullval=0.0;

    assert(val!=NULL);

    LONGLONG nelem = (LONGLONG) n2;
    float tmp[FITSTABLE_NG];

    int ii;
    for (ii=0; ii<n1;ii++){
        if(fits_read_col(fptr, TFLOAT, colnum, ii+1, 1, nelem ,&nullval,tmp, &anynul, status)) return;
        // int jj;
        // for (jj = 0; jj < n2; jj++)
        // 	// if (isfinite(tmp[jj]) || isinf(tmp[jj])) { // == NAN) { // || tmp[jj] == INFINITY) {
        // 	if (colnum == 6 || colnum == 7) {
        // 		if (tmp[jj] > 1.0)
        // 		{
        // 			// printf(" cosne[%d][%d][%d] = %f\n", colnum - 5, ii, jj, tmp[jj]);
        // 			tmp[jj] = 0.998;
        // 			n_cosnes++;
        // 		}
        // 	}
        // 	if (tmp[jj] != tmp[jj] || tmp[jj] == INFINITY) {
        // 		//printf("%f\n", tmp[jj]);
        // 		n_inf_nans = n_inf_nans + 1;
        // 		tmp[jj] = 0.;
        		

        // 	}
        cpy(val[ii], tmp, FITSTABLE_NG);

    }
}




//	nonkerr modifications
/** get a new and empty nk fits table (structure will be allocated)  */
fitsTable* new_fitsTable(int n_a, int n_mu0, int n_r, int n_g, int n_defpar, int* status){
	fitsTable* tab = (fitsTable*) malloc (sizeof(fitsTable));
	CHECK_MALLOC_RET_STATUS(tab,status,tab);

	// we know which dimensions the table should have
	tab->n_a   =  n_a;
	tab->n_mu0 =  n_mu0;
	tab->n_r   =  n_r;
	tab->n_g   =  n_g;
	tab->n_defpar = n_defpar;

	tab->a = NULL;
	tab->mu0 = NULL;
	tab->defpar = NULL;

	tab->arr=NULL;

	tab->arr = (fitsDat****) malloc (sizeof(fitsDat***)*tab->n_a);
	CHECK_MALLOC_RET_STATUS(tab->arr,status,tab);

	int ii; int jj; int kk;
	for (ii=0; ii<tab->n_a; ii++){

		tab->arr[ii] = (fitsDat***) malloc (sizeof(fitsDat**)*tab->n_defpar);
		CHECK_MALLOC_RET_STATUS(tab->arr[ii],status,tab);

		for (jj=0; jj<tab->n_defpar; jj++){
			tab->arr[ii][jj] = (fitsDat**) malloc (sizeof(fitsDat*)*tab->n_mu0);
			CHECK_MALLOC_RET_STATUS(tab->arr[ii][jj],status,tab);

			for(kk = 0; kk < tab->n_mu0; kk++){
				tab->arr[ii][jj][kk] = NULL;
				CHECK_STATUS_RET(*status,tab);
			}
		}

	}

	// printf(" new NK object is initialized\n");
	return tab;
}





static fitsDat* new_fitsDat(int nr, int ng, int* status){
	fitsDat* dat = (fitsDat*) malloc (sizeof(fitsDat));
	CHECK_MALLOC_RET_STATUS(dat,status,dat);

	dat->r = (double*) malloc( sizeof(double) * nr);
	CHECK_MALLOC_RET_STATUS(dat->r,status,dat);
	dat->gmin = (double*) malloc( sizeof(double) * nr);
	CHECK_MALLOC_RET_STATUS(dat->gmin,status,dat);
	dat->gmax = (double*) malloc( sizeof(double) * nr);
	CHECK_MALLOC_RET_STATUS(dat->gmax,status,dat);

	dat->cosne1 = (double**) malloc( sizeof(double*) * nr);
	CHECK_MALLOC_RET_STATUS(dat->cosne1,status,dat);
	dat->cosne2 = (double**) malloc( sizeof(double*) * nr);
	CHECK_MALLOC_RET_STATUS(dat->cosne2,status,dat);
	dat->trff1 = (double**) malloc( sizeof(double*) * nr);
	CHECK_MALLOC_RET_STATUS(dat->trff1,status,dat);
	dat->trff2 = (double**) malloc( sizeof(double*) * nr);
	CHECK_MALLOC_RET_STATUS(dat->trff2,status,dat);

	int ii;
	for (ii=0; ii<nr; ii++){
		dat->cosne1[ii] = (double*) malloc( sizeof(double) * ng);
		CHECK_MALLOC_RET_STATUS(dat->cosne1,status,dat);
		dat->cosne2[ii] = (double*) malloc( sizeof(double) * ng);
		CHECK_MALLOC_RET_STATUS(dat->cosne2,status,dat);
		dat->trff1[ii] = (double*) malloc( sizeof(double) * ng);
		CHECK_MALLOC_RET_STATUS(dat->trff1,status,dat);
		dat->trff2[ii] = (double*) malloc( sizeof(double) * ng);
		CHECK_MALLOC_RET_STATUS(dat->trff2,status,dat);
	}
	return dat;
}







/** load one single data extension from the relline table   */
static fitsDat* load_single_fitsDat(fitsfile* fptr, char* extname, int nhdu, int* status){

	// int extver = 0;
	// fits_movnam_hdu(fptr, BINARY_TBL, extname, extver ,status);
	int exttype;
	fits_movabs_hdu(fptr,nhdu,&exttype, status);
	if (*status!=EXIT_SUCCESS){
		printf(" *** error moving to extension %s\n",extname);
		return NULL;
	}

	// get the column id-number
	int colnum_r,colnum_gmin,colnum_gmax;
	int colnum_trff1, colnum_trff2;
	int colnum_cosne1, colnum_cosne2;
	if(fits_get_colnum(fptr, CASEINSEN, "r", &colnum_r, status)) return NULL;
	if(fits_get_colnum(fptr, CASEINSEN, "gmin", &colnum_gmin, status)) return NULL;
	if(fits_get_colnum(fptr, CASEINSEN, "gmax", &colnum_gmax, status)) return NULL;
	if(fits_get_colnum(fptr, CASEINSEN, "trff1", &colnum_trff1, status)) return NULL;
	if(fits_get_colnum(fptr, CASEINSEN, "trff2", &colnum_trff2, status)) return NULL;
	if(fits_get_colnum(fptr, CASEINSEN, "cosne1", &colnum_cosne1, status)) return NULL;
	if(fits_get_colnum(fptr, CASEINSEN, "cosne2", &colnum_cosne2, status)) return NULL;

	// check the number of rows (need to coincide with FITSTABLE_NR
	long n;
	if (fits_get_num_rows(fptr, &n, status)) return NULL;
	if (n != FITSTABLE_NR){
		TR_ERROR("inconsistent number of rows in tr table",status);
		printf("    -> expecting %i, but found %ld in extensions %s",FITSTABLE_NR,n,extname);
		return NULL;
	}

	// allocate the memory for the table
	fitsDat* dat = new_fitsDat(FITSTABLE_NR,FITSTABLE_NG,status);
	CHECK_STATUS_RET(*status,NULL);

	// now load the table column by column

    // (1) start with the 1D columns
    int anynul=0;
    double nullval=0.0;
    LONGLONG nelem = (LONGLONG) FITSTABLE_NR;
    float vall[FITSTABLE_NR];
    fits_read_col(fptr, TFLOAT, colnum_r, 1, 1, nelem ,&nullval,vall, &anynul, status);
    cpy(dat->r, vall, FITSTABLE_NR);
    CHECK_STATUS_RET(*status,dat);
    fits_read_col(fptr, TFLOAT, colnum_gmin, 1, 1, nelem ,&nullval,vall, &anynul, status);
    cpy(dat->gmin, vall, FITSTABLE_NR);
    CHECK_STATUS_RET(*status,dat);
    fits_read_col(fptr, TFLOAT, colnum_gmax, 1, 1, nelem ,&nullval,vall, &anynul, status);
    cpy(dat->gmax, vall, FITSTABLE_NR);
    CHECK_STATUS_RET(*status,dat);

    // (2) and finally the 2D columns

    load_single_fitsDat_2dcol(fptr,dat->trff2,FITSTABLE_NR,FITSTABLE_NG,colnum_trff2,status);
    CHECK_STATUS_RET(*status,dat);
    load_single_fitsDat_2dcol(fptr,dat->trff1,FITSTABLE_NR,FITSTABLE_NG,colnum_trff1,status);
    CHECK_STATUS_RET(*status,dat);
    load_single_fitsDat_2dcol(fptr,dat->cosne1,FITSTABLE_NR,FITSTABLE_NG,colnum_cosne1,status);
    CHECK_STATUS_RET(*status,dat);
    load_single_fitsDat_2dcol(fptr,dat->cosne2,FITSTABLE_NR,FITSTABLE_NG,colnum_cosne2,status);
    CHECK_STATUS_RET(*status,dat);

    return dat;

}

static void free_fitsDat(fitsDat* dat, int nr){
	if (dat!=NULL){
		int ii;
		for (ii=0; ii<nr; ii++){
			if (dat->cosne1 !=NULL) free(dat->cosne1[ii]);
			if (dat->cosne2 !=NULL) free(dat->cosne2[ii]);
			if (dat->trff1 !=NULL) free(dat->trff1[ii]);
			if (dat->trff2 !=NULL) free(dat->trff2[ii]);
		}
		free(dat->cosne1);
		free(dat->cosne2);
		free(dat->trff1);
		free(dat->trff2);

		free(dat->r);
		free(dat->gmin);
		free(dat->gmax);
	}
}




void free_fitsTable(fitsTable* tab){
	if(tab!=NULL){
		if (tab->arr!=NULL){
			int ii;
			for (ii=0; ii<tab->n_a; ii++){
				if (tab->arr[ii] !=NULL){
					int jj;
					for (jj=0; jj<tab->n_defpar; jj++){
						if(tab->arr[ii][jj] != NULL){
							int kk;
							for(kk=0; kk<tab->n_mu0; kk++){
								free_fitsDat(tab->arr[ii][jj][kk],tab->n_r);
								free(tab->arr[ii][jj][kk]);
							}
							free(tab->arr[ii][jj]);
						}
					}
					free(tab->arr[ii]);
				}
			}
			free(tab->arr);
		}
		free(tab->a);
		free(tab->mu0);
		free(tab);
	}
}










void read_tr_table(char* filename, fitsTable** inp_tab, int* status){

	fitsTable* tab = (*inp_tab);
	fitsfile *fptr=NULL;

	char* fullfilename=NULL;
	char* extname=NULL;

	// printf(" point 1\n");

	do{ // Errot handling loop
		if (tab != NULL && !check_caching_defpartype()){
			TR_ERROR("tr table already loaded",status);
			break;
		}

		if (tab != NULL){
			free_fitsTable(tab);
			tab = NULL;
		}


		tab = new_fitsTable(FITSTABLE_NA,FITSTABLE_NMU0,FITSTABLE_NR,FITSTABLE_NG,FITSTABLE_NDEFPAR,status);
		CHECK_STATUS_BREAK(*status);

		// should be set by previous routine
		assert(tab!=NULL);
		assert(tab->arr!=NULL);

		// get the full filename
		if (asprintf(&fullfilename, "%s/%s", get_tr_table_path() , get_filename()) == -1){
			TR_ERROR("failed to construct full paththe rel table",status);
			break;
		}

		// open the file
		if (fits_open_table(&fptr, fullfilename, READONLY, status)) {
			CHECK_TR_ERROR("opening of the tr table failed",status);
			printf("    either the full path given (%s) is wrong \n",fullfilename);
			// printf("    or you need to download the table ** %s **  from \n",filename);
			// printf("    http://www.physics.fudan.edu.cn/tps/people/bambi/Site/RELXILL_NK.html \n");
			break;
		}




		// nonkerr modifications  *** yet to be modified ***

		// printf("ready to read\n");


		int _nhdu_ = 2;
		int exttype;
		fits_movabs_hdu(fptr, _nhdu_, &exttype, status);
		if (*status!=EXIT_SUCCESS){
			printf(" *** error moving to spin and defpar table %s\n",extname);
			return;
		}
		int colnum_a, colnum_defpar;
		long _n_;
		if(fits_get_colnum(fptr, CASEINSEN, "a", &colnum_a, status)) { printf(" error at fits get colnum a"); break; }
		// printf(" colnum_a -- %d\n", colnum_a);

		if(fits_get_colnum(fptr, CASEINSEN, get_defparcolname(), &colnum_defpar, status)) { printf(" error at fits get colnum defpar"); break; }				// *** need to be checked more
		// printf(" colnum_defpar -- %d\n", colnum_defpar);

		if (fits_get_num_rows(fptr, &_n_, status)) return;

		if (_n_ != FITSTABLE_NA){
			// printf(" %d %d\n", _n_, FITSTABLE_NA);
			TR_ERROR("inconsistent number of rows in spin table",status);
			return;
		}

		// printf(" length of a, defpar tables -- %d\n", (int)_n_);

		tab->a = (float*) malloc(tab->n_a * sizeof(float));

		int anynul=0;
		double nullval=0.0;
		LONGLONG nelem = (LONGLONG) _n_;
		fits_read_col(fptr, TFLOAT, colnum_a, 1, 1, nelem ,&nullval,tab->a, &anynul, status);

		if (*status!=EXIT_SUCCESS){
			printf(" *** error reading spin tables %s\n",extname);
			return;
		}

		tab->defpar = (double**) malloc( sizeof(double*) * tab->n_a);
		int iii;
		for (iii = 0; iii < tab->n_a; iii++){
			tab->defpar[iii] = (double*) malloc( sizeof(double) * tab->n_defpar);
		}

		anynul=0;
		nullval=0.0;
   		nelem = (LONGLONG) FITSTABLE_NDEFPAR;
   		float vall[FITSTABLE_NDEFPAR];

		for (iii=0; iii<FITSTABLE_NA;iii++){
	        	if(fits_read_col(fptr, TFLOAT, colnum_defpar, iii+1, 1, nelem ,&nullval,vall, &anynul, status)) return;
	        	cpy(tab->defpar[iii], vall, FITSTABLE_NDEFPAR);
		}

		if (*status!=EXIT_SUCCESS){
			printf(" *** error reading defpar table %s\n",extname);
			return;
		}

	    _nhdu_ = _nhdu_ + 1;
		fits_movabs_hdu(fptr, _nhdu_, &exttype, status);
		if (*status!=EXIT_SUCCESS){
			printf(" *** error moving to mu0 table %s\n",extname);
			return;
		}
		int colnum_mu0;
		if(fits_get_colnum(fptr, CASEINSEN, "mu0", &colnum_mu0, status)) { printf(" error at fits get colnum mu0"); break; }
		
		// printf(" colnum_mu0 -- %d\n", colnum_mu0);

		if (fits_get_num_rows(fptr, &_n_, status)) return;
		if (_n_ != FITSTABLE_NMU0){
			TR_ERROR("inconsistent number of rows in mu0 table",status);
			return;
		}

		tab->mu0 = (float*) malloc(tab->n_mu0 * sizeof(float));
	
		anynul=0;
		nullval=0.0;
		nelem = (LONGLONG) _n_;
		//fits_read_col(fptr, TFLOAT, colnum, 1, 1, nelem ,&nullval,*val, &anynul, status);
		fits_read_col(fptr, TFLOAT, colnum_mu0, 1, 1, nelem ,&nullval,tab->mu0, &anynul, status);

		if (*status!=EXIT_SUCCESS){
			printf(" *** error reading mu0 table %s\n",extname);
			return;
		}
	

		//now load the full table (need to go through all extensions)
		int ii; int jj; int kk; int nhdu = 4;
		for (ii=0; ii<tab->n_a; ii++){
			for (jj=0; jj < tab->n_defpar; jj++){    								// ????????
				 for(kk = 0; kk < tab->n_mu0; kk++) {

					if (asprintf(&extname, "%i_%i_%i", ii+1,jj+1,kk+1) == -1){ 									// 		?????? extname  - where is it used? *** found
						TR_ERROR("failed to construct full path the rel nk table",status);
						break;
					}

					assert(tab->arr[ii][jj][kk]==NULL); 												// nonkerr modifications
					//int nhdu = (ii)*tab->n_defpar + (jj)*tab->n_mu0 + kk + 4; 							// nonkerr modifications
					tab->arr[ii][jj][kk] = load_single_fitsDat(fptr, extname, nhdu, status);
					nhdu = nhdu + 1;
					free(extname);
					if (*status!=EXIT_SUCCESS){
						TR_ERROR("failed to load data from the rel table into memory",status);
						break;
					}	
				}
			}
		}

	} while(0);

	// printf(" point 2\n");

	if (*status==EXIT_SUCCESS){
		// assigne the value
		(*inp_tab) = tab;
	} else {
		free_fitsTable(tab);
	}
	free(fullfilename);

	if (fptr!=NULL) {fits_close_file(fptr,status);}

	return;
}





static void interpol_elements(fitsTable* fitstab, int ind_a, int ind_mu0, double ifac_a, double ifac_mu0, double defpar, double x1, double x2, double x1l, double x1r, double x2l, double x2r){
//-----------------------------------------------------------------------------------------------------------------------------------//
//---------  In this function we calculate the interpolation kernel. Results are assigned to the global variable d[9]. --SN ---------//
//-----------------------------------------------------------------------------------------------------------------------------------//
    int order_max = 4, face_num = 6, node_num = 8;
    int order[] = {4, 4, 4, 4, 4, 4};
	int node1[4][6] = {{2, 0, 0, 7, 5, 7}, /* 1 */{6, 1, 3, 5, 6, 1}, /* 2 */ {1, 7, 4, 4, 2, 6}, /* 3 */ {0, 3, 2, 3, 4, 5} /* 4 */};
    
    double crd[8][3] = { \
		{0., 0., x2l},     /* 1 */ \
		{1., 0., x2r}, /* 2 */ \
		{0., 1., x2l},     /* 3 */ \
		{0., 0., x1l}, \
		{0., 1., x1l}, /* 5 */ \
		{1., 1., x1r}, /* 6 */ \
		{1., 1., x2r}, /* 7 */ \
		{1., 0., x1r} /* 8 */ \
	};
	//defparvoltotal = vol_interpol(coord, order_max, face_num, node, node_num, order);
	d[0] = POLYHEDRON_VOLUME_3D(crd, order_max, face_num, node1, node_num, order);
    //printf("dtotal = %f\n", d[0]);

	double crd1[8][3] = { \
		{ifac_a, ifac_mu0, defpar},
		{1., ifac_mu0, defpar},
		{ifac_a, 1., defpar},
		{ifac_a, ifac_mu0, x1},
		{ifac_a, 1., x1},
		{1., 1., x1r},
		{1., 1., defpar},
		{1., ifac_mu0, x1r}
	};
	// defparvol1 = vol_interpol(coord1, order_max, face_num, node, node_num, order);
	d[1] = POLYHEDRON_VOLUME_3D(crd1, order_max, face_num, node1, node_num, order);
	//printf("%f\n", d1);

	double crd2[8][3] = {
		{0., ifac_mu0, defpar},
		{ifac_a, ifac_mu0, defpar},
		{0., 1., defpar},
		{0., ifac_mu0, x1l},
		{0., 1., x1l},
		{ifac_a, 1., x1},
		{ifac_a, 1., defpar},
		{ifac_a, ifac_mu0, x1}
	};
	// defparvol2 = vol_interpol(coord2, order_max, face_num, node, node_num, order);
	d[2] = POLYHEDRON_VOLUME_3D(crd2, order_max, face_num, node1, node_num, order);
	//printf("%f\n", d2);

	double crd3[8][3] = {
		{ifac_a, 0., defpar},
		{1., 0., defpar},
		{ifac_a, ifac_mu0, defpar},
		{ifac_a, 0., x1},
		{ifac_a, ifac_mu0, x1},
		{1., ifac_mu0, x1r},
		{1., ifac_mu0, defpar},
		{1., 0., x1r}
	};
	// defparvol3 = vol_interpol(coord3, order_max, face_num, node, node_num, order);
	d[3] = POLYHEDRON_VOLUME_3D(crd3, order_max, face_num, node1, node_num, order);
	//printf("%f\n", d3);

	double crd4[8][3] = {
		{ifac_a, ifac_mu0, x2},
		{1., ifac_mu0, x2r},
		{ifac_a, 1., x2},
		{ifac_a, ifac_mu0, defpar},
		{ifac_a, 1., defpar},
		{1., 1., defpar},
		{1., 1., x2r},
		{1., ifac_mu0, defpar}
	};
	// defparvol4 = vol_interpol(coord4, order_max, face_num, node, node_num, order);
	d[4] = POLYHEDRON_VOLUME_3D(crd4, order_max, face_num, node1, node_num, order);
	//printf("%f\n", d4);

	double crd5[8][3] = {
		{ifac_a, 0., x2},
		{1., 0., x2r},
		{ifac_a, ifac_mu0, x2},
		{ifac_a, 0., defpar},
		{ifac_a, ifac_mu0, defpar},
		{1., ifac_mu0, defpar},
		{1., ifac_mu0, x2r},
		{1., 0., defpar}
	};
	// defparvol5 = vol_interpol(coord5, order_max, face_num, node, node_num, order);
	d[5] = POLYHEDRON_VOLUME_3D(crd5, order_max, face_num, node1, node_num, order);	
	//printf("%f\n", d5);

	double crd6[8][3] = {
		{0., 0., x2l},
		{ifac_a, 0., x2},
		{0., ifac_mu0, x2l},
		{0., 0., defpar},
		{0., ifac_mu0, defpar},
		{ifac_a, ifac_mu0, defpar},
		{ifac_a, ifac_mu0, x2},
		{ifac_a, 0., defpar}
	};
	// defparvol6 = vol_interpol(coord6, order_max, face_num, node, node_num, order);
	d[6] = POLYHEDRON_VOLUME_3D(crd6, order_max, face_num, node1, node_num, order);
    //printf("%f\n", d6);

	double crd7[8][3] = {
		{0., 0., defpar},
		{ifac_a, 0., defpar},
		{0., ifac_mu0, defpar},
		{0., 0., x1l},
		{0., ifac_mu0, x1l},
		{ifac_a, ifac_mu0, x1},
		{ifac_a, ifac_mu0, defpar},
		{ifac_a, 0., x1}
	};
	// defparvol7 = vol_interpol(coord7, order_max, face_num, node, node_num, order);
	d[7] = POLYHEDRON_VOLUME_3D(crd7, order_max, face_num, node1, node_num, order);
    //printf("%f\n", d7);

	double crd8[8][3] = {
		{0., ifac_mu0, x2l},
		{ifac_a, ifac_mu0, x2},
		{0., 1., x2l},
		{0., ifac_mu0, defpar},
		{0., 1., defpar},
		{ifac_a, 1., defpar},
		{ifac_a, 1., x2},
		{ifac_a, ifac_mu0, defpar}
	};
	// defparvol8 = vol_interpol(coord8, order_max, face_num, node, node_num, order);
	d[8] = POLYHEDRON_VOLUME_3D(crd8, order_max, face_num, node1, node_num, order);
    //printf("%lf %f %f %f %f %f %f %f %f\n", d[1], d[2], d[3], d[4], d[5], d[6], d[7], d[8], (d[1]+d[2]+d[3]+d[4]+d[5]+d[6]+d[7]+d[8]));
    
    return;
}


/** linear interpolation in 1 dimension **/
double interp_lin_1d(double ifac_r, double rlo, double rhi){
	return ifac_r*rhi + (1.0-ifac_r)*rlo;
}




// static void interpol_a_defpar_mu0(double a, double defpar, double mu0, relSysPar* sysPar, fitsTable* fitstab, double rin, double rout) {
static void interpol_a_defpar_mu0(double a, double defpar, double mu0, fitsTable* fitstab, FTable* Ftab) {
//-----------------------------------------------------------------------------------------------------------------------------------//
//---------  In this function we interpolate r, gmin, gmax, trff, and cosem. --SN ---------//
//-----------------------------------------------------------------------------------------------------------------------------------//
    int ii, jj, kk;
    // int ind_rmin = fitstab->n_r - 1, ind_rmax = 0;
    double unscaled_defpar_left_lim, unscaled_defpar_right_lim, unscaled_defpar_left, unscaled_defpar_right;
    int ind_a, ind_mu0, idefl=-1, idefr=-1, idl1, idl2, idr1, idr2;
    double ifac_a, ifac_mu0, x1, x2, x1l, x1r, x2l, x2r;
    int idefl1, idefl2, idefr1, idefr2;
    double ifac_dp1, ifac_dp2, tdp1, tdp2;
    
    /*printf("mu0 = %f\tincl = %f\n",mu0, acos(mu0)*180.0/M_PI);*/

    ind_a   = binary_search_float(fitstab->a,fitstab->n_a,a);
	ind_mu0 = binary_search_float(fitstab->mu0,fitstab->n_mu0,mu0);
    ifac_a   = (a-fitstab->a[ind_a])/(fitstab->a[ind_a+1]-fitstab->a[ind_a]);
	ifac_mu0 = (mu0-fitstab->mu0[ind_mu0])/(fitstab->mu0[ind_mu0+1]-fitstab->mu0[ind_mu0]);
    
	if(defpar < 0.0) {
        unscaled_defpar_left_lim = fitstab->defpar[ind_a][0];
        unscaled_defpar_right_lim = fitstab->defpar[ind_a+1][0];
        unscaled_defpar_left = (-1.0)*defpar*unscaled_defpar_left_lim;
        unscaled_defpar_right = (-1.0)*defpar*unscaled_defpar_right_lim;
    }
    else{
        unscaled_defpar_left_lim = fitstab->defpar[ind_a][29];
        unscaled_defpar_right_lim = fitstab->defpar[ind_a+1][29];
        unscaled_defpar_left = (1.0)*defpar*unscaled_defpar_left_lim;
        unscaled_defpar_right = (1.0)*defpar*unscaled_defpar_right_lim;
    }
    
    for (ii = 0; ii < FITSTABLE_NDEFPAR; ii++) {
		if (fitstab->defpar[ind_a][ii] <= unscaled_defpar_left) {
			idefl = ii;
        }
		if (fitstab->defpar[ind_a+1][ii] <= unscaled_defpar_right) {
			idefr = ii;
        }
    }
    
    //printf("INTERPOLATION INDEX: %d, %d, %d\n", ind_a, ind_mu0, idefl);
    
    if(idefl >= 0 && idefl < (FITSTABLE_NDEFPAR-1) && idefr >= 0 && idefr < (FITSTABLE_NDEFPAR-1)) {
        idl1 = 0; idl2 = 0; idr1 = 0; idr2 = 0;
    }
    else if(idefl == (FITSTABLE_NDEFPAR-1) && idefr == (FITSTABLE_NDEFPAR-1))
    {
        idefl--; idefr--;
        idl1 = 0; idl2 = 0; idr1 = 0; idr2 = 0;
    }
    else{
        printf("defpar outside the scope of the grid. idefl = %d, idefr = %d\n",idefl,idefr);
        exit(1);
    }
    
    //printf("INTERPOLATION INDEX: %d, %d, %d\n", ind_a, ind_mu0, idefl);

    idefl1 = idefl+idl1;
    idefl2 = idefl+1+idl2;
    idefr1 = idefr+idr1;
    idefr2 = idefr+1+idr2;

    ifac_dp1 = (unscaled_defpar_left - fitstab->defpar[ind_a][idefl1]) / (fitstab->defpar[ind_a][idefl2] - fitstab->defpar[ind_a][idefl1]);
    ifac_dp2 = (unscaled_defpar_right - fitstab->defpar[ind_a+1][idefr1]) / (fitstab->defpar[ind_a+1][idefr2] - fitstab->defpar[ind_a+1][idefr1]);
    
    assert(ifac_dp1 >= 0.0);
    assert(ifac_dp2 >= 0.0);
    
    if(defpar < 0){
        x1l = (-1.0)*fitstab->defpar[ind_a][idefl2]/unscaled_defpar_left_lim;
        x2l = (-1.0)*fitstab->defpar[ind_a][idefl1]/unscaled_defpar_left_lim;
        x1r = (-1.0)*fitstab->defpar[ind_a+1][idefr2]/unscaled_defpar_right_lim;
        x2r = (-1.0)*fitstab->defpar[ind_a+1][idefr1]/unscaled_defpar_right_lim;
    }
    else{
        x1l = (1.0)*fitstab->defpar[ind_a][idefl2]/unscaled_defpar_left_lim;
        x2l = (1.0)*fitstab->defpar[ind_a][idefl1]/unscaled_defpar_left_lim;
        x1r = (1.0)*fitstab->defpar[ind_a+1][idefr2]/unscaled_defpar_right_lim;
        x2r = (1.0)*fitstab->defpar[ind_a+1][idefr1]/unscaled_defpar_right_lim;
    }
    
    x1 = x1l + ifac_a*(x1r - x1l);
    x2 = x2l + ifac_a*(x2r - x2l);
    
    //printf("x2l = %f\t x2r = %f\t x1l = %f\t x1r = %f\nifac_dp1 = %f\t ifac_dp2 = %f\n", x2l, x2r, x1l, x1r, ifac_dp1, ifac_dp2);
    
    interpol_elements(fitstab, ind_a, ind_mu0, ifac_a, ifac_mu0, defpar, x1, x2, x1l, x1r, x2l, x2r);
    
    // printf(" point 9\n");
    
    for (ii=0; ii < fitstab->n_r; ii++){

    	tdp1 = interp_lin_1d(ifac_dp1,Ftab->arr[ind_a][idefl1]->F[ii],Ftab->arr[ind_a][idefl2]->F[ii]);
    	tdp2 = interp_lin_1d(ifac_dp2,Ftab->arr[ind_a+1][idefr1]->F[ii],Ftab->arr[ind_a+1][idefr2]->F[ii]);

    	F[ii] = interp_lin_1d(ifac_a, tdp1, tdp2);

		robs[ii] = interp_lin_3d_nk_float(
			fitstab->arr[ind_a][idefl+idl1][ind_mu0]->r[ii],
			fitstab->arr[ind_a+1][idefr+idr1][ind_mu0]->r[ii],
			fitstab->arr[ind_a][idefl+idl1][ind_mu0+1]->r[ii],
			fitstab->arr[ind_a][idefl+1+idl2][ind_mu0]->r[ii],
			fitstab->arr[ind_a][idefl+1+idl2][ind_mu0+1]->r[ii],
			fitstab->arr[ind_a+1][idefr+1+idr2][ind_mu0+1]->r[ii],
			fitstab->arr[ind_a+1][idefr+idr1][ind_mu0+1]->r[ii],
			fitstab->arr[ind_a+1][idefr+1+idr2][ind_mu0]->r[ii],
			d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7], d[8]);

		// printf(" %d r=%f, F=%e\n", ii, robs[ii], F[ii]);

	    gmin[ii] = interp_lin_3d_nk_float(
			fitstab->arr[ind_a][idefl+idl1][ind_mu0]->gmin[ii],
			fitstab->arr[ind_a+1][idefr+idr1][ind_mu0]->gmin[ii],
			fitstab->arr[ind_a][idefl+idl1][ind_mu0+1]->gmin[ii],
			fitstab->arr[ind_a][idefl+1+idl2][ind_mu0]->gmin[ii],
			fitstab->arr[ind_a][idefl+1+idl2][ind_mu0+1]->gmin[ii],
			fitstab->arr[ind_a+1][idefr+1+idr2][ind_mu0+1]->gmin[ii],
			fitstab->arr[ind_a+1][idefr+idr1][ind_mu0+1]->gmin[ii],
			fitstab->arr[ind_a+1][idefr+1+idr2][ind_mu0]->gmin[ii],
			d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7], d[8]);

		gmax[ii] = interp_lin_3d_nk_float(
			fitstab->arr[ind_a][idefl+idl1][ind_mu0]->gmax[ii],
			fitstab->arr[ind_a+1][idefr+idr1][ind_mu0]->gmax[ii],
			fitstab->arr[ind_a][idefl+idl1][ind_mu0+1]->gmax[ii],
			fitstab->arr[ind_a][idefl+1+idl2][ind_mu0]->gmax[ii],
			fitstab->arr[ind_a][idefl+1+idl2][ind_mu0+1]->gmax[ii],
			fitstab->arr[ind_a+1][idefr+1+idr2][ind_mu0+1]->gmax[ii],
			fitstab->arr[ind_a+1][idefr+idr1][ind_mu0+1]->gmax[ii],
			fitstab->arr[ind_a+1][idefr+1+idr2][ind_mu0]->gmax[ii],
			d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7], d[8]);




		    for (jj = 0; jj < fitstab->n_g; jj++) {
				
				trff1[ii][jj] = interp_lin_3d_nk_float(
					fitstab->arr[ind_a][idefl+idl1][ind_mu0]->trff1[ii][jj],
					fitstab->arr[ind_a+1][idefr+idr1][ind_mu0]->trff1[ii][jj],
					fitstab->arr[ind_a][idefl+idl1][ind_mu0+1]->trff1[ii][jj],
					fitstab->arr[ind_a][idefl+1+idl2][ind_mu0]->trff1[ii][jj],
					fitstab->arr[ind_a][idefl+1+idl2][ind_mu0+1]->trff1[ii][jj],
					fitstab->arr[ind_a+1][idefr+1+idr2][ind_mu0+1]->trff1[ii][jj],
					fitstab->arr[ind_a+1][idefr+idr1][ind_mu0+1]->trff1[ii][jj],
					fitstab->arr[ind_a+1][idefr+1+idr2][ind_mu0]->trff1[ii][jj],
					d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7], d[8]);
				
				trff2[ii][jj] = interp_lin_3d_nk_float(
					fitstab->arr[ind_a][idefl+idl1][ind_mu0]->trff2[ii][jj],
					fitstab->arr[ind_a+1][idefr+idr1][ind_mu0]->trff2[ii][jj],
					fitstab->arr[ind_a][idefl+idl1][ind_mu0+1]->trff2[ii][jj],
					fitstab->arr[ind_a][idefl+1+idl2][ind_mu0]->trff2[ii][jj],
					fitstab->arr[ind_a][idefl+1+idl2][ind_mu0+1]->trff2[ii][jj],
					fitstab->arr[ind_a+1][idefr+1+idr2][ind_mu0+1]->trff2[ii][jj],
					fitstab->arr[ind_a+1][idefr+idr1][ind_mu0+1]->trff2[ii][jj],
					fitstab->arr[ind_a+1][idefr+1+idr2][ind_mu0]->trff2[ii][jj],
					d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7], d[8]);
			
		        double c1, c2, c3, c4, c5, c6, c7, c8;
				
				c1 = fitstab->arr[ind_a][idefl+idl1][ind_mu0]->cosne1[ii][jj];
				c2 = fitstab->arr[ind_a+1][idefr+idr1][ind_mu0]->cosne1[ii][jj];
				c3 = fitstab->arr[ind_a][idefl+idl1][ind_mu0+1]->cosne1[ii][jj];
				c4 = fitstab->arr[ind_a][idefl+1+idl2][ind_mu0]->cosne1[ii][jj];
				c5 = fitstab->arr[ind_a][idefl+1+idl2][ind_mu0+1]->cosne1[ii][jj];
				c6 = fitstab->arr[ind_a+1][idefr+1+idr2][ind_mu0+1]->cosne1[ii][jj];
				c7 = fitstab->arr[ind_a+1][idefr+idr1][ind_mu0+1]->cosne1[ii][jj];
				c8 = fitstab->arr[ind_a+1][idefr+1+idr2][ind_mu0]->cosne1[ii][jj];		
				
				
				cosne1[ii][jj] = interp_lin_3d_nk_float(c1, c2, c3, c4, c5, c6, c7, c8, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7], d[8]);

				if (cosne1[ii][jj] > 1.) {
					cosne1[ii][jj] = 1.0 - pow(10,-6);
				}

				c1 = fitstab->arr[ind_a][idefl+idl1][ind_mu0]->cosne2[ii][jj];
				c2 = fitstab->arr[ind_a+1][idefr+idr1][ind_mu0]->cosne2[ii][jj];
				c3 = fitstab->arr[ind_a][idefl+idl1][ind_mu0+1]->cosne2[ii][jj];
				c4 = fitstab->arr[ind_a][idefl+1+idl2][ind_mu0]->cosne2[ii][jj];
				c5 = fitstab->arr[ind_a][idefl+1+idl2][ind_mu0+1]->cosne2[ii][jj];
				c6 = fitstab->arr[ind_a+1][idefr+1+idr2][ind_mu0+1]->cosne2[ii][jj];
				c7 = fitstab->arr[ind_a+1][idefr+idr1][ind_mu0+1]->cosne2[ii][jj];
				c8 = fitstab->arr[ind_a+1][idefr+1+idr2][ind_mu0]->cosne2[ii][jj];

				cosne2[ii][jj] = interp_lin_3d_nk_float(c1, c2, c3, c4, c5, c6, c7, c8, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7], d[8]);
				if (cosne2[ii][jj] > 1.) {
					cosne2[ii][jj] = 1.0 - pow(10,-6);
				}

		    }

    }
    
    return;
}






// static void interpol_fitsTable(relSysPar** sysPar_inp, double a, double defpar, double mu0, double rin, double rout,
// 		 int* status){
void get_interpolated_fitsTable(double a, double defpar, double mu0, int* status){
    
    // load tables
    // printf(" check_caching_defpartype=%d\n", check_caching_defpartype());
	if (tr_table==NULL || check_caching_defpartype()){
		read_tr_table(get_filename(),&tr_table,status); 			// nonkerr mods, uncomment before compiling
		CHECK_STATUS_VOID(*status);
		// printf(" point 3\n");
		read_Ftable(get_filenameF(),&F_table,status); 			// nonkerr mods, uncomment before compiling
		CHECK_STATUS_VOID(*status);
	}
	// printf(" cached_defpartype=%d\n", cached_defpartype);
	// set_cached_defpartype();
	fitsTable* tab = tr_table;
	assert(tab!=NULL);

	FTable* Ftab = F_table;
	assert(Ftab!=NULL);

	// printf(" point 7\n");
    interpol_a_defpar_mu0(a, defpar, mu0, tab, Ftab);

	return;
}






// nonkerr modifications
double interp_lin_3d_nk_float(double v1, double v2, double v3, double v4, double v5, double v6, double v7, double v8,
	double d, double d1, double d2, double d3, double d4, double d5, double d6, double d7, double d8){
	
	double res = ((v1 * d1 + v2 * d2 + v3 * d3 + v4 * d4 + d5 * v5 + d6 * v6 + d7 * v7 + d8 * v8) / d);
	return res;
}



/**  FLOAT search for value "val" in array "arr" (sorted ASCENDING!) with length n and
 	 return bin k for which arr[k]<=val<arr[k+1] **/
int binary_search_float(float* arr,int n,float val){

	int klo=0;
	int khi=n-1;
	int k=-1;
	while ( (khi-klo) > 1 ){
		k=(khi+klo)/2;
		if(arr[k]>val){
			khi=k;
		} else {
			klo=k;
		}
	}
	return klo;
}












void init_par_nkbb(kerrbbParam** krb_param, const double* inp_par, const int n_parameter, int* status){

	kerrbbParam* param = new_kerrbbParam(status);
	CHECK_STATUS_VOID(*status);

	assert(n_parameter == NUM_PARAM_NKBB);

	param->eta        = inp_par[0];
	param->a          = inp_par[1];
	param->incl       = cos(inp_par[2]*M_PI/180);
	param->Mbh        = inp_par[3];
	param->Mdd        = inp_par[4];
	param->Dbh        = inp_par[5];
	param->h          = inp_par[6];
	param->rflag      = inp_par[7];
	param->lflag      = inp_par[8];
	param->defpartype = inp_par[9];
	param->defpar     = inp_par[10];
	set_defpar(param->defpar);
	set_defpartype(param->defpartype);


	CHECK_STATUS_VOID(*status);

	*krb_param  = param;

	return;
}



// nonkerr modifications
/** KERRBB NONKERR MODEL FUNCTION**/
void tdnkbb(const double* ener0, const int n_ener0, double* photar, const double* parameter, const int n_parameter, int* status){

	kerrbbParam* krb_param = NULL;

	init_par_nkbb(&krb_param,parameter,n_parameter,status);
	CHECK_STATUS_VOID(*status);

	double* ener = (double*) ener0;
	nkbb_kernel(ener, photar, n_ener0, krb_param, status);
	CHECK_STATUS_VOID(*status);

	free_kerrbbParam(krb_param);

}


/***************************************************/
/**** XSPEC MODEL WRAPPER **************************/
/***************************************************/
/** XSPEC KERRBB NONKERR MODEL FUNCTION **/

void lmodnkbb(const double* ener0, const int n_ener0, const double* parameter, int ifl, double* photar, double* photer, const char* init){

	const int n_parameter = 11;
	int status = EXIT_SUCCESS;
	tdnkbb(ener0, n_ener0, photar, parameter, n_parameter, &status);

	if (status!=EXIT_SUCCESS)
	TR_ERROR(" Evaluating NKBB model failed",&status);
}




/* get a new kerrbb parameter structure and initialize it */
kerrbbParam* new_kerrbbParam(int* status){
	kerrbbParam* param = (kerrbbParam*) malloc(sizeof(kerrbbParam));
	if (param==NULL){
		TR_ERROR("memory allocation failed",status);
		return NULL;
	}
	
	param->eta        = PARAM_DEFAULT;
	param->a          = PARAM_DEFAULT;
	param->incl       = PARAM_DEFAULT;
	param->Mbh        = PARAM_DEFAULT;
	param->Mdd        = PARAM_DEFAULT;
	param->Dbh        = PARAM_DEFAULT;
	param->h          = PARAM_DEFAULT;
	param->rflag      = PARAM_DEFAULT;
	param->lflag      = PARAM_DEFAULT;
	param->defpartype = 1;
	param->defpar     = PARAM_DEFAULT;
	
	set_defpar(param->defpar);
	set_defpartype(param->defpartype);

	return param;
}

/* free kerrbb parameter */
void free_kerrbbParam(kerrbbParam* param){
	free(param);
}




// nkbb_kernel(ener0, photar, n_ener0, krb_param, status);
void nkbb_kernel(double* ener0, double* photar, int n_ener0, kerrbbParam* krb_param, int* status ){
    
    int ii, jj, kk;
	int truncated_index;
	get_interpolated_fitsTable(krb_param->a, krb_param->defpar, krb_param->incl, status);
	// printf(" cached_defpartype=%d\n", cached_defpartype);
	set_cached_defpartype();

	// set_cached_rel_param(param , &cached_rel_param_nk, status);
	// CHECK_STATUS_RET(*status,NULL);
    
    double g_star[FITSTABLE_NG] = {0.002, 0.00368523, 0.0056385, 0.00792911, 0.01065242, 0.01394334, 0.01799906, 0.02311944, 0.02978329, 0.03880341, 0.0516735,  0.07144532, 0.10533413, 0.15, 0.20384615, 0.25769231, 0.31153846, 0.36538462, 0.41923077, 0.47307692, 0.52692308, 0.58076923, 0.63461538, 0.68846154, 0.74230769, 0.79615385, 0.85, 0.89466587, 0.92855468, 0.9483265,  0.96119659, 0.97021671, 0.97688056, 0.98200094, 0.98605666, 0.98934758, 0.99207089, 0.9943615, 0.99631477, 0.998};
    
    double temparray[FITSTABLE_NG][2];
    double A1, A2, Upsilon1, Upsilon2, g_factor, temp1, temp2, temp3, sum, func_0, func_max;
    double mean_ener0;
    double t[FITSTABLE_NR], integrate_func[FITSTABLE_NR], integrate_value;
	double eta, delta, frac1, frac2;

	eta = krb_param->eta;
	delta = eta*robs[FITSTABLE_NR - 1];
	
    /*for (ii=0; ii < FITSTABLE_NR; ii++) {
	printf("%d r=%f, gmin=%f, gmax=%f, tr1[0]=%f, tr1[39]=%f, tr2[0]=%f, tr2[39]=%f, c1[0]=%f, c1[39]=%f, c2[0]=%f, c2[39]=%f, F[]=%e\n",
			ii,robs[ii],gmin[ii], gmax[ii], trff1[ii][0],trff1[ii][FITSTABLE_NG - 1], trff2[ii][0],trff2[ii][FITSTABLE_NG - 1], cosne1[ii][0], cosne1[ii][FITSTABLE_NG - 1], cosne2[ii][0], cosne2[ii][FITSTABLE_NG - 1], F[ii]);
    }*/

	//Setup for a truncated disk
	if(eta > 0.0)
	{
		for(ii = FITSTABLE_NR - 2; ii > 0; ii--)
		{
			if(robs[ii] - robs[FITSTABLE_NR - 1] > delta)
			{
				truncated_index = ii;
				//printf("eta-ii=%d\t%d\n", truncated_index + 1, truncated_index);
				frac1 = (robs[FITSTABLE_NR - 1] + delta - robs[truncated_index + 1])/(robs[truncated_index] - robs[truncated_index + 1]);
				frac2 = (robs[truncated_index] - robs[FITSTABLE_NR - 1] - delta)/(robs[truncated_index] - robs[truncated_index + 1]);

				F[truncated_index + 1] = F[truncated_index + 1]*frac2 + F[truncated_index]*frac1;
				gmin[truncated_index + 1] = gmin[truncated_index + 1]*frac2 + gmin[truncated_index]*frac1;
				gmax[truncated_index + 1] = gmax[truncated_index + 1]*frac2 + gmax[truncated_index]*frac1;

				for (jj = 0; jj < FITSTABLE_NG; jj++) {
					trff1[truncated_index + 1][jj] = trff1[truncated_index + 1][jj]*frac2 + trff1[truncated_index][jj]*frac1;
					trff2[truncated_index + 1][jj] = trff2[truncated_index + 1][jj]*frac2 + trff2[truncated_index][jj]*frac1;
					cosne1[truncated_index + 1][jj] = cosne1[truncated_index + 1][jj]*frac2 + cosne1[truncated_index][jj]*frac1;
					cosne2[truncated_index + 1][jj] = cosne2[truncated_index + 1][jj]*frac2 + cosne2[truncated_index][jj]*frac1;
				}

				robs[truncated_index + 1] = robs[FITSTABLE_NR - 1] + delta;
				break;
			}
		}

		if(truncated_index < FITSTABLE_NR - 2)
		{
			for(ii = FITSTABLE_NR - 1; ii > truncated_index + 1; ii--)
			{
				for(jj = 0; jj < FITSTABLE_NG; jj++)
				{
					trff1[ii][jj] = 0;
					trff2[ii][jj] = 0;
				}
			}
		}
	}

	for (ii = 0; ii < FITSTABLE_NR; ii++) {
        if(F[ii] < 0.0) {
            F[ii] = - F[ii];
            printf("Check Dimless F!\n");
        }
        
        for (jj = 0; jj < FITSTABLE_NG; jj++) {
            if(trff1[ii][jj] < 0.0)
                trff1[ii][jj] = 0.0;
            if(trff2[ii][jj] < 0.0)
                trff2[ii][jj] = 0.0;
        }
    }

	/*for (ii=0; ii < FITSTABLE_NR; ii++) {
	printf("%d r=%f, gmin=%f, gmax=%f, tr1[0]=%f, tr1[39]=%f, tr2[0]=%f, tr2[39]=%f, c1[0]=%f, c1[39]=%f, c2[0]=%f, c2[39]=%f, F[]=%e\n",
			ii,robs[ii],gmin[ii], gmax[ii], trff1[ii][0],trff1[ii][FITSTABLE_NG - 1], trff2[ii][0],trff2[ii][FITSTABLE_NG - 1], cosne1[ii][0], cosne1[ii][FITSTABLE_NG - 1], cosne2[ii][0], cosne2[ii][FITSTABLE_NG - 1], F[ii]);
    }*/
	
    
    A1 = 0.07205*pow(krb_param->Mbh, 2.0)/(pow(krb_param->h, 4.0)*pow(krb_param->Dbh, 2.0));
    A2 = 0.1331*sqrt(krb_param->Mbh)/((krb_param->h)*pow(krb_param->Mdd, 0.25));
    
    for (ii = 0; ii < FITSTABLE_NR; ii++)
        t[ii] = pow(F[ii], 0.25);
    
    if(krb_param->lflag > 0.0)
    {
        for(kk = 0; kk < n_ener0; kk++)
        {
            mean_ener0 = 0.5*(ener0[kk] + ener0[kk + 1]);
            
            for (ii = 0; ii < FITSTABLE_NR; ii++) {
                for(jj = 0; jj < FITSTABLE_NG; jj++)
                {
                    g_factor = g_star[jj]*(gmax[ii] - gmin[ii]) + gmin[ii];
                    temp1 = exp((A2*mean_ener0)/(g_factor*t[ii])) - 1.0;
                    
                    Upsilon1 = 0.5 + 0.75*cosne1[ii][jj];
                    temp2 = Upsilon1/temp1;
                    temp3 = temp2*Pi*robs[ii]*trff1[ii][jj]/g_factor;
                    temparray[jj][0] = temp3/sqrt(g_star[jj]*(1.0 - g_star[jj]));
                    
                    Upsilon2 = 0.5 + 0.75*cosne2[ii][jj];
                    temp2 = Upsilon2/temp1;
                    temp3 = temp2*Pi*robs[ii]*trff2[ii][jj]/g_factor;
                    temparray[jj][1] = temp3/sqrt(g_star[jj]*(1.0 - g_star[jj]));
                }
                
                temp1 = gmin[ii]*(exp(A2*mean_ener0/(gmin[ii]*t[ii])) - 1.0);
                temp2 = gmax[ii]*(exp(A2*mean_ener0/(gmax[ii]*t[ii])) - 1.0);
                func_0 = Pi*robs[ii]*0.25*(trff1[ii][0] + trff2[ii][0])*(1.0 + 0.75*cosne1[ii][0] + 0.75*cosne2[ii][0])/temp1;
                func_max = Pi*robs[ii]*0.25*(trff1[ii][FITSTABLE_NG - 1] + trff2[ii][FITSTABLE_NG - 1])*(1.0 + 0.75*cosne1[ii][FITSTABLE_NG - 1] + 0.75*cosne2[ii][FITSTABLE_NG - 1])/temp2;
                
                sum = 0.0894726*(func_0 + func_max)*2.0;
                
                for(jj = 0; jj < FITSTABLE_NG - 1; jj++)
                {
                    sum += 0.5*(temparray[jj][0] + temparray[jj + 1][0])*(g_star[jj + 1] - g_star[jj]) + 0.5*(temparray[jj][1] + temparray[jj + 1][1])*(g_star[jj + 1] - g_star[jj]);
                }
                
                integrate_func[ii] = sum;
            }
            
            integrate_value = 0.0;
            for (ii = 0; ii < FITSTABLE_NR - 1; ii++)
            {
                integrate_value += (integrate_func[ii] + integrate_func[ii + 1])*0.5*(robs[ii] - robs[ii + 1]);
            }
            
            photar[kk] = A1*pow(mean_ener0, 2.0)*integrate_value*(ener0[kk + 1] - ener0[kk]);
        }
    }
    else
    {
        for(kk = 0; kk < n_ener0; kk++)
        {
            mean_ener0 = 0.5*(ener0[kk] + ener0[kk + 1]);
            
            for (ii = 0; ii < FITSTABLE_NR; ii++) {
                for(jj = 0; jj < FITSTABLE_NG; jj++)
                {
                    g_factor = g_star[jj]*(gmax[ii] - gmin[ii]) + gmin[ii];
                    temp1 = exp((A2*mean_ener0)/(g_factor*t[ii])) - 1.0;
                    
                    temp2 = 1.0/temp1;
                    temp3 = temp2*Pi*robs[ii]*trff1[ii][jj]/g_factor;
                    temparray[jj][0] = temp3/sqrt(g_star[jj]*(1.0 - g_star[jj]));
                    
                    //temp2 = 1.0/temp1;
                    temp3 = temp2*Pi*robs[ii]*trff2[ii][jj]/g_factor;
                    temparray[jj][1] = temp3/sqrt(g_star[jj]*(1.0 - g_star[jj]));
                }
                
                func_0 = Pi*robs[ii]*0.5*(trff1[ii][0] + trff2[ii][0])/(gmin[ii]*(exp(A2*mean_ener0/(gmin[ii]*t[ii])) - 1.0));
                func_max = Pi*robs[ii]*0.5*(trff1[ii][FITSTABLE_NG - 1] + trff2[ii][FITSTABLE_NG - 1])/(gmax[ii]*(exp(A2*mean_ener0/(gmax[ii]*t[ii])) - 1.0));
                
                sum = 0.0894726*(func_0 + func_max)*2.0;
                
                for(jj = 0; jj < FITSTABLE_NG - 1; jj++)
                {
                    sum += 0.5*(temparray[jj][0] + temparray[jj + 1][0])*(g_star[jj + 1] - g_star[jj]) + 0.5*(temparray[jj][1] + temparray[jj + 1][1])*(g_star[jj + 1] - g_star[jj]);
                }
                
                integrate_func[ii] = sum;
            }
            
            
            integrate_value = 0.0;
            for (ii = 0; ii < FITSTABLE_NR - 1; ii++)
            {
                integrate_value += (integrate_func[ii] + integrate_func[ii + 1])*0.5*(robs[ii] - robs[ii + 1]);
            }
     
            photar[kk] = A1*pow(mean_ener0, 2.0)*integrate_value*(ener0[kk + 1] - ener0[kk]);
        }
    }
    
	return;
    
}
