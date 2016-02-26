#define ROTATE(am,i,j,k,l) g=am[i][j];h=am[k][l];am[i][j]=g-s*(h+g*tau);am[k][l]=h+s*(g-h*tau);
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

#define LMAX     1
#define NMAX     31
#define H        0.0010   // Accuracy parameter for MT// timestep for TTL
#define DTOUT    150.0
#define DTWRITE  150.0
#define BMAX     20000010
#define PI       3.1415927
#define ORDER    0      // Choice of order for integration, allowed 0 (2nd order), 1 (4th order), 2 (6th order)
#define RMAX	 100.0 //50 original	//Escape radius for removal of particle from simulation
#define ACCURACY 1.0E-6//1.0E-6


/* Symplectic integrator coefficients. */
#define BETA              1.2599210498948732

static double si_c[][8] = {{.5, .5, .0, .0, .0, .0, .0, .0},
    {1.0/2.0/(2.0-BETA), (1.0-BETA)/2.0/(2.0-BETA), (1.0-BETA)/2.0/(2.0-BETA), 1.0/2.0/(2.0-BETA), .0, .0, .0, .0},
    {.392256805238780, .510043411918458,-.471053385409758, .0687531682525198, .0687531682525198, -.471053385409758,.510043411918458, .392256805238780}};

static double si_d[][8] = {{1.0, .0, .0, .0, .0, .0, .0, .0},
    {1.0/(2.0-BETA), -BETA/(2.0-BETA), 1.0/(2.0-BETA), .0, .0, .0, .0, .0},
				{.784513610477560, .235573213359357, -1.17767998417887, 1.31518632068391, -1.17767998417887, .235573213359357,.784513610477560, 0.}};

static int si_coeffs[] = {2, 4, 8};
double tim[2] = {0.0, 0.0};

char *defv[] = {
    "in=sersic_10K",
    "tdis=0",         // Jacobi radius (0=off, 1=Van Hoerner, 2=Kuepper method)
    "gasexp=0",       // Gas expulsion (0=off, 1=Kroupa & Boily 2002)
    "gcevolu=0",      // GC size evolution (0=ff, 1=on)
    "scalemass=0.0",  // Scalefactor mass
    "scalesize=0.0",  // Scalefactor size
    "xcrit=0.5",      // Half mass/Jacobi radius
    "dfmethod=0",     // Dynamical friction method (0= off, 1= standard Chandrasekhar, 2=triaxial generalization)
    "rcap=0.0",       // Capture radius BH
    "tend=100.0",     // Simulation end time
    "msmbh=0.00",     // SMBH mass fraction
    "restart=N",      // New Run or Restart ?
    NULL,
};

double *t0,*mass,*x,*y,*z,*vx,*vy,*vz,*phif;
double *mgc,*rh0,*rh;
double *laddx,*laddy,*laddz,*tboxo,*tboxn,*ts;
int *box,*boy,*boz,*boxorb,*grid_cell_nr;


//GC Mass loss by Relaxation
double *rperi,*rapo,*va,*vp,*mgct0,*tgc,*tdisso;


//model variables
int rkselect = 1; //0: MT; 1: TTL Integrator
double scalemass,scalesize,xcrit,tinit,tscale;
double abar[LMAX][NMAX],narr[LMAX][LMAX];
double clmfull[LMAX][LMAX][NMAX],slmfull[LMAX][LMAX][NMAX];
double rcap,tend,msmbh,msmbh0;
double tout,twrite,etot0,ecorr;
double mcom,rxcom,rycom,rzcom,vxcom,vycom,vzcom,rcom,vcom;
int nbody=0,ncap=0;
int tdis,gcevolu,gasexp,ngcdis=0,nescape=0,dfmethod;
int nbox=0;
void read_nbody(),get_coeff0(),mt_drv(),read_restart(),get_host();
void get_coeff(),get_pot(),getforce_scf();
int check_smbh_capture();
int check_gc_disruption();
double get_densiti();
double cputim(),get_poti();
double get_phi_smbh();
void get_force_smbh(double *xi, double *a);
void com_correction();

//TTL Integrator
void get_ttl();
void  get_forcettl();

//Dynamical friction
void get_force_df();
void get_sigma();
double get_integral();
double qgaus();
void get_grid();
//int grid_cell_nr[BMAX];
#define GRIDCELLS 5
#define GRIDDEPTH 15  //15 orig
double grid_coordinates[98*GRIDDEPTH+1][4];
int grid_cell_occupation[98*GRIDDEPTH+1];
double sigma_1[98*GRIDDEPTH+1];
double sigma_2[98*GRIDDEPTH+1];
double sigma_3[98*GRIDDEPTH+1];
double alpha_euler[98*GRIDDEPTH+1];
double beta_euler[98*GRIDDEPTH+1];
double gamma_euler[98*GRIDDEPTH+1];
double e1[98*GRIDDEPTH+1][3];
double e2[98*GRIDDEPTH+1][3];
double e3[98*GRIDDEPTH+1][3];
double logkernel();

//NUMERICAL RECIPES
#define NR_END 1
#define FREE_ARG char*
void jacobi();
void eigsrt();
double factrl(),gammln();
void nrerror(char error_text[]);
#define BIGNUMBER 1.0e30


float *vector(long nl, long nh) {
    /* allocate a float vector with subscript range v[nl..nh] */
    
    float *v;
    
    v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
    if (!v) nrerror("allocation failure in vector()");
    
    return v-nl+NR_END;
    
}

void free_vector(float *v, long nl, long nh) {
    /* free a float vector allocated with vector() */
    
    free((FREE_ARG) (v+nl-NR_END));
    
}

void nrerror(char error_text[]) {
    /* Numerical Recipes standard error handler */
    
    fprintf(stderr,"Numerical Recipes run-time error...\n");
    fprintf(stderr,"%s\n",error_text);
    fprintf(stderr,"...now exiting to system...\n");
    
    exit(1);
    
}

