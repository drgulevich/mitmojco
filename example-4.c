/* 
	example-4.c:
	Fluxon in an Annular Josephson junction
*/
//==============================================//
//---------- by Dmitry R. Gulevich -------------//
//------ d.r.gulevich@metalab.ifmo.ru ----------//
//--- ITMO University, St Petersburg, Russia ---//
//==============================================//
/* MODEL DESCRIPTION:
     physical nodes: 1 to Nnodes
     shadow nodes: 0 and Nnodes+1
*/

#include <stdio.h>
#include <stdlib.h> // malloc
#include <complex.h> // complex numbers
#include <math.h> // sqrt, sin, cos
#include <string.h> // memcpy
#include <stdbool.h> // boolean datatype
#include <time.h> // timing
#include <omp.h> // OpenMP
#include <ctype.h> // character handling for getopt routine
#include <unistd.h> // getopt routine
#include "mitmojco.h" // MiTMoJCo header file

/* Model Parameters */
#define OMP_NUM_THREADS 4 //  number of OpenMP threads
#define SETTLING_TIME 200. // time interval upon which voltage starts being recorded
#define TMAX 1000. // integration time
#define AMP_FILE "amplitudes/BCS42_008.fit" // tunnel current amplitudes file
#define A_SUPP 0.7 // pair current suppression
#define KGAP 3.3 // normalized gap frequency (omega_g/omega_J)
#define WINDING 1
#define L0 40.0
#define DTREL 0.5
#define BETA 0.0
#define DTOUT 2.
#define Nnodes 1000

void fluxon(double *phi);
void ajj_set(double *phi, double *phi_old);
int find_fluxon(double *phi);
void ajj(double gamma_start, double gamma_finish, double gamma_step);

const double dx=L0/Nnodes;
const double xmin=-0.5*L0-0.5*(L0/Nnodes); // position of the left shadow node
const double dt=DTREL*L0/Nnodes;


/**************************************************************************************************
 NAME:	    main
 FUNCTION:  process command line arguments and execute simulation
 INPUTS:     command line arguments
**************************************************************************************************/
int main (int argc, char* argv[]) {

	// Check if parameters are consistent
	if(TMAX<SETTLING_TIME)
		printf("Error: TMAX should be larger than SETTLING_TIME\n");

	int c;
	opterr = 0;
	while ((c = getopt (argc, argv, "f:")) != -1)
	    switch (c)
	    {
	    case '?':
				if (isprint (optopt))
		          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
		        else
		          fprintf (stderr,
		                   "Unknown option character `\\x%x'.\n",
		                   optopt);
		        return 1;
	    default:
		        abort ();
	    }

	double gamma_start, gamma_finish, gamma_step = 0.01;
	
	printf("#=========================================================\n");
	printf("#---- Example 4: Fluxon in Annular Josephson Junction ----\n");
	printf("#=========================================================\n");

	if(argc-optind==1) {
		gamma_start = atof(argv[optind]);
		printf("# gamma: %.4f\n", gamma_start);
		gamma_finish=gamma_start;
	} else if(argc-optind==3) {
		gamma_start = atof(argv[optind]);
		gamma_finish = atof(argv[optind+1]);
		gamma_step = atof(argv[optind+2]);
		printf("# gamma: (%.4f, %.4f, %.4f)\n", gamma_start,gamma_finish,gamma_step);
	} else {
		printf("# Incorrect number of arguments.\n");
		printf("# Please, supply 1 or 3 arguments in the following order:\n");
		printf("# 1. gamma_start\n");
		printf("# 2. gamma_finish\n");
		printf("# 3. gamma_step\n");
		return -1;
	}

	omp_set_num_threads(OMP_NUM_THREADS);
	printf("# OMP threads: %d\n",omp_get_max_threads());
	
	printf("# L0: %.4f\n", L0);
//	printf("# Nnodes: %d\n", Nnodes);
	printf("# BETA: %.4f\n", BETA);
	printf("# SETTLING_TIME: %.2f\n", SETTLING_TIME);
	printf("# TMAX: %.2f\n", TMAX);
	
	ajj(gamma_start, gamma_finish, gamma_step);
	
return 0;
}


/**************************************************************************************************
 NAME:	    ajj
 FUNCTION:  start simulation of AJJ
 INPUTS: 
			gamma_start: initial value of the bias current
			gamma_finish: final value of the bias current
			gamma_step: step size for the bias current
**************************************************************************************************/
void ajj(double gamma_start, double gamma_finish, double gamma_step) {

	double *phi = malloc( (Nnodes+2) * sizeof(double) ); 	// superconducting phase difference array
	int skipinds[2] = {0,Nnodes+1}; // indices of shadow nodes

	/* Create tunnel current object (TunnelCurrentType pointer) by calling the constructor with arguments:		
		1. Tunnel current amplitudes file.
		2. Pair current suppression parameter (1 if no suppression).
		3. Normalized gap frequency (omega_g/omega_J).
		4. Integration time step. 
		5. Total number of nodes (size of array phi).
		6. Pointer to array phi.
		7. Number of shadow nodes to be skipped.
		8. Pointer to the array of indices of the shadow nodes.
	*/

	TunnelCurrentType *ajj_tunnel_current = mitmojco_create( AMP_FILE, A_SUPP, KGAP, dt, Nnodes+2, phi, 2, skipinds);

	/* Object member variables are accessed via the arrow operator '->' :
		ajj_tunnel_current->error: error flag (false if no errors occurred)
		ajj_tunnel_current->alphaN: damping due to a pure normal resistance
		ajj_tunnel_current->jbar: array of the tunnel currents
	*/

	/* Check that no errors ocurred */
	if( ajj_tunnel_current->error )
		return;

	int i, count;
	double gamma, t;
	double shift;
	double *phi_old = malloc( (Nnodes+2) * sizeof(double) );
	double *phi_new = malloc( (Nnodes+2) * sizeof(double) );
	double *gradphi = malloc( Nnodes * sizeof(double) );
	double x_initial, x_final, phi_fluxon_initial, phi_fluxon_final;
	int i_fluxon, ncycles;
	bool started;

	double dt2=dt*dt;
	double r2=dt2/(dx*dx);
	double q = BETA*dt/(dx*dx);
	double factor=1./(1 + 0.5 * ajj_tunnel_current->alphaN * dt);

	//--- Display info ---//
	printf("# dx: %.4f\n", dx);
	printf("# dt: %.4f\n", dt);

	//--- Initialization ---//
	ajj_set( phi, phi_old );
	mitmojco_init( ajj_tunnel_current );  // Initialize the tunnel current object

	//--- Start the gamma loop ---//
	gamma_step=fabs(gamma_step);
	gamma_step = (gamma_start<= gamma_finish)?gamma_step:(-gamma_step);
	printf("#-------------------------------------------------\n");
	printf("# 1 gamma\n");
	printf("# 2 velocity\n");

	clock_t cticks0 = clock(); // cpu timing
	double wtime0 = omp_get_wtime(); // wall clock timing

	for(gamma=gamma_start; ((gamma < gamma_finish + 0.5*gamma_step) && (gamma_step>0)) ||
			((gamma > gamma_finish + 0.5*gamma_step) && (gamma_step<0)); gamma+=gamma_step) {

		for(t=0.,count=0,started=false; t<TMAX; t+=dt, count++) {
	        if(t > SETTLING_TIME) {
				if(!started) {
					i_fluxon = find_fluxon(phi);
					x_initial = xmin + i_fluxon*dx;
					phi_fluxon_initial = phi[i_fluxon];
					started = true;
				}
			}

	      	mitmojco_update( ajj_tunnel_current ); // Update the tunnel current

			// Calculate phi_new: internal nodes
		    for(i=1;i<=Nnodes;i++)
		        phi_new[i] = factor*( (r2+q)*phi[i-1] + 2.*(1-r2-q)*phi[i] + (r2+q)*phi[i+1] +
		        		(-1+0.5* ajj_tunnel_current->alphaN*dt+2.*q)*phi_old[i] - q*phi_old[i-1] - q*phi_old[i+1] 
						+ dt2*(gamma - ajj_tunnel_current->jbar[i]) );

			phi_new[0] = phi_new[Nnodes] + 2.*M_PI*WINDING;
			phi_new[Nnodes+1] = phi_new[1] - 2.*M_PI*WINDING;
			
		   	// Update variables
	        memcpy(phi_old,phi,(Nnodes+2)*sizeof(double));
	        memcpy(phi,phi_new,(Nnodes+2)*sizeof(double));

		} // end of t loop

		i_fluxon = find_fluxon(phi);
		x_final = xmin + i_fluxon*dx;
		phi_fluxon_final = phi[i_fluxon];
	
		ncycles = round((phi_fluxon_final - phi_fluxon_initial)/(2.*M_PI));
		double distance = ncycles*L0 + x_final-x_initial;

		printf("%f %f\n", gamma, distance/(t-SETTLING_TIME));

		// Shift phase by multiple of 4*pi on all nodes
		shift = 4.*M_PI * floor(phi[0] /(4.*M_PI)); // 4.*M_PI because of sin(0.5*phi)
		for(i=0;i<Nnodes+2;i++) {
			phi[i] = phi[i] - shift;
			phi_old[i] = phi_old[i] - shift;
		}

	} // end of gamma loop

	clock_t cticks1 = clock();
	double wtime1 = omp_get_wtime();
	
	printf("#-------------------------------------------------\n");
	printf ("# CPU time:        %f\n", (float) (cticks1 - cticks0)/CLOCKS_PER_SEC);
	printf ("# Wall clock time: %f\n", wtime1-wtime0);

	free(phi);
	free(phi_old);
	free(gradphi);
	free(phi_new);

	mitmojco_free(ajj_tunnel_current); 	// Clear memory allocated for the object
}


/**************************************************************************************************
 NAME:	    fluxon
 INPUTS:    
            phi: pointer to the superconducting phase difference array
 FUNCTION:  
            set phi to fluxon profile
**************************************************************************************************/
void fluxon(double *phi) {
    int i;
    double x; 
    for(i=1;i<=Nnodes;i++) {
        x=xmin+i*dx;
        phi[i] = -WINDING*4.0*atan(exp(x));
    }
	phi[0] = phi[1];
	phi[Nnodes+1] = phi[Nnodes];
}

/**************************************************************************************************
 NAME:	    ajj_set
 INPUTS:    
            phi: pointer to the superconducting phase difference array
            phi_old: pointer to the superconducting phase difference array on the previous step            
 FUNCTION:  
            set initial conditions
**************************************************************************************************/
void ajj_set(double *phi, double *phi_old) {
	fluxon(phi);
	memcpy(phi_old,phi,(Nnodes+2)*sizeof(double));
}

/**************************************************************************************************
 NAME:	    find_fluxon
 INPUTS:    
            phi: pointer to the superconducting phase difference array
 FUNCTION:  
            find fluxon location
 RETURNS:
            node index of the center of fluxon
**************************************************************************************************/
int find_fluxon(double *phi) {
	double shift = 2.*M_PI * round(phi[1] /(2.*M_PI));
	int i;
	for(i=1;i<=Nnodes;i++)
		if(phi[i] < shift - M_PI || phi[i] > shift + M_PI) {
			break;
		}
	return i;
}

