/* 
	example-3.c:
	Sine-Gordon breather in long Josephson junction
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
#define TMAX 150. // integration time
#define AMP_FILE "amplitudes/BCS42_008.fit" // tunnel current amplitudes file
#define A_SUPP 0.7 // pair current suppression
#define KGAP 3.3 // normalized gap frequency (omega_g/omega_J)
#define L0 40.0
#define Nnodes 500
#define DTREL 0.5
#define BETA 0.0
//#define DTOUT 0.5
#define DTOUT 0.05

void breather(double *phi);

const double dx=L0/Nnodes;
const double xmin=-0.5*L0-0.5*(L0/Nnodes); // position of the left shadow node


/**************************************************************************************************
 NAME:	    main
 INPUT:     
 FUNCTION:  display info and simulation of long Josephson junction
**************************************************************************************************/
int main () {
	
	printf("#===========================================\n");
	printf("#----- Example 3: Sine-Gordon breather -----\n");
	printf("#===========================================\n");

	omp_set_num_threads(OMP_NUM_THREADS);
	printf("# OMP threads: %d\n",omp_get_max_threads());
	printf("# L0: %.4f\n", L0);
	printf("# Nnodes: %d\n", Nnodes);
	printf("# BETA: %.4f\n", BETA);
	printf("# TMAX: %.2f\n", TMAX);
	printf("# DTOUT: %.2f\n", DTOUT);

	double dt=DTREL*dx;	
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

	TunnelCurrentType *ljj_tunnel_current = mitmojco_create( AMP_FILE, A_SUPP, KGAP, dt, Nnodes+2, phi, 2, skipinds);

	/* Object member variables are accessed via the arrow operator '->' :
		ljj_tunnel_current->error: error flag (false if no errors occurred)
		ljj_tunnel_current->alphaN: damping due to a pure normal resistance
		ljj_tunnel_current->jbar: array of the tunnel currents
	*/

	/* Check that no errors ocurred */
	if( ljj_tunnel_current->error )
		return 1;

	int i, count;
	double t;
	double *phi_old = malloc( (Nnodes+2) * sizeof(double) );
	double *phi_new = malloc( (Nnodes+2) * sizeof(double) );

	FILE *fout_phi;
	fout_phi = fopen("breather.dat","w");

	double dt2=dt*dt;
	double r2=dt2/(dx*dx);
	double q = BETA*dt/(dx*dx);
	double factor=1./(1 + 0.5 * ljj_tunnel_current->alphaN * dt);

	int countout=(int)fmax(1.,DTOUT/dt);

	printf("# countout: %d\n", countout);

	//--- Display info ---//
	printf("# dx: %.4f\n", dx);
	printf("# dt: %.4f\n", dt);

	//--- Initialization ---//
	breather(phi);
	memcpy(phi_old,phi,(Nnodes+2)*sizeof(double));
	mitmojco_init( ljj_tunnel_current );  // Initialize the tunnel current object

	clock_t cticks0 = clock(); // cpu timing
	double wtime0 = omp_get_wtime(); // wall clock timing

	for(t=0.,count=0; t<TMAX; t+=dt, count++) {

		// Output phi into file
		if(count%countout==0) {
			for(i=1;i<=Nnodes;i++)
				fprintf(fout_phi,"%.4f ",phi[i]);
			fprintf(fout_phi,"\n");
		}

		mitmojco_update( ljj_tunnel_current ); // Update the tunnel current

		// Calculate phi_new: internal nodes
		for(i=1;i<=Nnodes;i++)
			phi_new[i] = factor*( (r2+q)*phi[i-1] + 2.*(1-r2-q)*phi[i] + (r2+q)*phi[i+1] +
					(-1+0.5* ljj_tunnel_current->alphaN*dt+2.*q)*phi_old[i] - q*phi_old[i-1] - q*phi_old[i+1] 
					+ dt2*( -ljj_tunnel_current->jbar[i]) );

		// Open boundary conditions
		phi_new[0] = phi_new[1];
		phi_new[Nnodes+1] = phi_new[Nnodes];
			
	   	// Update variables
		memcpy(phi_old,phi,(Nnodes+2)*sizeof(double));
		memcpy(phi,phi_new,(Nnodes+2)*sizeof(double));

	} // end of t loop

	clock_t cticks1 = clock();
	double wtime1 = omp_get_wtime();
	
	printf("#-------------------------------------------------\n");
	printf ("# CPU time:        %f\n", (float) (cticks1 - cticks0)/CLOCKS_PER_SEC);
	printf ("# Wall clock time: %f\n", wtime1-wtime0);

	free(phi);
	free(phi_old);
	free(phi_new);

	mitmojco_free(ljj_tunnel_current); 	// Clear memory allocated for the object

	fclose(fout_phi);

	return 0;
}


void breather(double *phi) {
	int i;
	double x, u=0.2, gfactor=1.0/sqrt(1.0 + u*u ); 
	for(i=1;i<=Nnodes;i++) {
		x=xmin+i*dx;
		phi[i] = 4.0*atan(cos( gfactor*u*0.0 )/(u*cosh( gfactor*x )));
	}
	phi[0] = phi[1];
	phi[Nnodes+1] = phi[Nnodes];
}
