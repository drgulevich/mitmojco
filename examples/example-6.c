/* 
	example-6.c: Gap resonance (see Ref. D. V. Averin, Gap resonance in 
	the classical dynamics of the current-biased Josephson tunnel junctions,
	arXiv:2111.07206)
*/
//==============================================//
//---------- by Dmitry R. Gulevich -------------//
//----------- drgulevich@gmail.com -------------//
//==============================================//
/* MODEL DESCRIPTION:
*/

#include <stdio.h>
#include <stdlib.h> // malloc
#include <complex.h> // complex numbers
#include <math.h> // sqrt, sin, cos
#include <string.h> // memcpy
#include <time.h> // timing
#include <omp.h> // OpenMP
#include <ctype.h> // character handling for getopt routine
#include <unistd.h> // getopt routine
#include <mitmojco/mitmojco.h> // MiTMoJCo main header file
#include <mitmojco/opt_filter.h> // MiTMoJCo optimum filtration

/* Model Parameters */
#define OMP_NUM_THREADS 1 //  number of OpenMP threads (use 1 for small Josephson contact)
#define SETTLING_TIME 200. // time interval upon which voltage starts being recorded
#define TMAX 1000. // integration time
#define DT 0.01 // time step
#define N_OPT_FILTER 5 // optimum filtration level (1 for direct summation)
#define AMP_FILE "../amplitudes/NbNb_4K2_008.fit" // tunnel current amplitudes file
#define A_SUPP 1.0 // pair current suppression
#define KGAP 0.5 // normalized gap frequency (omega_g/omega_J)

void sis_cbias(double gamma_start, double gamma_finish, double gamma_step);

/**************************************************************************************************
 NAME:      main
 INPUT:     command line arguments
 FUNCTION:  process command line arguments and execute simulation
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
					fprintf (stderr, "Unknown option character `\\x%x'.\n",	optopt);
				return 1;
		default:
				abort ();
		}

	printf("#==============================================\n");
	printf("#--- Example 6: Gap resonance ---\n");
	printf("#==============================================\n");

	double gamma_start, gamma_finish, gamma_step;
	
	if(argc-optind==1) {
		gamma_start = atof(argv[optind]);
		printf("# Bias current (gamma): %.4f\n", gamma_start);
		gamma_finish=gamma_start;
		gamma_step = 0.01;
	} else if(argc-optind==3) {
		gamma_start = atof(argv[optind]);
		gamma_finish = atof(argv[optind+1]);
		gamma_step = atof(argv[optind+2]);
		printf("# Bias current (gamma): (%.4f, %.4f, %.4f)\n", gamma_start,gamma_finish,gamma_step);
	} else {
		printf("# Incorrect number of arguments.\n");
		printf("# Please, supply 1 or 3 arguments in the following order:\n");
		printf("# 1. Value of bias current (gamma_start)\n");
		printf("# or\n");
		printf("# 1. Starting value of bias current (gamma_start)\n");
		printf("# 2. Final value of bias current (gamma_finish)\n");
		printf("# 3. Bias current step (gamma_step)\n");
		return -1;
	}

	omp_set_num_threads(OMP_NUM_THREADS);
	printf("# OMP threads: %d\n",omp_get_max_threads());
	printf("# Settling time (SETTLING_TIME): %.2f\n", SETTLING_TIME);
	printf("# Integration time (TMAX): %.2f\n", TMAX);
	printf("# Time step (DT): %.5f\n", DT);

	// Run simulation
	sis_cbias(gamma_start, gamma_finish, gamma_step);
	
	return 0;
}


/**************************************************************************************************
 NAME:      sis_cbias
 FUNCTION:  start simulation of SIS junction dynamics
 INPUTS:
            gamma_start: initial value of the bias current
            gamma_finish: final value of the bias current
            gamma_step: step size for the bias current
 COMMENTS:
            Both upward (gamma_start<gamma_finish) and downward (gamma_start>gamma_finish) current 
            sweeps are possible. 
**************************************************************************************************/
void sis_cbias(double gamma_start, double gamma_finish, double gamma_step) {

	double phi; // superconducting phase difference

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

	TunnelCurrentType *sis_tunnel_current = mitmojco_create( AMP_FILE, A_SUPP, KGAP, DT, 1, &phi, 0, NULL );

	/* Object member variables are accessed via the arrow operator '->' :
		sis_tunnel_current->error: error flag (false if no errors occurred)
		sis_tunnel_current->alphaN: damping due to a pure normal resistance
		sis_tunnel_current->jbar: array of the tunnel currents
	*/

	/* Check that no errors occurred */
	if( sis_tunnel_current->error )
		return;

	/* Create filter object (OptFilterType pointer) for voltage filtration.
		The constructor accepts one argument - the optimum filtration level
		(parameter n in Ref. A. A. Odintsov, V. K. Semenov and A. B. Zorin, 
		IEEE Trans. Magn. 23, 763 (1987)). n=1 corresponds to the arithmetic mean.
	*/
	OptFilterType *voltage_filter = opt_filter_create(N_OPT_FILTER);

	double gamma, t;
	double shift, Vdc;
	double phi_old, phi_new;
	double factor=1./(1+0.5 * sis_tunnel_current->alphaN * DT);


	/* Set the initial conditions */
	phi=0.;
	phi_old=0.;
	
	mitmojco_init( sis_tunnel_current );  // Initialize the tunnel current object

	/* Treating the upward (gamma_start<gamma_finish) and downward (gamma_start>gamma_finish) sweeps */
	gamma_step=fabs(gamma_step);
	gamma_step = (gamma_start<= gamma_finish)?gamma_step:(-gamma_step);

	printf("#----------- Starting the calculation ------------\n");
	printf("# Column 1: gamma\n");
	printf("# Column 2: Vdc\n");
	
	clock_t cticks0 = clock(); // cpu timing
	double wtime0 = omp_get_wtime(); // wall clock timing

	// Down
	for(gamma=gamma_start; ((gamma < gamma_finish + 0.5*gamma_step) && (gamma_step>0)) ||
			((gamma > gamma_finish + 0.5*gamma_step) && (gamma_step<0)); gamma+=gamma_step) {
	
		opt_filter_init(voltage_filter); // Initialize the filter object
	
		for(t=0.; t<TMAX; t+=DT) {
	
			mitmojco_update( sis_tunnel_current ); // Update the tunnel current

			phi_new = factor*( 2.*phi + ( -1 + 0.5 * sis_tunnel_current->alphaN * DT)*phi_old + DT*DT*(gamma - sis_tunnel_current->jbar[0]) );
	
			/* Make a record to the filter object */
			if(t > SETTLING_TIME)
				opt_filter_update(voltage_filter, phi_new-phi_old);
	
		   	/* Update variables */
			phi_old = phi;
			phi = phi_new;
		}
	
		/* DC voltage in units hbar*omega_J/e.	Factor 0.25 comes from the Josephson relation Vdc = 0.5*dphi/DT 
			and 2nd order discretization for th derivative,	dphi/dt = 0.5*(phi_new-phi_old)/dt */
		Vdc=0.25*opt_filter_result(voltage_filter)/DT;
	
		printf("%.4f %f\n",gamma, Vdc);  // Output
	
		/* Shift phase to avoid precision errors */
		shift = 4.*M_PI * floor(phi /(4.*M_PI));
		phi = phi - shift;
		phi_old = phi_old - shift;	
	}
	

	gamma_finish = gamma_start;
	gamma_start = gamma;
	gamma_step = -gamma_step;
	
	// Up again to find the hysteresis
	for(gamma=gamma_start; ((gamma < gamma_finish + 0.5*gamma_step) && (gamma_step>0)) ||
			((gamma > gamma_finish + 0.5*gamma_step) && (gamma_step<0)); gamma+=gamma_step) {
	
		opt_filter_init(voltage_filter); // Initialize the filter object
	
		for(t=0.; t<TMAX; t+=DT) {
	
			mitmojco_update( sis_tunnel_current ); // Update the tunnel current

			phi_new = factor*( 2.*phi + ( -1 + 0.5 * sis_tunnel_current->alphaN * DT)*phi_old + DT*DT*(gamma - sis_tunnel_current->jbar[0]) );
	
			/* Make a record to the filter object */
			if(t > SETTLING_TIME)
				opt_filter_update(voltage_filter, phi_new-phi_old);
	
		   	/* Update variables */
			phi_old = phi;
			phi = phi_new;
		}
	
		/* DC voltage in units hbar*omega_J/e.	Factor 0.25 comes from the Josephson relation Vdc = 0.5*dphi/DT 
			and 2nd order discretization for th derivative,	dphi/dt = 0.5*(phi_new-phi_old)/dt */
		Vdc=0.25*opt_filter_result(voltage_filter)/DT;
	
		printf("%.4f %f\n",gamma, Vdc);  // Output
	
		/* Shift phase to avoid precision errors */
		shift = 4.*M_PI * floor(phi /(4.*M_PI));
		phi = phi - shift;
		phi_old = phi_old - shift;	
	}	
	
	

	clock_t cticks1 = clock();
	double wtime1 = omp_get_wtime();
	
	printf("#-------------------------------------------------\n");
	printf ("# CPU time:        %f\n", (float) (cticks1 - cticks0)/CLOCKS_PER_SEC);
	printf ("# Wall clock time: %f\n", wtime1-wtime0);

	/* Clear memory allocated for the objects */
	mitmojco_free(sis_tunnel_current);
	opt_filter_free(voltage_filter);
}


