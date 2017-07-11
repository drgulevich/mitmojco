/* 
	example-2.c: Voltage-biased SIS junction 
*/
//==============================================//
//---------- by Dmitry R. Gulevich -------------//
//------ d.r.gulevich@metalab.ifmo.ru ----------//
//--- ITMO University, St Petersburg, Russia ---//
//==============================================//
/* MODEL DESCRIPTION:
AC driven voltage-biased SIS junction.
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
#include "mitmojco.h" // MiTMoJCo header file
#include "opt_filter.h" // optimum filtration

/* Model Parameters */
#define OMP_NUM_THREADS 1 //  number of OpenMP threads (use 1 for small Josephson contact)
#define SETTLING_TIME 200. // time interval upon which voltage starts being recorded
#define TMAX 500. // integration time
#define DT 0.01 // time step
#define N_OPT_FILTER 5 // optimum filtration level (1 for direct summation)
#define AMP_FILE "amplitudes/BCS42_008.fit" // tunnel current amplitudes file
#define A_SUPP 1.0 // pair current suppression
#define KGAP 3.3 // normalized gap frequency (omega_g/omega_J)

void sis_vbias(double Vac, double omega, double V_start, double V_finish, double V_step);

/**************************************************************************************************
 NAME:	    main
 INPUT:     command line arguments
 FUNCTION:  process command line arguments and execute simulation
**************************************************************************************************/
int main (int argc, char* argv[])
	{

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

	printf("#=================================================\n");
	printf("#---- Example 2: Voltage-biased SIS junction -----\n");
	printf("#=================================================\n");
	
	double Vac, omega, Vdc_start, Vdc_finish, Vdc_step;
	
	if(argc-optind==5) {
		Vac = atof(argv[optind]);
		omega = atof(argv[optind+1]);
		Vdc_start = atof(argv[optind+2]);
		Vdc_finish = atof(argv[optind+3]);
		Vdc_step = atof(argv[optind+4]);
	} else {
		printf("# Incorrect number of arguments.\n");
		printf("# Please, supply 5 arguments in the following order:\n");
		printf("# 1. Vac\n");
		printf("# 2. omega\n");
		printf("# 3. Vdc_start\n");
		printf("# 4. Vdc_finish\n");
		printf("# 5. Vdc_step\n");
		return -1;
	}
	
	printf("# Vac: %.3f\n", Vac);
	printf("# omega: %.3f\n", omega);
	printf("# Vdc_start: %.3f\n", Vdc_start);
	printf("# Vdc_finish: %.3f\n", Vdc_finish);
	printf("# Vdc_step: %.3f\n", Vdc_step);

	omp_set_num_threads(OMP_NUM_THREADS);
	printf("# OMP threads: %d\n",omp_get_max_threads());

	sis_vbias(Vac, omega, Vdc_start, Vdc_finish, Vdc_step);
	
	return 0;
}


/**************************************************************************************************
 NAME:	    sis_vbias
 FUNCTION:  start simulation of voltage-biased SIS junction
 INPUTS: 
 			Vac: amplitude of the drive
 			omega: angular frequency of the drive
			Vdc_start: initial value of DC voltage
			Vdc_finish: final value of DC voltage
			Vdc_step: step size for DC voltage
 COMMENTS:
			Both upward (Vdc_start<Vdc_finish) and downward (Vdc_start>Vdc_finish) voltage
			sweeps are possible. 
**************************************************************************************************/
void sis_vbias(double Vac, double omega, double Vdc_start, double Vdc_finish, double Vdc_step) {

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
		sis_tunnel_current->jbar: pointer to the tunnel current array
	*/

	/* Check that no errors ocurred */
	if( sis_tunnel_current->error )
		return;

	/* Create filter object (OptFilterType pointer) for current filtration.
		The constructor accepts one argument - the optimum filtration level
		(parameter n in Ref. A. A. Odintsov, V. K. Semenov and A. B. Zorin, 
		IEEE Trans. Magn. 23, 763 (1987)). n=1 corresponds to the arithmetic mean.
	*/
	OptFilterType *current_filter = opt_filter_create(N_OPT_FILTER);; 

	double phi0, Vdc, t;
	
	/* Set initial conditions. Josephson steps are absent for phi0=0. and present otherwise. */
	phi0=0.;

	/* Treating the upward (Vdc_start<Vdc_finish) and downward (Vdc_start>Vdc_finish) sweeps */
	Vdc_step=fabs(Vdc_step);
	Vdc_step = (Vdc_start<= Vdc_finish)?Vdc_step:(-Vdc_step);   	

	printf("# phi0: %f\n",phi0);
	printf("#----------- Starting the calculation ------------\n");
	printf("# Column 1: Vdc\n");
	printf("# Column 2: gamma\n");

	clock_t cticks0 = clock(); // cpu timing
	double wtime0 = omp_get_wtime(); // wall clock timing

	for(Vdc=Vdc_start; ((Vdc < Vdc_finish + 0.5*Vdc_step) && (Vdc_step>0)) ||
		((Vdc > Vdc_finish + 0.5*Vdc_step) && (Vdc_step<0)); Vdc+=Vdc_step) {

		mitmojco_init( sis_tunnel_current );  // Initialize the tunnel current object

		opt_filter_init(current_filter); // Initialize the filter object
		
		for(t=DT; t<TMAX; t+=DT) {
		
			/* In this example, SIS is driven by a single harmonics.
				In general, an arbitrary signal can stay here */
			double phidot = 2.*(Vdc+Vac*cos(omega*t));
			phi = phi0 + 2.*Vdc*t + 2.*Vac*sin(omega*t)/omega;

			mitmojco_update( sis_tunnel_current ); // Update the tunnel current

			/* Total tunnel current including the normal resistance term */
			double gamma = sis_tunnel_current->jbar[0] + sis_tunnel_current->alphaN * phidot;

			/* Make a record to the filter object */
			if(t > SETTLING_TIME)
				opt_filter_update(current_filter, gamma);			

		} 
		
		printf("%.3f %.6f\n",Vdc, opt_filter_result(current_filter) ); 	// Output
		
	}

	clock_t cticks1 = clock();
	double wtime1 = omp_get_wtime();

	printf("#-------------------------------------------------\n");
	printf ("# CPU time:        %f\n", (float) (cticks1 - cticks0)/CLOCKS_PER_SEC);
	printf ("# Wall clock time: %f\n", wtime1-wtime0);

	/* Clear memory allocated for the objects */
	mitmojco_free(sis_tunnel_current);
	opt_filter_free(current_filter);
}

