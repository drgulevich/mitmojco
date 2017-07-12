/* 
	example-5.c: Flux Flow Oscillator
*/
//==============================================//
//---------- by Dmitry R. Gulevich -------------//
//------ d.r.gulevich@metalab.ifmo.ru ----------//
//--- ITMO University, St Petersburg, Russia ---//
//==============================================//
/* MODEL DESCRIPTION:
Numerical calculation of current-voltage characteristics of FFO using a quasi-1d model, Eq.(15)-(16)
in arXiv:1704.03045, without load. Geometrical parameters are defined as follows:

                           L0                          
         <----------------------------------->        
            _______________________________           
    /|\    /                               \          
     |    /                                 \         
  W0 |   |                                   | Wend   
     |    \                                 /         
    \|/    \_______________________________/          
         <->                               <->        
         S0                                S0         

At hext>0 and gamma>0 fluxons enter through the injection end (at x = -L0/2) and leave through the radiation 
end (at x = L0/2). The quasi-one dimensional model (Eq.15 in arXiv:1704.03045) is discretized in space in 
(Nnodes+2) nodes: 1 to Nnodes are the physics nodes through which currents are assumed to flow, and two 
shadow nodes 0 and Nnodes+1 for treatment of the boundary conditions. The injection end (x = -L0/2) is 
half-step between the nodes 0 and 1 and the radiation end (x = L0/2) is half-step between the nodes Nnodes 
and Nnodes+1.*/

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

#define LAMBDA_J 5.5 // in microns
const double L0 = 400./LAMBDA_J; // Length
const double W0 = 16./LAMBDA_J; // Width
const double S0 = 40./LAMBDA_J; // Length of shapered region
const double Wend = 1.0/LAMBDA_J; // Width of the end
const double DX = 0.5/LAMBDA_J; // Spatial discretization

/* Model Parameters */
#define OMP_NUM_THREADS 4 //  number of OpenMP threads (use 1 for small Josephson contact)
#define SETTLING_TIME 200. // time interval upon which voltage starts being recorded
#define TMAX 1000. // integration time
#define N_OPT_FILTER 5 // optimum filtration level (1 for direct summation)
#define AMP_FILE "amplitudes/BCS42_008.fit" // tunnel current amplitudes file
#define A_SUPP 0.7 // pair current suppression
#define KGAP 3.3 // normalized gap frequency (omega_g/omega_J)
#define DTREL 0.5 // ratio (dt/dx)
#define BETA 0.02 // surface damping

int Nnodes;
double dx;
double xmin;
double *width;
double *dlogw;
double Area;

void ffo_configure();
void ffo_set_IC(double *phi, double *phi_old, double hext);
void set_Gammaeff(double gamma, double *Gammaeff);
void ffo_explicit(double hext, double gamma_start, double gamma_finish, double gamma_step);
void ffo_info();


/**************************************************************************************************
 NAME:      main
 INPUTS:    command line arguments
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

	printf("#===============================================\n");
	printf("#------- Example 5: Flux Flow Oscillator -------\n");
	printf("#===============================================\n");


	double hext, gamma_start, gamma_finish, gamma_step;
	
	if(argc-optind==2) {
		hext = atof(argv[optind]);
		gamma_start = atof(argv[optind+1]);
		printf("# Magnetic field (hext): %.4f\n", hext);
		printf("# Bias current (gamma): %.4f\n", gamma_start);
		gamma_finish=gamma_start;
		gamma_step = 0.01;
	} else if(argc-optind==4) {
		hext = atof(argv[optind]);
		gamma_start = atof(argv[optind+1]);
		gamma_finish = atof(argv[optind+2]);
		gamma_step = atof(argv[optind+3]);
		printf("# Magnetic field (hext): %.4f\n", hext);
		printf("# Bias current (gamma): (%.4f, %.4f, %.4f)\n", gamma_start,gamma_finish,gamma_step);
	} else {
		printf("# Incorrect number of arguments.\n");
		printf("# Please, supply 2 or 4 arguments in the following order:\n");
		printf("# 1. Magnetic field (hext)\n");
		printf("# 2. Value of bias current (gamma_start)\n");
		printf("# or\n");
		printf("# 1. Magnetic field (hext)\n");
		printf("# 2. Starting value of bias current (gamma_start)\n");
		printf("# 3. Final value of bias current (gamma_finish)\n");
		printf("# 4. Bias current step (gamma_step)\n");
		return -1;
	}
	
	// Set number of OpenMP threads
	omp_set_num_threads(OMP_NUM_THREADS);

	// Set FFO
	ffo_configure();

	// Display info
	ffo_info();

	// Run simulation
	ffo_explicit(hext, gamma_start, gamma_finish, gamma_step);
	
	return 0;
}


/**************************************************************************************************
 NAME:      ffo_configure
 FUNCTION:  initialize global variables which are constant throughout the simulation
**************************************************************************************************/
void ffo_configure() {

	// Initialize global variables
	Nnodes = round(L0/DX); // number of nodes
	dx=L0/Nnodes; // spatial discretization
	xmin=-0.5*L0-0.5*dx; // position of the left shadow node

	int i;
	double absWder = (W0-Wend)/S0;
	width  = malloc( (Nnodes+2) * sizeof(double) ); // width W(x)
	dlogw  = malloc( (Nnodes+2) * sizeof(double) ); // logarithmics derivative W'(x)/W(x)

	// Loop over physical nodes (1 to Nnodes)
	for(i=1;i<=Nnodes;i++) 
	{ 
		double x=xmin + i*dx;
		if(x < -0.5*L0 + S0) {
			width[i] = Wend + absWder*(x+0.5*L0);
			dlogw[i] = absWder/width[i];
		}
		else if(x > 0.5*L0 - S0) {
			width[i] = Wend - absWder*(x-0.5*L0);
			dlogw[i] = - absWder/width[i];
		}
		else {
			width[i] = W0;
			dlogw[i] = 0.;
		}
	}

	// Calculate FFO area
	Area =0.;
	for(i=1;i<=Nnodes;i++)
		Area += width[i];
	Area *= dx;
}


/**************************************************************************************************
 NAME:      ffo_info
 FUNCTION:  print FFO model parameters
**************************************************************************************************/
void ffo_info() {
	printf("# OpenMP threads: %d\n",omp_get_max_threads());
	printf("# Tapered ends FFO geometry:\n");
	printf("# 	L0: %.3f micron\n", L0*LAMBDA_J);
	printf("# 	W0: %.3f micron\n", W0*LAMBDA_J);
	printf("# 	S0: %.3f micron\n", S0*LAMBDA_J);
	printf("# 	Wend: %.3f micron\n", Wend*LAMBDA_J);
	printf("# 	Josephson penetration (LAMBDA_J): %.2f micron\n", LAMBDA_J);
	printf("# Surface damping (BETA): %.4f\n", BETA);
	printf("# Settling time (SETTLING_TIME): %.2f\n", SETTLING_TIME);
	printf("# Integration time (TMAX): %.2f\n", TMAX);
	printf("# FFO area: %.2f LAMBDA_J^2\n",Area);
}


/**************************************************************************************************
 NAME:      ffo_set_IC
 FUNCTION:  set initial state of superconducting phase difference.
 INPUTS:
            phi: superconducting phase difference 
            phi_old: superconducting phase difference for the previous step
            hext: external magnetic field
 COMMENTS:
            The initial state may be an approximate solution or a random guess. A realistic guess 
            is desirable to minimize the transient dynamics	and guarantee establishing of the flux 
            flow regime. In this implementation a linear function satisfying the boundary 
            conditions at the FFO ends is used to initialize phi and phi_old.
**************************************************************************************************/
void ffo_set_IC(double *phi, double *phi_old, double hext) {
	int i;
	for(i=0;i<Nnodes+2;i++) {
		double x=xmin + i*dx;
		phi[i] = -hext*x;
		phi_old[i] = phi[i];
	}
}


/**************************************************************************************************
 NAME:      set_Gammaeff
 FUNCTION:  set the effective bias current at each physical node
 INPUTS: 
            gamma: average bias current defined by Eq.(18) in arXiv:1704.03045.
            Gammaeff: pointer to the effective bias current array defined by Eq.(17) in arXiv:1704.03045.
 COMMENTS:
            The model assumes a homogeneos feed of the bias current to the FFO. The homogeneous feed
            results in an inhomogeneous effective bias current for a FFO with variable width.
**************************************************************************************************/
void set_Gammaeff(double gamma, double *Gammaeff) {
	int i;
	for(i=1;i<=Nnodes;i++) {
		double h2gamma_i = gamma*Area/L0;
		Gammaeff[i] = h2gamma_i/width[i];
	}
}


/**************************************************************************************************
 NAME:      ffo_explicit
 FUNCTION:  start simulation of FFO dynamics using 2nd order central differences
 INPUTS:
            hext: magnetic field
            gamma_start: initial value of the bias current
            gamma_finish: final value of the bias current
            gamma_step: step size for the bias current
 COMMENTS:
            Both upward (gamma_start<gamma_finish) and downward (gamma_start>gamma_finish) current 
            sweeps are possible. 
**************************************************************************************************/
void ffo_explicit(double hext, double gamma_start, double gamma_finish, double gamma_step) {

	double dt=DTREL*dx; // time step

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

	TunnelCurrentType *ffo_tunnel_current = mitmojco_create( AMP_FILE, A_SUPP, KGAP, dt, Nnodes+2, phi, 2, skipinds);

	/* Object member variables are accessed via the arrow operator '->' :
		ffo_tunnel_current->error: error flag (false if no errors occurred)
		ffo_tunnel_current->alphaN: damping due to a pure normal resistance
		ffo_tunnel_current->jbar: array of the tunnel currents
	*/

	/* Check that no errors ocurred */
	if( ffo_tunnel_current->error )
		return;

	//--- Display info ---//
	printf("# Number of physical nodes (Nnodes): %d\n", Nnodes);
	printf("# dx: %.5f\n", dx);
	printf("# dt: %.5f\n", dt);

	int i;
	double gamma, t;
	double shift, Vdc;
	double dt2=dt*dt;
	double r2=dt2/(dx*dx);
	double q = BETA*dt/(dx*dx);
	double factor=1./(1 + 0.5 * ffo_tunnel_current->alphaN * dt);
	double coeff0 = 0.5*BETA*dt/dx;
	double coeff1 = coeff0 + 0.5*dt2/dx;
	double *phi_old = malloc( (Nnodes+2) * sizeof(double) );
	double *phi_new = malloc( (Nnodes+2) * sizeof(double) );
	double *jbar = malloc( (Nnodes+2) * sizeof(double) );
	double *Gammaeff = malloc( (Nnodes+2) * sizeof(double) );

	// Upper bound for the flux flow voltage
	double Vdc_max = 0.8*KGAP;

	/* Create filter object (OptFilterType pointer) for voltage filtration.
		The constructor accepts one argument - the optimum filtration level
		(parameter n in Ref. A. A. Odintsov, V. K. Semenov and A. B. Zorin, 
		IEEE Trans. Magn. 23, 763 (1987)). n=1 corresponds to the arithmetic mean.
	*/
	OptFilterType *voltage_filter = opt_filter_create(N_OPT_FILTER);

	gamma_step=fabs(gamma_step);
	gamma_step = (gamma_start<= gamma_finish)?gamma_step:(-gamma_step);
	printf("#----------- Starting the calculation ------------\n");
	printf("# Column 1: gamma\n");
	printf("# Column 2: Vdc\n");

	double wtime0 = omp_get_wtime();
	clock_t cticks0 = clock();

	// Set initial conditions
	ffo_set_IC(phi, phi_old, hext);
	
	// Initialize MiTMoJCo
	mitmojco_init( ffo_tunnel_current );  // Initialize the tunnel current object
	
	for(gamma=gamma_start; (((gamma < gamma_finish + 0.5*gamma_step) && (gamma_step>0)) ||
			((gamma > gamma_finish + 0.5*gamma_step) && (gamma_step<0))) && Vdc <= Vdc_max; gamma+=gamma_step) {
	
		set_Gammaeff(gamma,Gammaeff);
	
		opt_filter_init(voltage_filter); // Initialize the filter object
	
		for(t=0.; t<TMAX; t+=dt) {
		
			mitmojco_update( ffo_tunnel_current ); // Update the tunnel current
	
			// Loop over physical nodes
			for(i=1;i<=Nnodes;i++)
				phi_new[i] = factor*( (r2+q-dlogw[i]*coeff1)*phi[i-1] + 2.*(1-r2-q)*phi[i] + (r2+q+dlogw[i]*coeff1)*phi[i+1] +
						(-1+0.5*ffo_tunnel_current->alphaN*dt+2.*q)*phi_old[i] + (-q+dlogw[i]*coeff0)*phi_old[i-1] + (-q-dlogw[i]*coeff0)*phi_old[i+1] 
						+ dt2*(Gammaeff[i] - ffo_tunnel_current->jbar[i] + dlogw[i]*hext) );
	
			// Update two shadow nodes
			phi_new[0] = phi_new[1] + dx*hext;
			phi_new[Nnodes+1] = phi_new[Nnodes] - dx*hext;
	
			// Make a record to the filter object: 1st physical node at the injection end is chosen
			if(t > SETTLING_TIME)
				opt_filter_update(voltage_filter, phi_new[1]-phi_old[1]);
	
		   	// Update variables
			memcpy(phi_old,phi,(Nnodes+2)*sizeof(double));
			memcpy(phi,phi_new,(Nnodes+2)*sizeof(double));
		}
	

		/* DC voltage in units hbar*omega_J/e.	Factor 0.25 comes from the Josephson relation 
			Vdc = 0.5*dphi/dt, and 2nd order discretization for the derivative dphi/dt */
		Vdc=0.25*opt_filter_result(voltage_filter)/dt;		
	
		// Output
		printf("%.4f %f\n",gamma, Vdc);
	
		// Shift phase to avoid precision error
		shift = 4.*M_PI * floor(phi[0] /(4.*M_PI));
		for(i=0;i<Nnodes+2;i++) {
			phi[i] = phi[i] - shift;
			phi_old[i] = phi_old[i] - shift;
		}
	
	} // end of gamma loop
	
	double wtime1 = omp_get_wtime();
	clock_t cticks1 = clock();
	
	printf("#-------------------------------------------------\n");
	printf ("# CPU time:        %f\n", (float) (cticks1 - cticks0)/CLOCKS_PER_SEC);
	printf ("# Wall clock time: %f\n", wtime1-wtime0);
	
	free(phi);
	free(phi_old);
	free(phi_new);
	free(jbar);
	free(Gammaeff);

	/* Clear memory allocated for the objects */
	mitmojco_free( ffo_tunnel_current );
	opt_filter_free(voltage_filter);
}


