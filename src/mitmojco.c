/* 
	MiTMoJCo
*/
//==============================================//
//---------- by Dmitry R. Gulevich -------------//
//------ d.r.gulevich@metalab.ifmo.ru ----------//
//--- ITMO University, St Petersburg, Russia ---//
//==============================================//
#include <stdio.h>
#include <stdlib.h> // malloc
#include <complex.h> // complex numbers
#include <math.h> // sqrt, sin, cos
#include <omp.h>
#include "mitmojco/mitmojco.h"

typedef struct {

	// interface to public struct
	TunnelCurrentType PubInterface;

	// private
	int Nnodes;
	int Nskip;
	int *skipinds;
	int Nexps;
	double kgap_over_Rejptilde0;
	double *sin05phi_old;
	double *cos05phi_old;
	double complex *C;
	double complex *D;
	double complex *cexp_z;
	double complex *alpha0;
	double complex *alpha1;
	double complex *F;
	double complex *G;

	double **phi_ptrs;
	double **jbar_ptrs;

} TunnelCurrentType_Private;

/* private methods */
int mitmojco_nfloats_in_file(const char *filename);
void mitmojco_loadamps(const char *filename, double complex *p, double complex *A, double complex *B);


/**************************************************************************************************
 NAME:      mitmojco_nfloats_in_file
 FUNCTION:  count number of floats in a file
 INPUTS:    filename: pointer to file name
 RETURNS:   number of floats
 COMMENTS:  private method (used internally by other MiTMoJCo methods)
**************************************************************************************************/
int mitmojco_nfloats_in_file(const char *filename) {
	int i=0;
	double value;
	FILE* f;
	if ((f = fopen(filename, "r")) != NULL) {
		while(fscanf(f,"%lf", &value) != EOF)
			i++;
		fclose(f);
		return i;
	}
	return -1;
}


/**************************************************************************************************
 NAME:      mitmojco_loadamps
 FUNCTION:  load tunnel current amplitudes from file (internal use only)
 INPUTS:    filename: pointer to tunnel current amplitude file name
			p, A, B: pointers to the fitting parameters
 COMMENTS:  private method (used internally by other MiTMoJCo methods)
**************************************************************************************************/
void mitmojco_loadamps(const char *filename, double complex *p, double complex *A, double complex *B) {
	double rep, imp, reA, imA, reB, imB;
	int i=0;
	FILE *f = fopen(filename, "r");
   	while(fscanf(f,"%lf", &rep ) != EOF && fscanf(f,"%lf", &imp ) != EOF && 
   		fscanf(f,"%lf", &reA ) != EOF && fscanf(f,"%lf", &imA ) != EOF &&
   		fscanf(f,"%lf", &reB ) != EOF && fscanf(f,"%lf", &imB ) != EOF) {
   			p[i] = rep + I*imp;
			A[i] = reA + I*imA;
			B[i] = reB + I*imB;
			i++;
   	}
	printf("# Amplitudes loaded from file \"%s\"\n",filename);
	fclose(f);
}


/**************************************************************************************************
 NAME:      mitmojco_create
 FUNCTION:  constructor of tunnel current amplitude object
 INPUTS:    filename: pointer to tunnel current amplitude file name
			a_supp: pair current suppression parameter
			kgap: normalized gap frequency
			dt: integration time step
			Ntotal: size of the array phi
			phi: pointer to the superconducting difference
			Nskip: number of shadow nodes to skip
			skipinds: pointer to the array of integers to skip
 RETURNS:   tunnel current amplitude object (pointer to TunnelCurrentType struct)
**************************************************************************************************/
TunnelCurrentType* mitmojco_create(const char *filename, double a_supp, double kgap, double dt, int Ntotal, double *phi, int Nskip, int *skipinds) {

	TunnelCurrentType_Private *s = malloc(sizeof(TunnelCurrentType_Private));

	printf("#------------------- MiTMoJCo --------------------\n");
	int Nvalues = mitmojco_nfloats_in_file(filename);
	if(Nvalues < 0) {
		printf("# Error: file %s not found\n",filename);
		s->PubInterface.error = true;
	} else if(Nvalues == 0 || Nvalues % 6 != 0) {
		printf("# Error: incorrect format of file %s\n",filename);
		s->PubInterface.error = true;
	}
	else {
		// Number of fitting exponents
		s->Nexps = Nvalues/6;

		// Declare fit parameters
		double complex *p = malloc( s->Nexps * sizeof(double complex) );
		double complex *A = malloc( s->Nexps * sizeof(double complex) );
		double complex *B = malloc( s->Nexps * sizeof(double complex) );

		// Load tunnel current amplitudes
		mitmojco_loadamps(filename, p, A, B);

		s->phi_ptrs = malloc( Ntotal * sizeof(double*) ); // array of pointers to phi at physical nodes
		s->jbar_ptrs = malloc( Ntotal * sizeof(double*) ); // array of pointers to jbar at physical nodes
		s->PubInterface.jbar = malloc( Ntotal * sizeof(double) );
		int i,j;
		int count = 0;
		bool accept;
		for(i=0;i<Ntotal;i++) {

			accept = true;

			// Skip shadow nodes
			for(j=0;j<Nskip;j++)
				if(i==skipinds[j])
					accept = false;

			// Accept physical nodes
			if(accept) {
				s->phi_ptrs[count] = &phi[i];
				s->jbar_ptrs[count] = &(s->PubInterface.jbar[i]);
				count++;	
			}
		}
		s->Nnodes = count; // Number of physical nodes
		s->PubInterface.error = false;
		s->PubInterface.self = s; // Keep record of self pointer
		s->sin05phi_old = malloc( s->Nnodes * sizeof(double) );
		s->cos05phi_old = malloc( s->Nnodes * sizeof(double) );
		s->C = malloc( s->Nexps * sizeof(double complex) );
		s->D = malloc( s->Nexps * sizeof(double complex) );
		s->cexp_z = malloc( s->Nexps * sizeof(double complex) );
		s->alpha0 = malloc( s->Nexps * sizeof(double complex) );
		s->alpha1 = malloc( s->Nexps * sizeof(double complex) );
		s->F = malloc( s->Nnodes * s->Nexps * sizeof(double complex) );
		s->G = malloc( s->Nnodes * s->Nexps * sizeof(double complex) );

		int n;
		for(n=0;n<(s->Nexps);n++) {
			double complex z=p[n]*kgap*dt;
			double complex ez=cexp(z);
			s->cexp_z[n] = ez;
			s->alpha0[n] = -ez + (ez - 1.)/z;
			s->alpha1[n] = 1. + (1. - ez)/z;
			s->C[n]=(a_supp*A[n]+B[n])/(-kgap*p[n]);
			s->D[n]=(a_supp*A[n]-B[n])/(-kgap*p[n]);
		}

		double Rejptilde0=0.;
		for(n=0;n<(s->Nexps);n++)
			Rejptilde0 -= creal(A[n]/p[n]);
		Rejptilde0 *= a_supp;
	
		s->PubInterface.Rejptilde0 = Rejptilde0;
		s->kgap_over_Rejptilde0 = kgap/Rejptilde0;
		s->PubInterface.alphaN = 1./(2.*kgap*Rejptilde0);

		printf("# MiTMoJCo started with parameters:\n");
		printf("# 	Normalized gap frequency (kgap): %.3f\n", kgap);
		printf("# 	Pair current suppression (a_supp): %.2f\n", a_supp);
		printf("# 	Rejp~(0): %.5f\n", s->PubInterface.Rejptilde0);
		printf("# 	alphaN: %.5f\n", s->PubInterface.alphaN);
		printf("# 	Number of fitting exponents (Nexps): %d\n", s->Nexps);

		free(p);
		free(A);
		free(B);
	}

	return &(s->PubInterface);
}


/**************************************************************************************************
 NAME:      mitmojco_init
 FUNCTION:  initialize the state of the system described by the arrays F, G, sin05phi_old and cos05phi_old
 INPUTS:    tunnel current amplitude object (pointer to TunnelCurrentType struct)
 COMMENTS:  the assumption is made that the system stayed infinitely long in the given initial state,
			that is, voltage=0 at t<0.
**************************************************************************************************/
void mitmojco_init( TunnelCurrentType *object ) {

	TunnelCurrentType_Private *s = object->self;

	int i, n;
	int Nexps = s->Nexps;
	for(i=0;i<(s->Nnodes);i++) { 
		double sin05phi_i=sin( 0.5*(*s->phi_ptrs[i]) );
		double cos05phi_i=cos( 0.5*(*s->phi_ptrs[i]) );
		s->sin05phi_old[i]=sin05phi_i;
		s->cos05phi_old[i]=cos05phi_i;
		for(n=0;n<Nexps;n++) {
			s->F[i*Nexps+n] = cos05phi_i + 0.*I;
			s->G[i*Nexps+n] = sin05phi_i + 0.*I;
		}
	}
}


/**************************************************************************************************
 NAME:      mitmojco_update
 FUNCTION:  update the state of the system (arrays F, G, sin05phi_old and cos05phi_old);
			update values of the tunnel current at the physical nodes
 INPUTS:    tunnel current amplitude object (pointer to TunnelCurrentType struct)
**************************************************************************************************/
void mitmojco_update( TunnelCurrentType *object ) {

	TunnelCurrentType_Private *s = object->self;

	int i, n;
	int Nexps = s->Nexps;
	#pragma omp parallel default(none) private(i,n) shared(s, Nexps)
	{
		#pragma omp for schedule(static)
		for(i=0;i<(s->Nnodes);i++) {
			double ReCFs=0.;
			double ReDGs=0.;
			double sin05phi_i = sin( 0.5*(*s->phi_ptrs[i]) );
			double cos05phi_i = cos( 0.5*(*s->phi_ptrs[i]) );
			for(n=0;n<Nexps;n++) {
				s->F[i*Nexps+n] = s->cexp_z[n]*s->F[i*Nexps+n] + s->alpha0[n]*s->cos05phi_old[i] + s->alpha1[n]*cos05phi_i;
				s->G[i*Nexps+n] = s->cexp_z[n]*s->G[i*Nexps+n] + s->alpha0[n]*s->sin05phi_old[i] + s->alpha1[n]*sin05phi_i;
				ReCFs += creal( s->C[n] * s->F[i*Nexps+n] );
				ReDGs += creal( s->D[n] * s->G[i*Nexps+n] );
			}
			(*s->jbar_ptrs[i]) = s->kgap_over_Rejptilde0*( sin05phi_i*ReCFs + cos05phi_i*ReDGs );
			s->sin05phi_old[i]=sin05phi_i;
			s->cos05phi_old[i]=cos05phi_i;
		}
	}
}


/**************************************************************************************************
 NAME:      mitmojco_free
 FUNCTION:  clear the allocated memory
 INPUTS:    tunnel current amplitude object (pointer to TunnelCurrentType struct)
**************************************************************************************************/
void mitmojco_free( TunnelCurrentType *object ) {
	TunnelCurrentType_Private *s = object->self;

	free(s->PubInterface.jbar);
	free(s->sin05phi_old);
	free(s->cos05phi_old);
	free(s->C);
	free(s->D);
	free(s->cexp_z);
	free(s->alpha0);
	free(s->alpha1);
	free(s->F);
	free(s->G);
	free(s);
}
