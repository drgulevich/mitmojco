/* 
	MiTMoJCo header file
*/
//==============================================//
//---------- by Dmitry R. Gulevich -------------//
//------ d.r.gulevich@metalab.ifmo.ru ----------//
//--- ITMO University, St Petersburg, Russia ---//
//==============================================//
#ifndef MITMOJCO_HEADER
#define MITMOJCO_HEADER

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h> // boolean data type
#include <complex.h> // complex numbers

typedef struct {
	double *sin05phi_old; // Nnodes entries
	double *cos05phi_old; // Nnodes entries
	double _Complex *F; // Nnodes*Nexps entries
	double _Complex *G; // Nnodes*Nexps entries
} MemState; // struct containing previous evolution information 

typedef struct {
	const char *filename;
	double a_supp;
	double kgap;
	double dt;
	int Ntotal;
	double *phi;
	int Nskip;
	int *skipinds;
	int Nnodes; // number of active nodes, Nnodes = Ntotal-Nskip
	int Nexps; // number of exponentials used in fitting
	MemState memstate; // struct containing previous evolution information 
	double Rejptilde0; // normalized critical current
	double alphaN; // damping due to the normal resistance
	double *jbar; // pointer to the reduced current density
	void *self; // pointer to private struct (for internal use only)
	bool error; // error handling
} TunnelCurrentType;

/* methods */
TunnelCurrentType* mitmojco_create(const char *filename, double a_supp, double kgap, double dt, int Ntotal, double *phi, int Nskip, int *skipinds);
extern void mitmojco_init( TunnelCurrentType *object );
extern void mitmojco_update( TunnelCurrentType *object );
extern void mitmojco_free( TunnelCurrentType *object );

#ifdef __cplusplus
}
#endif

#endif
