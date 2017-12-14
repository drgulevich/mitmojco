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

typedef struct {
	bool error; // error handling
	double Rejptilde0; // normalized critical current
	double alphaN; // damping due to the normal resistance
	double *jbar; // pointer to the reduced current density
	void *self; // pointer to private struct
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
