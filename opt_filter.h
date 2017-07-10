/* 
    Header file for optimal filtration
*/
//==============================================//
//---------- by Dmitry R. Gulevich -------------//
//------ d.r.gulevich@metalab.ifmo.ru ----------//
//--- ITMO University, St Petersburg, Russia ---//
//==============================================//
#ifndef OPT_FILTER_HEADER
#define OPT_FILTER_HEADER

typedef struct {
	int n;
    double a;
    double *y;
	void *self; // pointer to private struct
} OptFilterType;

extern OptFilterType* opt_filter_create( int n );
extern void opt_filter_init( OptFilterType *object );
extern void opt_filter_update( OptFilterType *object, double signal );
extern double opt_filter_result( const OptFilterType *object );
extern void opt_filter_free( OptFilterType *object );

#endif
