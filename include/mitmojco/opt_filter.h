/* 
	Header file for optimal filtration
*/
//==============================================//
//---------- by Dmitry R. Gulevich -------------//
//----------- drgulevich@gmail.com -------------//
//==============================================//
#ifndef OPT_FILTER_HEADER
#define OPT_FILTER_HEADER

#ifdef __cplusplus
extern "C" {
#endif

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

#ifdef __cplusplus
}
#endif

#endif
