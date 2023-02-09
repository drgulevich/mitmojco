/* 
	Optimal filtration: implementation of the algorithm outlined by Eqs.(26)-(28) in Ref. 
   [1] A. A. Odintsov, V. K. Semenov and A. B. Zorin, IEEE Trans. Magn. 23, 763 (1987).
*/
//==============================================//
//---------- by Dmitry R. Gulevich -------------//
//----------- drgulevich@gmail.com -------------//
//==============================================//
#include <math.h>
#include <stdlib.h> // malloc
#include "mitmojco/opt_filter.h"

typedef struct {

	// interface to public struct
	OptFilterType PubInterface;

	// private
	int iter;
	double inva;

} OptFilterType_Private;


/**************************************************************************************************
 NAME:      opt_filter_create
 FUNCTION:  constructor of filter object
 INPUTS:    optimum filtration level (parameter n in Ref. [1])
 RETURNS:   filter object (pointer to OptFilterType struct)
**************************************************************************************************/
OptFilterType* opt_filter_create(int n) {
	OptFilterType_Private *s = malloc(sizeof(OptFilterType_Private));
	s->PubInterface.self = s; // Keep record of self pointer
	s->PubInterface.n = n;
	s->PubInterface.a = 1./(pow(2.,1./n)-1.); // optimum for sinusoidal signal
	s->PubInterface.y = malloc( (n+1) * sizeof(double) );
	s->inva = 1./(s->PubInterface.a);
	s->iter = 0;
	return &(s->PubInterface);
}


/**************************************************************************************************
 NAME:      opt_filter_init
 FUNCTION:  initialize filter object
 INPUTS:    filter object (pointer to OptFilterType struct)
**************************************************************************************************/
void opt_filter_init(OptFilterType *object) {
	OptFilterType_Private *s = object->self;
	s->iter = 0;
}


/**************************************************************************************************
 NAME:      opt_filter_update
 FUNCTION:  record new value of the signal and update the filter
 INPUTS:    filter object (pointer to OptFilterType struct)
**************************************************************************************************/
void opt_filter_update(OptFilterType *object, double signal) {
	OptFilterType_Private *s = object->self;

	int m;
	if((s->iter)==0) {
		for(m=0;m<=(object->n);m++)
			object->y[m] = signal;
	}
	else {
		object->y[0] = signal;
		double factor = s->inva * s->iter;
		double prefactor = 1./(1.+factor);
		for(m=1;m<=(object->n);m++)
			object->y[m] = prefactor*( factor * object->y[m] + object->y[m-1] );
	}
	s->iter++;
}


/**************************************************************************************************
 NAME:      opt_filter_result
 FUNCTION:  return the optimum filtration result
 INPUTS:    filter object (pointer to OptFilterType struct)
 RETURNS:   optimum filtration result
**************************************************************************************************/
double opt_filter_result(const OptFilterType *object) {
	return object->y[object->n];
}


/**************************************************************************************************
 NAME:      opt_filter_free
 FUNCTION:  clear the allocated memory
 INPUTS:    filter object (pointer to OptFilterType struct)
**************************************************************************************************/
void opt_filter_free( OptFilterType *object ) {
	OptFilterType_Private *s = object->self;
	free(s->PubInterface.y);
	free(s);
}
