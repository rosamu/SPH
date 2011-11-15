/*
 * EOSBASE.cpp
 *
 *  Created on: Jul 7, 2011
 *      Author: tony
 */

#include "EOSBASE.h"
#include "math.h"

EOS_BASE::EOS_BASE(STORAGE *store) {
	h_SWL = store->h_SWL;
	B = store->B;
	gamma = store->gamma;
	pinf = store->pinf;
	einf = store->einf;
	p0 = store->p0;
	coef = store->coef;
	rho0 = store->rho0;
	cs0 = sqrt(gamma*B/rho0);
}

EOS_BASE::~EOS_BASE() {
	// TODO Auto-generated destructor stub
}
