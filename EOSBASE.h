/*
 * EOSBASE.h
 *
 *  Created on: Jul 7, 2011
 *      Author: tony
 */

#ifndef EOSBASE_H_
#define EOSBASE_H_
#include "STORAGE.h"
class EOS_BASE {
public:
	EOS_BASE(STORAGE*);
	virtual ~EOS_BASE();
	//method
	virtual void eos(double rho,double TE,double &p,double &cs)=0;
	virtual void ThermalEnergy(double,double&,double)=0;
	//data
    double h_SWL;
    double B;
    double gamma;
    double pinf;
    double einf;
    double p0;
    double coef;
    double rho0;
    double cs0;
    STORAGE store;
};

#endif /* EOSBASE_H_ */
