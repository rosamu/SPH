/*
 * EOSSPOLY.cpp
 *
 *  Created on: Jul 7, 2011
 *      Author: tony
 */

#include "EOSSPOLY.h"
#include "math.h"

EOS_SPOLY::EOS_SPOLY(STORAGE *store)
:EOS_BASE(store){}

EOS_SPOLY::~EOS_SPOLY() {
	// TODO Auto-generated destructor stub
}

void EOS_SPOLY::eos(double rho,double TE,double &p,double &cs)
{
	p = (gamma - 1.0)*(TE + einf) * rho - gamma*pinf;
	cs = sqrt(gamma*(p + pinf)/rho);
	return;
}

void EOS_SPOLY::ThermalEnergy(double rho,double &TE,double p)
{
	TE = (p+gamma*pinf)/(gamma-1.0)/rho - einf;
	return;
}
