/*
 * EOSTAIT.cpp
 *
 *  Created on: Jul 7, 2011
 *      Author: tony
 */

#include "EOSTAIT.h"
#include "math.h"

EOS_TAIT::EOS_TAIT(STORAGE *store):EOS_BASE(store){
	// TODO Auto-generated constructor stub

}

EOS_TAIT::~EOS_TAIT() {
	// TODO Auto-generated destructor stub
}
void EOS_TAIT::eos(double rho,double TE,double &p,double &cs)
{
    p =  B * ( pow((rho/rho0),(int)gamma) - 1);
    cs = cs0*pow((rho/rho0),3);
}
