/*
 * EOSTAIT.h
 *
 *  Created on: Jul 7, 2011
 *      Author: tony
 */

#ifndef EOSTAIT_H_
#define EOSTAIT_H_

#include "EOSBASE.h"

class EOS_TAIT: public EOS_BASE {
public:
	EOS_TAIT(STORAGE *store);
	virtual ~EOS_TAIT();
	//method
	virtual void eos(double,double,double&,double&);
	virtual void ThermalEnergy(double rho,double &TE,double p){TE = 0;};
};

#endif /* EOSTAIT_H_ */
