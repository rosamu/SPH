/*
 * EOSSPOLY.h
 *
 *  Created on: Jul 7, 2011
 *      Author: tony
 */

#ifndef EOSSPOLY_H_
#define EOSSPOLY_H_

#include "EOSBASE.h"

class EOS_SPOLY: public EOS_BASE {
public:
	EOS_SPOLY(STORAGE *store);
	virtual ~EOS_SPOLY();
	virtual void eos(double,double,double&,double&);
	virtual void ThermalEnergy(double,double&,double);
};

#endif /* EOSSPOLY_H_ */
