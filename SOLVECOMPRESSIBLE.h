/*
 * SOLVECOMPRESSIBLE.h
 *
 *  Created on: Jul 8, 2011
 *      Author: tony
 */

#ifndef SOLVECOMPRESSIBLE_H_
#define SOLVECOMPRESSIBLE_H_

#include "SOLVEBASE.h"

namespace std {

class SOLVE_COMPRESSIBLE: public SOLVE_BASE {
public:
	SOLVE_COMPRESSIBLE(STORAGE*,EOS_BASE*);
	virtual ~SOLVE_COMPRESSIBLE();
	//method
	virtual void step()=0;
	virtual void initialization();
	virtual void ac();
	virtual void correct();
	virtual void variable_time_step();
	virtual void allocation();
	virtual void dellocation();
	void parallelinitialization();
	//data
	double dt2;
	double ddt_p;
	double ddt_c;
	double dt_new;
	int nstep_DBC;
	int iDBC;
};

} /* namespace std */
#endif /* SOLVECOMPRESSIBLE_H_ */
