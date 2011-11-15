/*
 * ACINTERIOR.h
 *
 *  Created on: Jul 8, 2011
 *      Author: tony
 */

#ifndef ACINTERIOR_H_
#define ACINTERIOR_H_

#include "ACBASE.h"

namespace std {

class AC_INTERIOR: public AC_BASE {
public:
	AC_INTERIOR(STORAGE*,SOLVE_BASE*);
	virtual ~AC_INTERIOR();
	//method
	virtual void doit();
	virtual void self(int,int,int);
	virtual void celij(int,int,int,int,int);
};

} /* namespace std */
#endif /* ACINTERIOR_H_ */
