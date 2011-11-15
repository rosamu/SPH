/*
 * ACMLS.h
 *
 *  Created on: Jul 7, 2011
 *      Author: tony
 */

#ifndef ACPREMLS_H_
#define ACPREMLS_H_

#include "ACBASE.h"

class AC_MLS: public AC_BASE {
public:
	AC_MLS(STORAGE*,SOLVE_BASE*);
	virtual ~AC_MLS();
	virtual void doit();
	virtual void self(int,int,int);
	virtual void celij(int,int,int,int,int);
};

#endif /* ACPREMLS_H_ */
