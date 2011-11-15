/*
 * ACBASE.h
 *
 *  Created on: Jul 7, 2011
 *      Author: tony
 */

#ifndef ACBASE_H_
#define ACBASE_H_
#include "STORAGE.h"
#include "SOLVE.h"

class AC_BASE {
public:
	AC_BASE(STORAGE*,SOLVE_BASE*);
	virtual ~AC_BASE();
	virtual void doit()=0;
	virtual void self(int,int,int)=0;
	virtual void celij(int,int,int,int,int)=0;
	//data
	STORAGE *store;
	SOLVE_BASE *solve;
	vector<double> pr;
	//for MLS
	int i_MLS_part;
	int i_ELE_part;
};

#endif /* ACBASE_H_ */
