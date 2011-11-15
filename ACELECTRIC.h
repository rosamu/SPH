/*
 * ACELECTRONIC.h
 *
 *  Created on: Aug 1, 2011
 *      Author: tony
 */

#ifndef ACELECTRIC_H_
#define ACELECTRIC_H_

#include "ACBASE.h"

namespace std {

class AC_ELECTRIC: public AC_BASE {
public:
	AC_ELECTRIC(STORAGE*,SOLVE_BASE*);
	virtual ~AC_ELECTRIC();
	//method
	virtual void doit();
	virtual void self(int,int,int);
	virtual void celij(int,int,int,int,int);
};

} /* namespace std */
#endif /* ACELECTRIC_H_ */
