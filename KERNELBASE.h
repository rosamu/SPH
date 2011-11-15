/*
 * KERNELBASE.h
 *
 *  Created on: Jul 8, 2011
 *      Author: tony
 */

#ifndef KERNELBASE_H_
#define KERNELBASE_H_
#include "STORAGE.h"

namespace std {

class KERNEL_BASE {
public:
	KERNEL_BASE(STORAGE *store);
	virtual ~KERNEL_BASE();
	//method
	virtual void kernel(double,double,double,int,int,int,int,double)=0;
	STORAGE *store;
	//data
	double a1,a2,aa,a24,b1,c1,d1,e1,f1,c2,deltap,Wdeltap,od_Wdeltap,Wab;
};

} /* namespace std */
#endif /* KERNELBASE_H_ */
