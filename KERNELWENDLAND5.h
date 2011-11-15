/*
 * KERNELWENDLAND5.h
 *
 *  Created on: Jul 19, 2011
 *      Author: tony
 */

#ifndef KERNELWENDLAND5_H_
#define KERNELWENDLAND5_H_

#include "KERNELBASE.h"

namespace std {

class KERNEL_WENDLAND5: public std::KERNEL_BASE {
public:
	KERNEL_WENDLAND5(STORAGE *store);
	virtual ~KERNEL_WENDLAND5();
	//method
	virtual void kernel(double,double,double,int,int,int,int,double);
	//data
	double Awen, Bwen;
};

} /* namespace std */
#endif /* KERNELWENDLAND5_H_ */
