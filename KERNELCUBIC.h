/*
 * KERNELCUBIC.h
 *
 *  Created on: Jul 8, 2011
 *      Author: tony
 */

#ifndef KERNELCUBIC_H_
#define KERNELCUBIC_H_

#include "KERNELBASE.h"

namespace std {

class KERNEL_CUBIC: public std::KERNEL_BASE {
public:
	KERNEL_CUBIC(STORAGE *store);
	virtual ~KERNEL_CUBIC();
	//method
	virtual void kernel(double,double,double,int,int,int,int,double);
};

} /* namespace std */
#endif /* KERNELCUBIC_H_ */
