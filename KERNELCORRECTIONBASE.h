/*
 * KERNELCORRECTIONBASE.h
 *
 *  Created on: Jul 8, 2011
 *      Author: tony
 */

#ifndef KERNELCORRECTIONBASE_H_
#define KERNELCORRECTIONBASE_H_
#include "STORAGE.h"

namespace std {

class KERNEL_CORRECTION_BASE {
public:
	KERNEL_CORRECTION_BASE(STORAGE *store);
	virtual ~KERNEL_CORRECTION_BASE();
	//method
	virtual void kernel_correction(int,int)=0;
	STORAGE *store;
};

} /* namespace std */
#endif /* KERNELCORRECTIONBASE_H_ */
