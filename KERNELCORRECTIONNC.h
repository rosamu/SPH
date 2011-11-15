/*
 * KERNELCORRECTIONNC.h
 *
 *  Created on: Jul 8, 2011
 *      Author: tony
 */

#ifndef KERNELCORRECTIONNC_H_
#define KERNELCORRECTIONNC_H_

#include "KERNELCORRECTIONBASE.h"

namespace std {

class KERNEL_CORRECTION_NC: public std::KERNEL_CORRECTION_BASE {
public:
	KERNEL_CORRECTION_NC(STORAGE *store);
	virtual ~KERNEL_CORRECTION_NC();
	//method
	virtual void kernel_correction(int,int);
};

} /* namespace std */
#endif /* KERNELCORRECTIONNC_H_ */
