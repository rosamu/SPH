/*
 * KERNELCORRECTIONNC.cpp
 *
 *  Created on: Jul 8, 2011
 *      Author: tony
 */

#include "KERNELCORRECTIONNC.h"

namespace std {

KERNEL_CORRECTION_NC::KERNEL_CORRECTION_NC(STORAGE *store):KERNEL_CORRECTION_BASE(store) {
	// TODO Auto-generated constructor stub

}

KERNEL_CORRECTION_NC::~KERNEL_CORRECTION_NC() {
	// TODO Auto-generated destructor stub
}
void KERNEL_CORRECTION_NC::kernel_correction(int i,int j)
{
	store->frxi = store->frx;
	store->fryi = store->fry;
	store->frzi = store->frz;

	store->frxj = store->frx;
	store->fryj = store->fry;
	store->frzj = store->frz;
}

} /* namespace std */
