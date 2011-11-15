/*
 * KERNELCUBIC.cpp
 *
 *  Created on: Jul 8, 2011
 *      Author: tony
 */

#include "KERNELCUBIC.h"
#include "math.h"

namespace std {

KERNEL_CUBIC::KERNEL_CUBIC(STORAGE *store):KERNEL_BASE(store) {
	a1 = 0.25/atan(1);
	a2 = a1/(store->h*store->h*store->h);
	aa = a1/(store->h*store->h*store->h*store->h);

	a24 = 0.25*a2;

	c1 = -3*aa;
	d1 = 9*aa/4;
	c2 = -3*aa/4;
	deltap=1/1.5;
	Wdeltap=a2*(1.-1.5*deltap*deltap+0.75*deltap*deltap*deltap);
	od_Wdeltap=1./Wdeltap;
	store->adh = a2;
}

KERNEL_CUBIC::~KERNEL_CUBIC() {
	// TODO Auto-generated destructor stub
}
void KERNEL_CUBIC::kernel(double drx,double dry,double drz,int i,int j,int j1,int j2,double drr2)
{
	store->index_tensile= 1;

	//c     Cubic 3D Kernel

	double rad = sqrt(drr2);
	double qq  = rad/store->h;

	if(drr2 > store->h*store->h)
	{
		double tmq = 2.-qq;
		Wab = a24*tmq*tmq*tmq;

		double fac = c2*tmq*tmq /rad;
		store->frx = fac * drx;
		if (store->dim == 3)
			store->fry = fac * dry;
		store->frz = fac * drz;
	}
	else
	{
		Wab = a2*( 1. - 1.5*qq*qq + 0.75*qq*qq*qq );

		double fac = (c1*qq + d1*qq*qq) /rad;
		store->frx = fac * drx;
		if (store->dim == 3)
			store->fry = fac * dry;
		store->frz = fac * drz;
	}


	return;
}

} /* namespace std */
