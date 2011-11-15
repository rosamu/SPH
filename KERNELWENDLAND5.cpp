/*
 * KERNELWENDLAND5.cpp
 *
 *  Created on: Jul 19, 2011
 *      Author: tony
 */

#include "KERNELWENDLAND5.h"
#include "math.h"

namespace std {

KERNEL_WENDLAND5::KERNEL_WENDLAND5(STORAGE *store):KERNEL_BASE(store) {
	Awen=0.41778/(pow(store->h,3));
	Bwen=-2.08891/(pow(store->h,4));

	deltap=0.5;
	Wdeltap= Awen*(pow((1.0-deltap*0.5),4.0))*(2.0*deltap+1.0);
	od_Wdeltap=1./Wdeltap;
	store->adh=Awen;
}

KERNEL_WENDLAND5::~KERNEL_WENDLAND5() {
	// TODO Auto-generated destructor stub
}

void KERNEL_WENDLAND5::kernel(double drx,double dry,double drz,int i,int j,int j1,int j2,double drr2)
{
	store->index_tensile= 0;

	//c     Cubic 3D Kernel

	double rad = sqrt(drr2);
	double qq  = rad/store->h;
	double wqq=2*qq+1.0;
	double wqq1=1-0.5*qq;
	double wqq2=wqq1*wqq1;
	double wqq3=wqq2*wqq1;
	double wqq4=wqq3*wqq1;
	Wab= Awen*wqq*wqq4;
	double fac =Bwen*qq*wqq3 /rad;
	store->frx = fac * drx;
	if (store->dim == 3)
		store->fry = fac * dry;
	store->frz = fac * drz;



	return;
}

} /* namespace std */
