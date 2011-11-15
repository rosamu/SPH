/*
 * SOLVECOMPRESSIBLE.cpp
 *
 *  Created on: Jul 8, 2011
 *      Author: tony
 */

#include "SOLVECOMPRESSIBLE.h"
#include "AC.h"
#include "math.h"

namespace std {

SOLVE_COMPRESSIBLE::SOLVE_COMPRESSIBLE(STORAGE *store,EOS_BASE *eos):SOLVE_BASE(store,eos),nstep_DBC(0) {
	// TODO Auto-generated constructor stub

}

SOLVE_COMPRESSIBLE::~SOLVE_COMPRESSIBLE(){
	// TODO Auto-generated destructor stub
}



void SOLVE_COMPRESSIBLE::initialization()
{
	vector<double>::size_type i;
	double TE;
	store->grx = store->gry = 0;
	store->grz = -store->g;
	store->eps = 0.5;
	ddt_p = dt;
	ddt_c = dt;
	for (i = 0;i<store->np;i++)
	{
				equation_of_state->ThermalEnergy(store->rho[i],store->TE[i],store->p[i]);
		//		store->fld.TE.push_back(TE);
				equation_of_state->eos(store->rho[i],store->TE[i],store->p[i],store->cs[i]);
		if (store->i_EoS == 1)
		{
			double ZZmax;
			if (store->zp[i] <=2)
				ZZmax = 2;
			else
				ZZmax = 5;
			store->rho[i]=store->rho0*pow((1.0+store->rho0*store->g*(ZZmax-store->zp[i])/store->B),(1/store->gamma));
			store->p[i]=store->B * ( pow((store->rho[i]/store->rho0),(int)(store->gamma)) - 1);
		}
		
		//	}
		//	for (size_t i=0;i<store->np;i++)
		//	{
		//	}
		//	store->dudx_CSPH.assign(store->np,0);
		//	store->dudy_CSPH.assign(store->np,0);
		//	store->dudz_CSPH.assign(store->np,0);
		//	store->dvdx_CSPH.assign(store->np,0);
		//	store->dvdy_CSPH.assign(store->np,0);
		//	store->dvdz_CSPH.assign(store->np,0);
		//	store->dwdx_CSPH.assign(store->np,0);
		//	store->dwdy_CSPH.assign(store->np,0);
		//	store->dwdz_CSPH.assign(store->np,0);
		//	store->dTEdx_CSPH.assign(store->np,0);
		//	store->dTEdy_CSPH.assign(store->np,0);
		//	store->dTEdz_CSPH.assign(store->np,0);
		//	store->drhodx_CSPH.assign(store->np,0);
		//	store->drhody_CSPH.assign(store->np,0);
		//	store->drhodz_CSPH.assign(store->np,0);
		//	store->swpdot.xpdot.assign(store->np,0);
		//	store->swpdot.ypdot.assign(store->np,0);
		//	store->swpdot.zpdot.assign(store->np,0);
	}
}
void SOLVE_COMPRESSIBLE::ac()
{
	AC_BASE *ac;
	//	computeElectricFieldPotential();
	//	store->i_Print_J = 1;
	//	store->print_out(itime,time);
	if (store->iRiemannSolver != 1)
		ac = new AC_INTERIOR(store,this);
	else
	{}
	ac->doit();
	delete ac;
}

void SOLVE_COMPRESSIBLE::correct()
{

	vector<double>::size_type i;
	for (i = store->nb;i<store->np;i++)
	{
		store->updot[i] += store->grx;
		if (store->dim == 3)
			store->vpdot[i] += store->gry;
		store->wpdot[i] += store->grz;
	}
	//c
	//c  ...  account for XSPH
	//c

	for (i = store->nb;i<store->np;i++)
	{
		store->xpdot[i] = store->up[i] + store->xcor[i];
		if (store->dim == 3)
			store->ypdot[i] = store->vp[i] + store->ycor[i];
		store->zpdot[i] = store->wp[i] + store->zcor[i];
	}

}

void SOLVE_COMPRESSIBLE::variable_time_step()
{
	double Fa_max=-1000;
	double cs_max=-1000;
	double vel_max = -1000;
	vector<double>::size_type i_fa = -1000;
	vector<double>::size_type nstart_VTS;
	if(store->IBC == 1)
		nstart_VTS = 10;
	else
		nstart_VTS = store->nb;
	vector<double>::size_type i;
	for (i = nstart_VTS;i <store->np; i++)
	{
		double Fa = sqrt(store->ax[i]*store->ax[i]+store->ay[i]*store->ay[i]+store->az[i]*store->az[i]);
		if (Fa >= Fa_max)
		{
			i_fa = i;
			Fa_max = Fa;
		}
		if (store->cs[i] >= cs_max)
			cs_max = store->cs[i];
		double vel = sqrt(store->up[i]*store->up[i]+store->vp[i]*store->vp[i]+store->wp[i]*store->wp[i]);
		if (vel > vel_max)
			vel_max = vel;
	}

	//c      Fa_sqrt=((ax(i_fa)+grx)**2+(ay(i_fa)+gry)**2 +
	//c     +         (az(i_fa)+grz)**2)**0.25
	//c      dt_1=sqrt(h)/Fa_sqrt
	double dt_1=0,dt_2=0;
	if (i_fa < 0)
		assert(false);
	else
	{
		i = i_fa;
		dt_1=sqrt(store->h)/sqrt(sqrt((store->ax[i]+store->grx)*(store->ax[i]+store->grx)+(store->ay[i]+store->gry)*(store->ay[i]+store->gry) +(store->az[i]+store->grz)*(store->az[i]+store->grz)));
	}
	//    !dt_2=h/(cs_max+h*visc_dt)
	//    if(cs_max.gt.0.0)then
	//      dt_2=h/(cs_max+h*visc_dt)
	//    else
	//      dt_2=dt_1
	//    endif
	if (cs_max > 0)
		dt_2=store->h/(cs_max+store->h*visc_dt);
	else
			dt_2 = dt_1;
//	dt_2 = store->dx / vel_max;
	dt_new=store->CFL_number*min(dt_1,dt_2);

	return;
}

void SOLVE_COMPRESSIBLE::allocation()
{
}

void SOLVE_COMPRESSIBLE::dellocation()
{
}

void SOLVE_COMPRESSIBLE::parallelinitialization()
{
	double temp;
	for (size_t i = 0; i < store->np; i++)
	{
		equation_of_state->eos(store->rho[i],store->TE[i],store->p[i],store->cs[i]);
	}
}
} /* namespace std */
