/*
 * SOLVEPREDICTORCORRECTOR.cpp
 *
 *  Created on: Jul 7, 2011
 *      Author: tony
 */

#include "SOLVEPREDICTORCORRECTOR.h"

SOLVE_PREDICTOR_CORRECTOR::SOLVE_PREDICTOR_CORRECTOR(STORAGE *store,EOS_BASE *eos):SOLVE_COMPRESSIBLE(store,eos)
{
	// TODO Auto-generated constructor stub

}

SOLVE_PREDICTOR_CORRECTOR::~SOLVE_PREDICTOR_CORRECTOR() {
	// TODO Auto-generated destructor stub
}

void SOLVE_PREDICTOR_CORRECTOR::step()
{
	MPI_Barrier(MPI_COMM_WORLD);
	dt2 = 0.5*dt;

	//	Predictor- Corrector Method

	//  ...  compute acceleration
	//if (store->i_densityFilter > 0 && itime/store->ndt_FilterPerform*store->ndt_FilterPerform == itime)
	//	densityFilter();
	//store->divide(store->nb,store->np,2);
	ac();


	//c
	//c  ... compute corrections to:
	//c         (a) rate of change of velocity due to
	//c             (i)  body forces, and
	//c             (ii) boundary forces
	//c
	//c         (b) rate of change of position (free-surface, XSPH)
	//c
	//c
	correct();


	//c	To calculate time step
	if (store->ivar_dt == 1)
	{
		if (store->iRiemannSolver != 1)
			variable_time_step();
		ddt_p=dt_new;
	}

	//c
	//c  ... predictor
	//c
	vector<double>::size_type i = 0;
	for (i = 0;i<store->np;i++)
	{
		store->xpo[i]=store->xp[i];
		store->ypo[i]=store->yp[i];
		store->zpo[i]=store->zp[i];
		store->upo[i]=store->up[i];
		store->vpo[i]=store->vp[i];
		store->wpo[i]=store->wp[i];
		store->TEo[i]=store->TE[i];
		store->rhoo[i]=store->rho[i];
	}
	for (i= store->nb;i<store->np;i++)
	{
		store->xp[i] = store->xpo[i] + dt2*store->xpdot[i];
		store->yp[i] = store->ypo[i] + dt2*store->ypdot[i];
		store->zp[i] = store->zpo[i] + dt2*store->zpdot[i];
		store->up[i] = store->upo[i] + dt2*store->updot[i];
		store->vp[i] = store->vpo[i] + dt2*store->vpdot[i];
		store->wp[i] = store->wpo[i] + dt2*store->wpdot[i];
	}
	for (i= store->nb;i<store->np;i++)
	{
		store->rho[i] = store->rhoo[i] + dt2*store->rhodot[i];
		store->TE[i]  = store->TEo[i]  + dt2*store->TEdot[i];
		equation_of_state->eos(store->rho[i],store->TE[i],store->p[i],store->cs[i]);
	}


	if(store->IBC == 2)//	!Hughes & Grahams 2009 correction
	{
		nstep_DBC=nstep_DBC+1;
		iDBC=0;
		if(nstep_DBC%store->ndt_DBCPerform == 0)
		{
			for (i = 0;i<store->nbb1;i++)
			{
				store->rho[i] = store->rhoo[i] + dt2*store->rhodot[i];
				if (store->rho[i]<store->rho0) store->rho[i] = store->rho0;
				store->TE[i]  = store->TEo[i]  + dt2*store->TEdot[i];
				equation_of_state->eos(store->rho[i],store->TE[i],store->p[i],store->cs[i]);
			}
			iDBC=1;
			nstep_DBC=0;
		}
	}
	dellocation();

	//c
	//c  ...  corrector
	//c
	MPI_Barrier(MPI_COMM_WORLD);
	store->check_limits();
	store->i_SR_step = 1;
	store->SendRecvBuffer();
	if (store->NumOfCpu > 1)
		parallelinitialization();
	allocation();
	//store->divide(store->nb,store->np-1,2);
	ac();
	correct();


	if (store->ivar_dt==1)
	{
		if (store->iRiemannSolver != 1)
			variable_time_step();
		ddt_c=dt_new;
		dt=min(ddt_p,ddt_c);// !Time step to be used in next loop
	}
	TimeSync();

	for (i = store->nb; i<store->np;i++)
	{
		store->rho[i] = store->rhoo[i] + dt2*store->rhodot[i];
		store->TE[i]  = store->TEo[i]  + dt2*store->TEdot[i];
	}

	if( iDBC == 1) //	!Hughes & Grahams 2009 correction
	{
		for (i = 0; i<store->nbb1;i++)
		{
			store->rho[i] = store->rhoo[i] + dt2*store->rhodot[i];
			if (store->rho[i]<store->rho0) store->rho[i] = store->rho0;
			store->TE[i]  = store->TEo[i]  + dt2*store->TEdot[i];
		}
	}


	for (i = store->nb;i<store->np;i++)
	{
		store->xp[i] = store->xpo[i] + dt2*store->xpdot[i];
		store->yp[i] = store->ypo[i] + dt2*store->ypdot[i];
		store->zp[i] = store->zpo[i] + dt2*store->zpdot[i];
		store->up[i] = store->upo[i] + dt2*store->updot[i];
		store->vp[i] = store->vpo[i] + dt2*store->vpdot[i];
		store->wp[i] = store->wpo[i] + dt2*store->wpdot[i];
	}

	//	!-- Perform final integration correction --
	for (i = store->nb;i<store->np;i++)
	{
		store->rho[i] = 2*store->rho[i] - store->rhoo[i];
		store->TE[i] = 2*store->TE[i] - store->TEo[i];
		equation_of_state->eos(store->rho[i],store->TE[i],store->p[i],store->cs[i]);
	}


	if(iDBC == 1)  //!Hughes & Grahams 2009 correction
	{
		for (i = 0; i<store->nbb1;i++)
		{
			store->rho[i] = 2*store->rho[i] - store->rhoo[i];
			if (store->rho[i]<store->rho0) store->rho[i] = store->rho0;
			store->TE[i]  = 2*store->TE[i]  - store->TEo[i];
			equation_of_state->eos(store->rho[i],store->TE[i],store->p[i],store->cs[i]);
		}
	}

	for (i = store->nb;i<store->np;i++)
	{
		store->xp[i] = 2*store->xp[i] - store->xpo[i];
		store->yp[i] = 2*store->yp[i] - store->ypo[i];
		store->zp[i] = 2*store->zp[i] - store->zpo[i];
		store->up[i] = 2*store->up[i] - store->upo[i];
		store->vp[i] = 2*store->vp[i] - store->vpo[i];
		store->wp[i] = 2*store->wp[i] - store->wpo[i];
	}
	MPI_Barrier(MPI_COMM_WORLD);
	store->check_limits();
	store->i_SR_step = 2;
	store->SendRecvBuffer();
	if (store->NumOfCpu > 1)
		parallelinitialization();
	MPI_Barrier(MPI_COMM_WORLD);
}
