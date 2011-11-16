/*
 * ACINTERIOR.cpp
 *
 *  Created on: Jul 8, 2011
 *      Author: tony
 */

#include "ACINTERIOR.h"
#include "KERNEL.h"
#include "math.h"

namespace std {

AC_INTERIOR::AC_INTERIOR(STORAGE *store, SOLVE_BASE *solve):AC_BASE(store,solve){
	// TODO Auto-generated constructor stub

}

AC_INTERIOR::~AC_INTERIOR() {
	// TODO Auto-generated destructor stub
}
void AC_INTERIOR::doit()
{

	store->bigUdot   = 0.0;
	if (store->dim == 3)
		store->bigVdot   = 0.0;
	store->bigWdot   = 0.0;
	store->X_Friction   = 0.0;
	if (store->dim == 3)
		store->Y_Friction   = 0.0;
	store->Z_Friction   = 0.0;
	store->nb_inFriction = 0;
	for (size_t i=0;i<store->np;i++)
	{
		store->xpdot[i] = 0;
		store->ypdot[i] = 0;
		store->zpdot[i] = 0;
		store->updot[i] = 0;
		store->vpdot[i] = 0;
		store->wpdot[i] = 0;
		store->TEdot[i] = 0;
		store->rhodot[i] = 0;
		store->xpdot[i]=0;
		store->ypdot[i]=0;
		store->zpdot[i]=0;
		store->updot[i]=0;
		store->vpdot[i]=0;
		store->wpdot[i]=0;
		store->TEdot[i]=0;
		store->rhodot[i]=0;
		store->ax[i]=0;
		store->ay[i]=0;
		store->az[i]=0;
		store->ar[i]=0;
		store->ux[i]=0;
		store->vx[i]=0;
		store->wx[i]=0;
		store->aTE[i]=0;
		store->xcor[i]=0;
		store->ycor[i]=0;
		store->zcor[i]=0;
		store->bflag[i] = 0;
	}
	//	store->dudx_CSPHo = store->dudx_CSPH;
	//	store->dudy_CSPHo = store->dudy_CSPH;
	//	store->dudz_CSPHo = store->dudz_CSPH;
	//	store->dvdx_CSPHo = store->dvdx_CSPH;
	//	store->dvdy_CSPHo = store->dvdy_CSPH;
	//	store->dvdz_CSPHo = store->dvdz_CSPH;
	//	store->dwdx_CSPHo = store->dwdx_CSPH;
	//	store->dwdy_CSPHo = store->dwdy_CSPH;
	//	store->dwdz_CSPHo = store->dwdz_CSPH;
	//	store->dTEdx_CSPHo = store->dTEdx_CSPH;
	//	store->dTEdy_CSPHo = store->dTEdy_CSPH;
	//	store->dTEdz_CSPHo = store->dTEdz_CSPH;
	//	store->drhodx_CSPHo = store->drhodx_CSPH;
	//	store->drhody_CSPHo = store->drhody_CSPH;
	//	store->drhodz_CSPHo = store->drhodz_CSPH;
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
	pr.clear();
	pr.assign(store->np,0);
	for (size_t i = 0;i<store->np;i++)
	{
		pr[i] = (store->p[i]/(store->rho[i]*store->rho[i]));
	}

	int ly2=0; //!initial zeroing of variable
	for (int kind_p1=1;kind_p1<=2;kind_p1++)
	{
		int ini_kind_p2 = kind_p1%2+1;
		if (store->dim == 3)
		{
			int local_ncn = store->local_ncn,local_ncm = store->local_ncm,local_ncl=store->local_ncl;
			//			int lzmax = store->local_ncl_interior[1],lymax = store->local_ncm_interior[1],lxmax = store->local_ncn_interior[1];
			for (int lz = 1;lz<=local_ncl;lz++)
				for (int ly = 1;ly<=local_ncm;ly++)
					for (int lx = 1;lx<=local_ncn;lx++)
					{
						size_t j1 = (lx-1) + (ly-1)*local_ncn +(lz-1)*local_ncn*local_ncm;
						//					if(store->cells[j1].nc[kind_p1-1] > 0)
						if ((kind_p1 == 1)?(store->nc_k1[j1]>0):(store->nc_k2[j1]>0))
						{
							//if the cell is not empty, then
							//loop over it and over neighboring
							//cells
							//periodic condition in x-y plane
							if (store->i_periodicOBs[0] ==1 || store->i_periodicOBs[1] == 1)
							{
								int lx2 = lx+1;
								if(lx2<local_ncn)
									celij(j1,j1+1,kind_p1,ini_kind_p2,lx2);     //!East
								else if (lx2 == local_ncn && store->i_periodicOBs[0])
								{
									celij(j1,j1+1,kind_p1,ini_kind_p2,lx2);
									store->i_periodicx = 1;
									celij(j1,(ly-1)*local_ncn +(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,lx2);
								}
								else if (store->i_periodicOBs[0])
								{
									store->i_periodicx = 1;
									celij(j1,(ly-1)*local_ncn +(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,lx2);
								}
								int ly2 = ly+1;
								if(ly2<local_ncm || (ly2 == local_ncm && !store->i_periodicOBs[1]))
								{
									celij(j1,j1+local_ncn+(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,ly2);   //!North
									lx2 = lx-1;
									if(lx2>=1) celij(j1,j1+local_ncn-1+(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,lx2);   //!N-West
									else if (store->i_periodicOBs[0])
									{
										store->i_periodicx = 1;
										celij(j1,(ly+1)*local_ncn-1+(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,lx2);
										store->i_periodicx = 1;
										celij(j1,(ly+1)*local_ncn-2+(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,lx2);
									}

									lx2 = lx+1;
									if(lx2<local_ncn) celij(j1,j1+local_ncn+1+(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,lx2); //!N-East
									else if (lx2 == local_ncn && store->i_periodicOBs[0])
									{
										celij(j1,j1+local_ncn+1+(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,lx2);
										store->i_periodicx = 1;
										celij(j1,ly*local_ncn+(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,lx2);
									}
									else if (store->i_periodicOBs[0])
									{
										store->i_periodicx = 1;
										celij(j1,ly*local_ncn+(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,lx2);
									}
								}
								else if (ly2 == local_ncm && store->i_periodicOBs[1])
								{
									celij(j1,j1+local_ncn+(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,ly2);   //!North
									store->i_periodicy = 1;
									celij(j1,lx-1+(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,ly2);
									lx2 = lx-1;
									if(lx2>1)
									{
										celij(j1,j1+local_ncn-1+(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,lx2);   //!N-West
										store->i_periodicy = 1;
										celij(j1,lx-2+(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,ly2);
									}
									else if (lx2 == 1)
									{
										celij(j1,j1+local_ncn-1+(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,lx2);
										store->i_periodicy = 1;
										celij(j1,(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,lx2);
									}
									else if (store->i_periodicOBs[0])
									{
										store->i_periodicx = 1;
										celij(j1,local_ncm*local_ncn-1+(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,lx2);
										store->i_periodicx = 1;
										celij(j1,local_ncm*local_ncn-2+(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,lx2);
										store->i_periodicx = 1;
										store->i_periodicy = 1;
										celij(j1,local_ncn-1+(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,lx2);
										store->i_periodicx = 1;
										store->i_periodicy = 1;
										celij(j1,local_ncn-2+(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,lx2);
									}

									lx2 = lx+1;
									if(lx2<local_ncn)
									{
										celij(j1,j1+local_ncn+1+(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,lx2); //!N-East
										store->i_periodicy = 1;
										celij(j1,lx+(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,lx2);
									}
									else if (lx2==local_ncn)
									{
										celij(j1,j1+local_ncn+1,kind_p1,ini_kind_p2,lx2); //!N-East
										store->i_periodicy = 1;
										celij(j1,local_ncn-1+(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,lx2);
										if (store->i_periodicOBs[0])
										{
											store->i_periodicx = 1;
											celij(j1,(local_ncm-1)*local_ncn+(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,lx2);

											store->i_periodicx = 1;
											store->i_periodicy = 1;
											celij(j1,0+(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,lx2);
										}
									}
									else if (store->i_periodicOBs[0])
									{
										store->i_periodicx = 1;
										celij(j1,(local_ncm-1)*local_ncn,kind_p1,ini_kind_p2,lx2);
										store->i_periodicx = 1;
										store->i_periodicy = 1;
										celij(j1,0+(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,lx2);
									}
								}
								else if (ly2 > local_ncm && store->i_periodicOBs[1])
								{
									store->i_periodicy = 1;
									celij(j1,lx-1+(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,ly2);   //!North
									lx2 = lx-1;
									if(lx2>=1)
									{
										store->i_periodicy = 1;
										celij(j1,lx-2+(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,lx2);   //!N-West
									}
									else if (store->i_periodicOBs[0])
									{
										store->i_periodicx = 1;
										store->i_periodicy = 1;
										celij(j1,local_ncn-1+(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,lx2);
										store->i_periodicx = 1;
										store->i_periodicy = 1;
										celij(j1,local_ncn-2+(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,lx2);

									}

									lx2 = lx+1;
									if(lx2<local_ncn)
									{
										store->i_periodicy = 1;
										celij(j1,lx+(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,lx2); //!N-East
									}
									if(lx2==local_ncn)
									{
										store->i_periodicy = 1;
										celij(j1,local_ncn-1+(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,lx2); //!N-East
										if (store->i_periodicOBs[0])
										{
											store->i_periodicx = 1;
											store->i_periodicy = 1;
											celij(j1,(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,lx2); //!N-East
										}
									}
									else if (store->i_periodicOBs[0])
									{
										store->i_periodicx = 1;
										store->i_periodicy = 1;
										celij(j1,(lz-1)*local_ncn*local_ncm,kind_p1,ini_kind_p2,lx2);
									}
								}
								//-- Cells in the next XY  sheet --
								int lz2=lz+1;
								if(lz2<=local_ncl)
								{

									//- Same row -
									celij(j1,j1+local_ncn*local_ncm,kind_p1,ini_kind_p2,ly2); //!Up

									lx2=lx-1;
									if(lx2>=1) celij(j1,j1+local_ncn*local_ncm-1,kind_p1,ini_kind_p2,ly2);   //!Up & West

									lx2=lx+1;
									if(lx2<=local_ncn) celij(j1,j1+local_ncn*local_ncm+1,kind_p1,ini_kind_p2,ly2);   //!Up & East

									//- Next row -
									ly2=ly+1;
									if(ly2<=local_ncm)
									{
										celij(j1,j1+local_ncn*local_ncm+local_ncn,kind_p1,ini_kind_p2,ly2);   //!Up & North

										lx2=lx-1;
										if(lx2>=1)
											celij(j1,j1+local_ncn*local_ncm+local_ncn-1,kind_p1,ini_kind_p2,ly2);   //!Up & North-West

										lx2=lx+1;
										if(lx2<=local_ncn)
											celij(j1,j1+local_ncn*local_ncm+local_ncn+1,kind_p1,ini_kind_p2,ly2);   //!Up & North-East
									}

									//- Previous row -
									ly2=ly-1;
									if(ly2>=1)
									{
										celij(j1,j1+local_ncn*local_ncm-local_ncn,kind_p1,ini_kind_p2,ly2);   //!Up & South

										lx2=lx-1;
										if(lx2>=1) celij(j1,j1+local_ncn*local_ncm-local_ncn-1,kind_p1,ini_kind_p2,ly2);   //!Up & SouthWest

										lx2=lx+1;
										if(lx2<=local_ncn) celij(j1,j1+local_ncn*local_ncm-local_ncn+1,kind_p1,ini_kind_p2,ly2);   //!Up & SouthEast
									}
								}
							}
							else
							{
								int lx2 = lx+1;
								if(lx2<=local_ncn)
									celij(j1,j1+1,kind_p1,ini_kind_p2,ly2);     //!East

								int ly2 = ly+1;
								if(ly2<=local_ncm)
								{
									celij(j1,j1+local_ncn,kind_p1,ini_kind_p2,ly2);   //!North
									lx2 = lx-1;
									if(lx2>=1) celij(j1,j1+local_ncn-1,kind_p1,ini_kind_p2,ly2);   //!N-West

									lx2 = lx+1;
									if(lx2<=local_ncn) celij(j1,j1+local_ncn+1,kind_p1,ini_kind_p2,ly2);   //!N-East
								}

								//-- Cells in the next XY  sheet --
								int lz2=lz+1;
								if(lz2<=local_ncl)
								{

									//- Same row -
									celij(j1,j1+local_ncn*local_ncm,kind_p1,ini_kind_p2,ly2); //!Up

									lx2=lx-1;
									if(lx2>=1) celij(j1,j1+local_ncn*local_ncm-1,kind_p1,ini_kind_p2,ly2);   //!Up & West

									lx2=lx+1;
									if(lx2<=local_ncn) celij(j1,j1+local_ncn*local_ncm+1,kind_p1,ini_kind_p2,ly2);   //!Up & East

									//- Next row -
									ly2=ly+1;
									if(ly2<=local_ncm)
									{
										celij(j1,j1+local_ncn*local_ncm+local_ncn,kind_p1,ini_kind_p2,ly2);   //!Up & North

										lx2=lx-1;
										if(lx2>=1)
											celij(j1,j1+local_ncn*local_ncm+local_ncn-1,kind_p1,ini_kind_p2,ly2);   //!Up & North-West

										lx2=lx+1;
										if(lx2<=local_ncn)
											celij(j1,j1+local_ncn*local_ncm+local_ncn+1,kind_p1,ini_kind_p2,ly2);   //!Up & North-East
									}

									//- Previous row -
									ly2=ly-1;
									if(ly2>=1)
									{
										celij(j1,j1+local_ncn*local_ncm-local_ncn,kind_p1,ini_kind_p2,ly2);   //!Up & South

										lx2=lx-1;
										if(lx2>=1) celij(j1,j1+local_ncn*local_ncm-local_ncn-1,kind_p1,ini_kind_p2,ly2);   //!Up & SouthWest

										lx2=lx+1;
										if(lx2<=local_ncn) celij(j1,j1+local_ncn*local_ncm-local_ncn+1,kind_p1,ini_kind_p2,ly2);   //!Up & SouthEast
									}
								}
							}
						}

					}
		}
		else
		{
			int local_ncn = store->local_ncn,local_ncl=store->local_ncl;
			for (int lz = 1;lz<=local_ncl;lz++)
				for (int lx = 1;lx<=local_ncn;lx++)
				{
					size_t j1 = lx + (lz-1)*local_ncn;
					j1--;
					//					if(store->cells[j1].nc[kind_p1-1] > 0)
					if ((kind_p1 == 1)?(store->nc_k1[j1]>0):(store->nc_k2[j1]>0))
					{
						//if the cell is not empty, then
						//loop over it and over neighboring
						//cells

						int lx2 = lx+1;
						if(lx2<=local_ncn)
							celij(j1,j1+1,kind_p1,ini_kind_p2,lx2);     //!East

						int lz2 = lz+1;
						if(lz2<=local_ncl)
						{
							celij(j1,j1+local_ncn,kind_p1,ini_kind_p2,lz2);   //!North
							lx2 = lx-1;
							if(lx2>=1) celij(j1,j1+local_ncn-1,kind_p1,ini_kind_p2,lx2);   //!N-West

							lx2 = lx+1;
							if(lx2<=local_ncn) celij(j1,j1+local_ncn+1,kind_p1,ini_kind_p2,lx2);   //!N-East
						}

					}
				}
		}
		//		for (size_t j1 = 0;j1<store->nx;j1++)
		//			if(store->cells[j1].nc[kind_p1-1] > 0)
		int local_ncn = store->local_ncn,local_ncm = store->local_ncm,local_ncl=store->local_ncl;
		int lzmax = store->local_ncl_interior[1],lymax = store->local_ncm_interior[1],lxmax = store->local_ncn_interior[1];
		for (int lz = store->local_ncl_interior[0]+1;lz<=lzmax;lz++)
			for (int ly = store->local_ncm_interior[0]+1;ly<=lymax;ly++)
				for (int lx = store->local_ncn_interior[0]+1;lx<=lxmax;lx++)
				{
					size_t j1 = (lx-1) + (ly-1)*local_ncn +(lz-1)*local_ncn*local_ncm;
					//		for (size_t j1 = 0;j1<store->nx;j1++)
					if ((kind_p1 == 1)?(store->nc_k1[j1]>0):(store->nc_k2[j1]>0))
						self(j1,kind_p1,ini_kind_p2);
				}
	}

	for (size_t i = 0;i<store->np;i++)
	{
		store->updot[i] = store->ax[i];
		store->vpdot[i] = store->ay[i];
		store->wpdot[i] = store->az[i];
		store->rhodot[i] = store->ar[i];
		store->TEdot[i] = store->aTE[i];
		store->xcor[i]=store->eps*store->ux[i];
		if (store->dim == 3)
			store->ycor[i]=store->eps*store->vx[i];
		store->zcor[i]=store->eps*store->wx[i];
	}


	return;
}
void AC_INTERIOR::self(int j1,int kind_p1,int ini_kind_p2)
{
	for(int kind_p2=ini_kind_p2;kind_p2<=2;kind_p2++)
	{
		if (kind_p1 > kind_p2) continue;
		vector<size_t>::iterator itjstart = ((kind_p2 == 1)?(store->ibox_k1[j1].begin()):(store->ibox_k2[j1].begin()));
		vector<size_t>::iterator iti=((kind_p1 == 1)?(store->ibox_k1[j1].begin()):(store->ibox_k2[j1].begin()));
		for(;iti!=((kind_p1 == 1)?(store->ibox_k1[j1].end()):(store->ibox_k2[j1].end()));++iti)
		{
			size_t i;
			i=*iti;
			if(kind_p1 == kind_p2) itjstart++;
			vector<size_t>::iterator itj;
			for (itj = itjstart;itj!=((kind_p2 == 1)?(store->ibox_k1[j1].end()):(store->ibox_k2[j1].end()));++itj)
			{
				size_t j;
				j=*itj;
				double drx = store->xp[i] - store->xp[j];
				double dry = (store->dim == 3)?(store->yp[i] - store->yp[j]):0;
				double drz = store->zp[i] - store->zp[j];

				double rr2 = drx*drx + dry*dry + drz*drz;

				if(rr2 < 4 * store->h * store->h && rr2 > 1.e-18)
				{
					double dux = store->up[i] - store->up[j];
					double duy = (store->dim == 3)?(store->vp[i] - store->vp[j]):0;
					double duz = store->wp[i] - store->wp[j];

					//            Calculating kernel & Normalized Kernel Gradient
					solve->kernel_func->kernel(drx,dry,drz,i,j,j1,j1,rr2);
					solve->kernel_correction_func->kernel_correction(i,j);

					//  ...  average density

					double robar  = 0.5*( store->rho[i] + store->rho[j] );
					double cbar  = 0.5*( store->cs[i] + store->cs[j] );

					//  ...  inner product rv
					//
					double dot = drx*dux + dry*duy + drz*duz;

					//	  Used to calculate the time step due to viscosity

					solve->visc_dt=std::max(dot/(rr2 + 0.01*store->h*store->h),solve->visc_dt);
					//
					//  ...  pressure and viscous force (Monaghan 1992; Ann. Rev.
					//			        Astron. Astrop. 30. Formula 3.3)
					//         pm(j) is mass of particle j
					//
					double p_v = pr[i] + pr[j];


					//  	Tensile correction (Monaghan , JCP.  159 (2000) 290- 311)
					//	Only to be activated with cubic spline kernel

					if (store->index_tensile == 1)
					{

						//              ____ Tensile correction

						double fab=solve->kernel_func->Wab*solve->kernel_func->od_Wdeltap;    //!NOTE: We'll use a non-normalized
						fab=fab*fab;           //!kernel to calculate tensile correction
						fab=fab*fab;           //!It's the (Wab/Wdeltap)**4  of Monaghan's paper
						double Ra,Rb;
						if (store->p[i] > 0)
							Ra= 0.01 *pr[i];
						else
							Ra= 0.2 *fabs(pr[i]);

						if (store->p[j]>0)
							Rb= 0.01 *pr[j];
						else
							Rb= 0.2 *fabs(pr[j]);

						p_v = p_v+ (Ra+Rb)*fab;

					}
					store->ax[i] = store->ax[i] - store->pm[j] * p_v * store->frxi;
					if (store->dim == 3)
						store->ay[i] = store->ay[i] - store->pm[j] * p_v * store->fryi;
					store->az[i] = store->az[i] - store->pm[j] * p_v * store->frzi;

					store->ax[j] = store->ax[j] + store->pm[i] * p_v * store->frxj;
					if (store->dim == 3)
						store->ay[j] = store->ay[j] + store->pm[i] * p_v * store->fryj;
					store->az[j] = store->az[j] + store->pm[i] * p_v * store->frzj;

					if (store->iRiemannSolver != 1)
					{
						store->gradients_calc(i,j,dux,duy,duz);
					}

					store->viscosity(dot,drx,dry,drz,dux,duy,duz,rr2,cbar,robar,1/robar,i,j,j1,j1);

					//   ...         Thermal Energy
					//		(Monaghan, JCP 110 (1994) 399- 406)


					double term1i=0.5 * p_v *( store->frxi*dux+store->fryi*duy+store->frzi*duz);
					double term1j=0.5 * p_v *( store->frxj*dux+store->fryj*duy+store->frzj*duz);

					store->aTE[i]=store->aTE[i]+store->pm[j]*(term1i + store->term2i);
					store->aTE[j]=store->aTE[j]+store->pm[i]*(term1j + store->term2j);


					//  ...  density acceleration (Monaghan 1992; Ann. Rev. Astron. Astrop. 30. Formula 3.9)
					//       using the derivative of the kernel, not the kernel itself

					double dot2i = dux*store->frxi + duy*store->fryi + duz*store->frzi;
					double dot2j = dux*store->frxj + duy*store->fryj + duz*store->frzj;
					store->ar[i] = store->ar[i] + store->pm[j]*dot2i;
					store->ar[j] = store->ar[j] + store->pm[i]*dot2j;

					//
					//  ...  XSPH correction (Monaghan 1994;  J. Comp. Phys. 110. Formula 2.6)
					//
					double pmj_Wab_over_rhobar = store->pm[j]*solve->kernel_func->Wab/robar;
					store->ux[i] = store->ux[i] - dux*pmj_Wab_over_rhobar;  //!pm(j) * dux * Wab / robar ! (2.6)
					if (store->dim == 3)
						store->vx[i] = store->vx[i] - duy*pmj_Wab_over_rhobar;  //!pm(j) * duy * Wab / robar ! (2.6)
					store->wx[i] = store->wx[i] - duz*pmj_Wab_over_rhobar;  //!pm(j) * duz * Wab / robar

					double pmi_Wab_over_rhobar = store->pm[i]*solve->kernel_func->Wab/robar;
					store->ux[j] = store->ux[j] + dux*pmi_Wab_over_rhobar;   //!pm(i) * dux * Wab / robar
					if (store->dim == 3)
						store->vx[j] = store->vx[j] + duy*pmi_Wab_over_rhobar;  //!pm(i) * duy * Wab / robar
					store->wx[j] = store->wx[j] + duz*pmi_Wab_over_rhobar;  //!pm(i) * duz * Wab / robar

				}
			}
		}
	}
}



void AC_INTERIOR::celij(int j1,int j2,int kind_p1,int ini_kind_p2,int ly2)
{
	for (int kind_p2 = ini_kind_p2;kind_p2<=2;kind_p2++)
	{
		if((kind_p2 == 1)?(store->nc_k1[j2]>0):(store->nc_k2[j2]>0))
		{
			vector<size_t>::iterator iti=((kind_p1 == 1)?(store->ibox_k1[j1].begin()):(store->ibox_k2[j1].begin()));
			for(;iti!=((kind_p1 == 1)?(store->ibox_k1[j1].end()):(store->ibox_k2[j1].end()));++iti)
			{
				size_t i;
				i = *iti;
				vector<size_t>::iterator itj=((kind_p2 == 1)?(store->ibox_k1[j2].begin()):(store->ibox_k2[j2].begin()));
				for (;itj!=((kind_p2 == 1)?(store->ibox_k1[j2].end()):(store->ibox_k2[j2].end()));++itj)
				{
					size_t j;
					j = *itj;
					double drx = store->xp[i] - store->xp[j];
					double dry = (store->dim == 3)?(store->yp[i] - store->yp[j]):0;
					double drz = store->zp[i] - store->zp[j];
					if (store->i_periodicOBs[0] && store->i_periodicx == 1)
					{
						drx = (drx>0)?(drx - (store->localvlx[1] - store->localvlx[0])):(drx + (store->localvlx[1] - store->localvlx[0]));
					}
					if (store->i_periodicOBs[1] && store->i_periodicy == 1)
					{
						dry = (dry>0)?(dry - (store->localvly[1] - store->localvly[0])):(dry + (store->localvly[1] - store->localvly[0]));
					}
					double rr2 = drx*drx + dry*dry + drz*drz;

					if(rr2 < 4 * store->h * store->h && rr2 > 1.e-18)
					{
						double dux = store->up[i] - store->up[j];
						double duy = (store->dim == 3)?(store->vp[i] - store->vp[j]):0;
						double duz = store->wp[i] - store->wp[j];

						//            Calculating kernel & Normalized Kernel Gradient
						solve->kernel_func->kernel(drx,dry,drz,i,j,j1,j2,rr2);
						solve->kernel_correction_func->kernel_correction(i,j);


						//  ...  average density

						double robar  = 0.5*( store->rho[i] + store->rho[j] );
						double cbar  = 0.5*( store->cs[i] + store->cs[j] );

						//  ...  inner product rv
						//
						double dot = drx*dux + dry*duy + drz*duz;

						//	  Used to calculate the time step due to viscosity

						solve->visc_dt=std::max(dot/(rr2 + 0.01*store->h*store->h),solve->visc_dt);
						//
						//  ...  pressure and viscous force (Monaghan 1992; Ann. Rev.
						//			        Astron. Astrop. 30. Formula 3.3)
						//         pm(j) is mass of particle j
						//
						double p_v = pr[i] + pr[j];


						//  	Tensile correction (Monaghan , JCP.  159 (2000) 290- 311)
						//	Only to be activated with cubic spline kernel

						if (store->index_tensile == 1)
						{

							//              ____ Tensile correction

							double fab=solve->kernel_func->Wab*solve->kernel_func->od_Wdeltap;    //!NOTE: We'll use a non-normalized
							fab=fab*fab;           //!kernel to calculate tensile correction
							fab=fab*fab;           //!It's the (Wab/Wdeltap)**4  of Monaghan's paper
							double Ra,Rb;
							if (store->p[i] > 0)
								Ra= 0.01 *pr[i];
							else
								Ra= 0.2 *fabs(pr[i]);

							if (store->p[j]>0)
								Rb= 0.01 *pr[j];
							else
								Rb= 0.2 *fabs(pr[j]);

							p_v = p_v+ (Ra+Rb)*fab;

						}

						store->ax[i] = store->ax[i] - store->pm[j] * p_v * store->frxi;
						if (store->dim == 3)
							store->ay[i] = store->ay[i] - store->pm[j] * p_v * store->fryi;
						store->az[i] = store->az[i] - store->pm[j] * p_v * store->frzi;

						store->ax[j] = store->ax[j] + store->pm[i] * p_v * store->frxj;
						if (store->dim == 3)
							store->ay[j] = store->ay[j] + store->pm[i] * p_v * store->fryj;
						store->az[j] = store->az[j] + store->pm[i] * p_v * store->frzj;
						if (store->iRiemannSolver != 1)
						{
							store->gradients_calc(i,j,dux,duy,duz);
						}

						store->viscosity(dot,drx,dry,drz,dux,duy,duz,rr2,cbar,robar,1/robar,i,j,j1,j2);

						//   ...         Thermal Energy
						//		(Monaghan, JCP 110 (1994) 399- 406)


						double term1i=0.5 * p_v *( store->frxi*dux+store->fryi*duy+store->frzi*duz);
						double term1j=0.5 * p_v *( store->frxj*dux+store->fryj*duy+store->frzj*duz);

						store->aTE[i]=store->aTE[i]+store->pm[j]*(term1i + store->term2i);
						store->aTE[j]=store->aTE[j]+store->pm[i]*(term1j + store->term2j);


						//  ...  density acceleration (Monaghan 1992; Ann. Rev. Astron. Astrop. 30. Formula 3.9)
						//       using the derivative of the kernel, not the kernel itself

						double dot2i = dux*store->frxi + duy*store->fryi + duz*store->frzi;
						double dot2j = dux*store->frxj + duy*store->fryj + duz*store->frzj;
						store->ar[i] = store->ar[i] + store->pm[j]*dot2i;
						store->ar[j] = store->ar[j] + store->pm[i]*dot2j;

						//
						//  ...  XSPH correction (Monaghan 1994;  J. Comp. Phys. 110. Formula 2.6)
						//
						double pmj_Wab_over_rhobar = store->pm[j]*solve->kernel_func->Wab/robar;
						store->ux[i] = store->ux[i] - dux*pmj_Wab_over_rhobar;  //!pm(j) * dux * Wab / robar ! (2.6)
						if (store->dim == 3)
							store->vx[i] = store->vx[i] - duy*pmj_Wab_over_rhobar;  //!pm(j) * duy * Wab / robar ! (2.6)
						store->wx[i] = store->wx[i] - duz*pmj_Wab_over_rhobar;  //!pm(j) * duz * Wab / robar

						double pmi_Wab_over_rhobar = store->pm[i]*solve->kernel_func->Wab/robar;
						store->ux[j] = store->ux[j] + dux*pmi_Wab_over_rhobar;   //!pm(i) * dux * Wab / robar
						if (store->dim == 3)
							store->vx[j] = store->vx[j] + duy*pmi_Wab_over_rhobar;  //!pm(i) * duy * Wab / robar
						store->wx[j] = store->wx[j] + duz*pmi_Wab_over_rhobar;  //!pm(i) * duz * Wab / robar

					}
				}
			}
		}
	}
	if (store->i_periodicx == 1) store->i_periodicx = 0;
	if (store->i_periodicy == 1) store->i_periodicy = 0;
}

} /* namespace std */
