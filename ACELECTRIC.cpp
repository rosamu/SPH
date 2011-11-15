/*
 * ACELECTRONIC.cpp
 *
 *  Created on: Aug 1, 2011
 *      Author: tony
 */

#include "ACELECTRIC.h"
#include "KERNEL.h"
#include "math.h"

namespace std {

AC_ELECTRIC::AC_ELECTRIC(STORAGE *store,SOLVE_BASE *solve):AC_BASE(store,solve)  {
	// TODO Auto-generated constructor stub

}

AC_ELECTRIC::~AC_ELECTRIC() {
	// TODO Auto-generated destructor stub
}

void AC_ELECTRIC::doit()
{
	int ly2=0; //!initial zeroing of variable
	int kind_p1 = 2;
	int ini_kind_p2=kind_p1%2+1;// !Gives the alternative value to kind_p2 of kind_p1
	int ncn = store->ncn,ncm = store->ncm,ncl=store->ncl;
	for (int lz = 1;lz<=ncl;lz++)
		for (int ly = 1;ly<=ncm;ly++)
			for (int lx = 1;lx<=ncn;lx++)
			{
				size_t j1 = lx + (ly-1)*ncn +(lz-1)*ncn*ncm;
				j1--;
				//if(store->cells[j1].nc[kind_p1-1] > 0)
				if ((kind_p1 == 1)?(store->nc_k1[j1]>0):(store->nc_k2[j1]>0))
				{
					//c                            ! if the cell is not empty, then
					//c                            ! loop over it and over neighboring
					//c                            ! cells
					int lx2 = lx+1;
					if(lx2<=ncn)
						celij(j1,j1+1,kind_p1,ini_kind_p2,ly2);     //!East

					int ly2 = ly+1;
					if(ly2<=ncm)
					{
						celij(j1,j1+ncn,kind_p1,ini_kind_p2,ly2);   //!North
						lx2 = lx-1;
						if(lx2>=1) celij(j1,j1+ncn-1,kind_p1,ini_kind_p2,ly2);   //!N-West

						lx2 = lx+1;
						if(lx2<=ncn) celij(j1,j1+ncn+1,kind_p1,ini_kind_p2,ly2);   //!N-East
					}

					//c          -- Cells in the next XY  sheet --
					int lz2=lz+1;
					if(lz2<=ncl)
					{

						//          !- Same row -
						celij(j1,j1+ncn*ncm,kind_p1,ini_kind_p2,ly2); //!Up

						lx2=lx-1;
						if(lx2>=1) celij(j1,j1+ncn*ncm-1,kind_p1,ini_kind_p2,ly2);   //!Up & West

						lx2=lx+1;
						if(lx2<=ncn) celij(j1,j1+ncn*ncm+1,kind_p1,ini_kind_p2,ly2);   //!Up & East

						//           !- Next row -
						ly2=ly+1;
						if(ly2<=ncm)
						{
							celij(j1,j1+ncn*ncm+ncn,kind_p1,ini_kind_p2,ly2);   //!Up & North

							lx2=lx-1;
							if(lx2>=1)
								celij(j1,j1+ncn*ncm+ncn-1,kind_p1,ini_kind_p2,ly2);   //!Up & North-West

							lx2=lx+1;
							if(lx2<=ncn)
								celij(j1,j1+ncn*ncm+ncn+1,kind_p1,ini_kind_p2,ly2);   //!Up & North-East
						}

						//           !- Previous row -
						ly2=ly-1;
						if(ly2>=1)
						{
							celij(j1,j1+ncn*ncm-ncn,kind_p1,ini_kind_p2,ly2);   //!Up & South

							lx2=lx-1;
							if(lx2>=1) celij(j1,j1+ncn*ncm-ncn-1,kind_p1,ini_kind_p2,ly2);   //!Up & SouthWest

							lx2=lx+1;
							if(lx2<=ncn) celij(j1,j1+ncn*ncm-ncn+1,kind_p1,ini_kind_p2,ly2);   //!Up & SouthEast
						}
					}
				}
			}
	for (size_t j1 = 0;j1<store->nx;j1++)
	  if ((kind_p1 == 1)?(store->nc_k1[j1]>0):(store->nc_k2[j1]>0))
	    self(j1,kind_p1,ini_kind_p2);
	return;
}
void AC_ELECTRIC::self(int j1,int kind_p1,int ini_kind_p2)
{
	int kind_p2 = 2;
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
			solve->kernel_func->kernel(drx,dry,drz,i,j,j1,j1,rr2);
			solve->kernel_correction_func->kernel_correction(i,j);

			if(rr2 < 4 * store->h * store->h && rr2 > 1.e-18)
			{
				if (i_ELE_part == 1)
				{
					double ubi[3],ubj[3];
					store->GetUxB(ubi,i);
					store->GetUxB(ubj,j);
					double temp[3];
					temp[0] = ubi[0]/(store->rho[i]*store->rho[i]) + ubj[0]/(store->rho[j]*store->rho[j]);
					temp[1] = ubi[1]/(store->rho[i]*store->rho[i]) + ubj[1]/(store->rho[j]*store->rho[j]);
					temp[2] = ubi[2]/(store->rho[i]*store->rho[i]) + ubj[2]/(store->rho[j]*store->rho[j]);
					store->dUxB[i]+=store->rho[i]*store->pm[j]*(temp[0]*store->frxi+temp[1]*store->fryi+temp[2]*store->frzi);
					store->dUxB[j]+=-store->rho[j]*store->pm[i]*(temp[0]*store->frxj+temp[1]*store->fryj+temp[2]*store->frzj);
				}
				else if (i_ELE_part == 2)
				{
					double index;
					index = 4 / (store->rho[i] + store->rho[j]);
					index = index * store->frxi*drx+store->fryi*dry+store->frzi*drz / (rr2 + 0.01*store->h*store->h);
					store->m_pSolver->Add_A(i, i, index * store->pm[j]);
					store->m_pSolver->Add_A(i, j, -index * store->pm[j]);
					store->m_pSolver->Add_A(j, j, index * store->pm[i]);
					store->m_pSolver->Add_A(j, i, -index * store->pm[i]);
				}
				else
				{
					double index;
					index=(store->Phi[i]/store->rho[i]+store->Phi[j]/store->rho[j]);
					store->GPx[i]+=store->pm[j]*index*store->frxi;
					store->GPy[i]+=store->pm[j]*index*store->fryi;
					store->GPz[i]+=store->pm[j]*index*store->frzi;
					store->GPx[j]+=-store->pm[i]*index*store->frxj;
					store->GPy[j]+=-store->pm[i]*index*store->fryj;
					store->GPz[j]+=-store->pm[i]*index*store->frzj;
				}
			}
		}
	}
}
void AC_ELECTRIC::celij(int j1,int j2,int kind_p1,int ini_kind_p2,int ly2)
{
	double kind_p2 = 2;
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

				double rr2 = drx*drx + dry*dry + drz*drz;

				if(rr2 < 4 * store->h * store->h && rr2 > 1.e-18)
				{
					solve->kernel_func->kernel(drx,dry,drz,i,j,j1,j1,rr2);
					solve->kernel_correction_func->kernel_correction(i,j);

					if(rr2 < 4 * store->h * store->h && rr2 > 1.e-18)
					{
						if (i_ELE_part == 1)
						{
							double ubi[3],ubj[3];
							store->GetUxB(ubi,i);
							store->GetUxB(ubj,j);
							double temp[3];
							temp[0] = ubi[0]/(store->rho[i]*store->rho[i]) + ubj[0]/(store->rho[j]*store->rho[j]);
							temp[1] = ubi[1]/(store->rho[i]*store->rho[i]) + ubj[1]/(store->rho[j]*store->rho[j]);
							temp[2] = ubi[2]/(store->rho[i]*store->rho[i]) + ubj[2]/(store->rho[j]*store->rho[j]);
							store->dUxB[i]+=store->rho[i]*store->pm[j]*(temp[0]*store->frxi+temp[1]*store->fryi+temp[2]*store->frzi);
							store->dUxB[j]+=-store->rho[j]*store->pm[i]*(temp[0]*store->frxj+temp[1]*store->fryj+temp[2]*store->frzj);
						}
						else if (i_ELE_part == 2)
						{
							double index;
							index = 4 / (store->rho[i] + store->rho[j]);
							index = index * store->frxi*drx+store->fryi*dry+store->frzi*drz / (rr2 + 0.01*store->h*store->h);
							store->m_pSolver->Add_A(i, i, index * store->pm[j]);
							store->m_pSolver->Add_A(i, j, -index * store->pm[j]);
							store->m_pSolver->Add_A(j, j, index * store->pm[i]);
							store->m_pSolver->Add_A(j, i, -index * store->pm[i]);
						}
						else
						{
							double index;
							index=(store->Phi[i]/store->rho[i]+store->Phi[j]/store->rho[j]);
							store->GPx[i]+=store->pm[j]*index*store->frxi;
							store->GPy[i]+=store->pm[j]*index*store->fryi;
							store->GPz[i]+=store->pm[j]*index*store->frzi;
							store->GPx[j]+=-store->pm[i]*index*store->frxj;
							store->GPy[j]+=-store->pm[i]*index*store->fryj;
							store->GPz[j]+=-store->pm[i]*index*store->frzj;
						}
					}
				}
			}
		}
	}
}
} /* namespace std */
