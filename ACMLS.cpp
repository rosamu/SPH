/*
 * ACMLS.cpp
 *
 *  Created on: Jul 7, 2011
 *      Author: tony
 */

#include "ACMLS.h"
#include "KERNEL.h"
#include "math.h"
#include "iostream"

AC_MLS::AC_MLS(STORAGE *store,SOLVE_BASE *solve):AC_BASE(store,solve){
	// TODO Auto-generated constructor stub

}

AC_MLS::~AC_MLS() {
	// TODO Auto-generated destructor stub
}

void AC_MLS::doit()
{
	int ly2=0; //!initial zeroing of variable

	for (int kind_p1=1;kind_p1<=2;kind_p1++)
	{
		int ini_kind_p2=kind_p1%2+1;// !Gives the alternative value to kind_p2 of kind_p1
		if (store->dim == 3)
		{
			int local_ncn = store->local_ncn,local_ncm = store->local_ncm,local_ncl=store->local_ncl;
			//			int lzmax = store->local_ncl_interior[1],lymax = store->local_ncm_interior[1],lxmax = store->local_ncn_interior[1];
			for (int lz = 1;lz<=local_ncl;lz++)
				for (int ly = 1;ly<=local_ncm;ly++)
					for (int lx = 1;lx<=local_ncn;lx++)
					{
						size_t j1 = (lx-1) + (ly-1)*local_ncn +(lz-1)*local_ncn*local_ncm;
						//if(store->cells[j1].nc[kind_p1-1] > 0)
						if ((kind_p1 == 1)?(store->nc_k1[j1]>0):(store->nc_k2[j1]>0))
						{
							//c                            ! if the cell is not empty, then
							//c                            ! loop over it and over neighboring
							//c                            ! cells
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

							//c          -- Cells in the next XY  sheet --
							int lz2=lz+1;
							if(lz2<=local_ncl)
							{

								//          !- Same row -
								celij(j1,j1+local_ncn*local_ncm,kind_p1,ini_kind_p2,ly2); //!Up

								lx2=lx-1;
								if(lx2>=1) celij(j1,j1+local_ncn*local_ncm-1,kind_p1,ini_kind_p2,ly2);   //!Up & West

								lx2=lx+1;
								if(lx2<=local_ncn) celij(j1,j1+local_ncn*local_ncm+1,kind_p1,ini_kind_p2,ly2);   //!Up & East

								//           !- Next row -
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

								//           !- Previous row -
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
		else
		{
			int local_ncn = store->local_ncn,local_ncl=store->local_ncl;
			for (int lz = 1;lz<=local_ncl;lz++)
				for (int lx = 1;lx<=local_ncn;lx++)
				{
					size_t j1 = lx + (lz-1)*local_ncn;
					j1--;
					//if(store->cells[j1].nc[kind_p1-1] > 0)
					if ((kind_p1 == 1)?(store->nc_k1[j1]>0):(store->nc_k2[j1]>0))
					{
						//c                            ! if the cell is not empty, then
						//c                            ! loop over it and over neighboring
						//c                            ! cells
						int lx2 = lx+1;
						if(lx2<=local_ncn)
							celij(j1,j1+1,kind_p1,ini_kind_p2,ly2);     //!East

						int lz2 = lz+1;
						if(lz2<=local_ncl)
						{
							celij(j1,j1+local_ncn,kind_p1,ini_kind_p2,ly2);   //!North
							lx2 = lx-1;
							if(lx2>=1) celij(j1,j1+local_ncn-1,kind_p1,ini_kind_p2,ly2);   //!N-West

							lx2 = lx+1;
							if(lx2<=local_ncn) celij(j1,j1+local_ncn+1,kind_p1,ini_kind_p2,ly2);   //!N-East
						}
					}
				}
		}
		//		for (size_t j1 = 0;j1<store->nx;j1++)
		//			if ((kind_p1 == 1)?(store->nc_k1[j1]>0):(store->nc_k2[j1]>0))
		//				self(j1,kind_p1,ini_kind_p2);
		int local_ncn = store->local_ncn,local_ncm = store->local_ncm,local_ncl=store->local_ncl;
		int lzmax = store->local_ncl_interior[1],lymax = store->local_ncm_interior[1],lxmax = store->local_ncn_interior[1];
		for (int lz = store->local_ncl_interior[0]+1;lz<=lzmax;lz++)
			for (int ly = store->local_ncm_interior[0]+1;ly<=lymax;ly++)
				for (int lx = store->local_ncn_interior[0]+1;lx<=lxmax;lx++)
				{
					size_t j1 = (lx-1) + (ly-1)*local_ncn +(lz-1)*local_ncn*local_ncm;
					if ((kind_p1 == 1)?(store->nc_k1[j1]>0):(store->nc_k2[j1]>0))
						self(j1,kind_p1,ini_kind_p2);
				}
	}
	return;
}
void AC_MLS::self(int j1, int kind_p1, int ini_kind_p2)
{
	for (int kind_p2 = ini_kind_p2;kind_p2<=2;kind_p2++)
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
				if(i>=store->nbb && j>=store->nbb)
				{

					double drx = store->xp[i] - store->xp[j];
					double drz = store->zp[i] - store->zp[j];
					double dry = (store->dim == 3)?(store->yp[i] - store->yp[j]):(0);

					double rr2 = drx*drx + dry*dry + drz*drz;

					if(rr2<4*store->h*store->h && rr2>pow(10,-18))
					{

						//						c               Calculating kernel & Normalized Kernel

						solve->kernel_func->kernel(drx,dry,drz,i,j,j1,j1,rr2);

						double V_j = store->pm[j] / store->rho[j];
						double V_i = store->pm[i] / store->rho[i];
						if(i_MLS_part == 1)//   !1st sweep for MLS
						{
							//!--- i-Matrix ---
							store->A_MLS[i][0][0] += solve->kernel_func->Wab*V_j;
							store->A_MLS[i][0][1] += drx*solve->kernel_func->Wab*V_j;
							store->A_MLS[i][0][2] += dry*solve->kernel_func->Wab*V_j;
							store->A_MLS[i][0][3] += drz*solve->kernel_func->Wab*V_j;

							store->A_MLS[i][1][1] += drx*drx*solve->kernel_func->Wab*V_j;
							store->A_MLS[i][1][2] += drx*dry*solve->kernel_func->Wab*V_j;
							store->A_MLS[i][1][3] += drx*drz*solve->kernel_func->Wab*V_j;

							store->A_MLS[i][2][2] += dry*dry*solve->kernel_func->Wab*V_j;
							store->A_MLS[i][2][3] += dry*drz*solve->kernel_func->Wab*V_j;

							store->A_MLS[i][3][3] += drz*drz*solve->kernel_func->Wab*V_j;

							//!--- j-Matrix ---
							store->A_MLS[j][0][0] += solve->kernel_func->Wab*V_i;
							store->A_MLS[j][0][1] += drx*solve->kernel_func->Wab*V_i;
							store->A_MLS[j][0][2] += dry*solve->kernel_func->Wab*V_i;
							store->A_MLS[j][0][3] += drz*solve->kernel_func->Wab*V_i;

							store->A_MLS[j][1][1] += drx*drx*solve->kernel_func->Wab*V_i;
							store->A_MLS[j][1][2] += drx*dry*solve->kernel_func->Wab*V_i;
							store->A_MLS[j][1][3] += drx*drz*solve->kernel_func->Wab*V_i;

							store->A_MLS[j][2][2] += dry*dry*solve->kernel_func->Wab*V_i;
							store->A_MLS[j][2][3] += dry*drz*solve->kernel_func->Wab*V_i;

							store->A_MLS[j][3][3] += drz*drz*solve->kernel_func->Wab*V_i;
						}
						else //    !2nd Sweep for MLS
						{
							double Wab_MLSi = (store->beta0_MLS[i]+ store->beta1x_MLS[i]*drx+ store->beta1y_MLS[i]*dry+ store->beta1z_MLS[i]*drz)*solve->kernel_func->Wab;
							double Wab_MLSj = (store->beta0_MLS[j]+ store->beta1x_MLS[j]*drx+ store->beta1y_MLS[j]*dry+ store->beta1z_MLS[j]*drz)*solve->kernel_func->Wab;
							store->rho_sum_MLS[i] = store->rho_sum_MLS[i] +store->pm[j]*Wab_MLSi;
							store->rho_sum_MLS[j] = store->rho_sum_MLS[j] +store->pm[i]*Wab_MLSj;
						}

					}
				}
			}
		}
	}
}
void AC_MLS::celij(int j1, int j2,int kind_p1, int ini_kind_p2,int ylevel)
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
					if(i>=store->nbb && j>=store->nbb)
					{
						double drx = store->xp[i] - store->xp[j];
						double drz = store->zp[i] - store->zp[j];
						double dry = (store->dim == 3)?(store->yp[i] - store->yp[j]):(0);

						double rr2 = drx*drx + dry*dry + drz*drz;

						if(rr2<4*store->h*store->h && rr2>pow(10,-18))
						{
							//						c               Calculating kernel & Normalized Kernel

							solve->kernel_func->kernel(drx,dry,drz,i,j,j1,j1,rr2);

							double V_j = store->pm[j] / store->rho[j];
							double V_i = store->pm[i] / store->rho[i];
							if(i_MLS_part == 1)//   !1st sweep for MLS
							{
								//!--- i-Matrix ---
								store->A_MLS[i][0][0] += solve->kernel_func->Wab*V_j;
								store->A_MLS[i][0][1] += drx*solve->kernel_func->Wab*V_j;
								store->A_MLS[i][0][2] += dry*solve->kernel_func->Wab*V_j;
								store->A_MLS[i][0][3] += drz*solve->kernel_func->Wab*V_j;

								store->A_MLS[i][1][1] += drx*drx*solve->kernel_func->Wab*V_j;
								store->A_MLS[i][1][2] += drx*dry*solve->kernel_func->Wab*V_j;
								store->A_MLS[i][1][3] += drx*drz*solve->kernel_func->Wab*V_j;

								store->A_MLS[i][2][2] += dry*dry*solve->kernel_func->Wab*V_j;
								store->A_MLS[i][2][3] += dry*drz*solve->kernel_func->Wab*V_j;

								store->A_MLS[i][3][3] += drz*drz*solve->kernel_func->Wab*V_j;

								//!--- j-Matrix ---
								store->A_MLS[j][0][0] += solve->kernel_func->Wab*V_i;
								store->A_MLS[j][0][1] += drx*solve->kernel_func->Wab*V_i;
								store->A_MLS[j][0][2] += dry*solve->kernel_func->Wab*V_i;
								store->A_MLS[j][0][3] += drz*solve->kernel_func->Wab*V_i;

								store->A_MLS[j][1][1] += drx*drx*solve->kernel_func->Wab*V_i;
								store->A_MLS[j][1][2] += drx*dry*solve->kernel_func->Wab*V_i;
								store->A_MLS[j][1][3] += drx*drz*solve->kernel_func->Wab*V_i;

								store->A_MLS[j][2][2] += dry*dry*solve->kernel_func->Wab*V_i;
								store->A_MLS[j][2][3] += dry*drz*solve->kernel_func->Wab*V_i;

								store->A_MLS[j][3][3] += drz*drz*solve->kernel_func->Wab*V_i;
							}
							else //    !2nd Sweep for MLS
							{
								double Wab_MLSi = (store->beta0_MLS[i]+ store->beta1x_MLS[i]*drx+ store->beta1y_MLS[i]*dry+ store->beta1z_MLS[i]*drz)*solve->kernel_func->Wab;
								double Wab_MLSj = (store->beta0_MLS[j]+ store->beta1x_MLS[j]*drx+ store->beta1y_MLS[j]*dry+ store->beta1z_MLS[j]*drz)*solve->kernel_func->Wab;
								store->rho_sum_MLS[i] = store->rho_sum_MLS[i] +store->pm[j]*Wab_MLSi;
								store->rho_sum_MLS[j] = store->rho_sum_MLS[j] +store->pm[i]*Wab_MLSj;
							}
						}
					}
				}
			}
		}
	}
}
