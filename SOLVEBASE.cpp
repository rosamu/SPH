/*
 * SOLVEBASE.cpp
 *
 *  Created on: Jul 6, 2011
 *      Author: tony
 */

#include "SOLVEBASE.h"
#include "math.h"
#include "AC.h"
#include <iostream>

SOLVE_BASE::SOLVE_BASE(STORAGE *store,EOS_BASE *eos):store(store),equation_of_state(eos),itime(0),time(0),dt(store->dt),visc_dt(0){
	// TODO Auto-generated constructor stub
	store->nprint = 0;
	store->next_print_time = store->out;
	store->i_output = 0;
}

SOLVE_BASE::~SOLVE_BASE() {
	// TODO Auto-generated destructor stub
}

void SOLVE_BASE::ModifyDtForPrint()
{
	while (itime == 1 && store->next_print_time - time < store->dt)
		store->next_print_time += store->out;
	if (time + dt >= store->next_print_time && time+dt < store->tmax)
	{
		dt = store->next_print_time - time + pow(10,-16);
		store->next_print_time += store->out;
		if (store->next_print_time>store->tmax)
			store->next_print_time=store->tmax;
		store->i_output = 1;
		store->nprint++;
	}
	else if (time+dt>=store->tmax)
	{
		dt = store->tmax-time+pow(10,-16);
		store->i_output = 1;
		store->nprint++;
	}
}

void SOLVE_BASE::densityFilter()
{

	AC_BASE *ac;
	if (store->i_densityFilter == 2)
		ac = new AC_MLS(store,this);
	else
		assert(false);
	double A_input_matrix[4][4], A_LU_equivalent_matrix[4][4];
	double b_input_vector[4], x_solution[4];
	double A_MLS_inverse[4][4];


	vector<double>::size_type nstart_MLS = store->nb;   //!nstart
	//c     -- Zeroing Variables --
	vector<double>::size_type i;
	//	MYARRAY myarray;
	//	for (int j =0;j<4;j++)
	//		for (int k = 0;k<4;k++)
	//			myarray.data[j][k] = 0;
	//	A_MLS.assign(np,myarray);
	for (size_t i=0;i<store->np;i++)
	{
		for (int j =0;j<4;j++)
			for (int k = 0;k<4;k++)
				store->A_MLS[i][j][k] = 0;
		store->rho_sum_MLS[i] = 0;
		store->sum_wab[i] = 0;
		store->rho_sum[i] = 0;
		store->beta0_MLS[i] = 0;
		store->beta1x_MLS[i] = 0;
		store->beta1y_MLS[i] = 0;
		store->beta1z_MLS[i] = 0;
	}

	//	vector<double>::size_type i_check = nb;
	int itime_check = 10*10;   //!0     !10e8
	itime_check = itime_check*itime_check;
	itime_check = itime_check*itime_check;
	ac->i_MLS_part = 1;
	ac->doit();    //!Perform the correction matrix construction

	//    !construct LU-decomposition matrix
	int i_singularMatrix_count = 0;
	int length = 4;
	for (i=nstart_MLS;i<store->np;i++)
	{
		//      !temporary matrix
		//      !Self contribution
		store->A_MLS[i][0][0] = store->A_MLS[i][0][0]+store->adh*store->pm[i]/store->rho[i];  //!self contribution
		double A_MLS_11 = store->A_MLS[i][0][0];
		//!Account for A_MLS Symmetry
		store->A_MLS[i][1][0] = store->A_MLS[i][0][1];
		store->A_MLS[i][2][0] = store->A_MLS[i][0][2];
		store->A_MLS[i][3][0] = store->A_MLS[i][0][3];
		store->A_MLS[i][2][1] = store->A_MLS[i][1][2];
		store->A_MLS[i][3][1] = store->A_MLS[i][1][3];
		store->A_MLS[i][3][2] = store->A_MLS[i][2][3];


		for (int j=0;j<4;j++)
			for (int k=0;k<4;k++)
				A_input_matrix[j][k] = store->A_MLS[i][j][k];

		//      !Amatrix_temp(1,1) = temp  !sum_wab(i)
		int i_Singular;
		int i_A_MLS_step = 1;
		LU_decomposition(A_input_matrix,b_input_vector,length,A_LU_equivalent_matrix,i_A_MLS_step,i_Singular,x_solution,i);
		//      !i_Singular = 1
		if(i_Singular>0)
		{
			//        !print*,'singular MLS matrix!'
			//        !print*,'i',i
			//        !print*,'xp, yp, rhop ',xp(i),yp(i),rhop(i)
			//        !print*,'sum_wab(i)',sum_wab(i)
			//        !print*,'Amatrix_temp(1,1)',Amatrix_temp(1,1)
			//        !stop
			if(A_MLS_11 == 0.0)
			{
				cout<<"A_MLS_11 = 0"<<endl;
				assert(false);
			}
			i_singularMatrix_count = i_singularMatrix_count +1;
			store->beta0_MLS[i] = 1.0/A_MLS_11;  //!1.0/temp   !1.0/sum_wab(i)
			store->beta1x_MLS[i]= 0.0;
			store->beta1y_MLS[i]= 0.0;
			store->beta1z_MLS[i]= 0.0;
			//        !print*,'Singular i, xp, zp ',i,xp(i),zp(i)
		}
		else
		{
			//        !Find inverse matrix by solving system n times
			//        !on n successive unit vectors (Gerald & Wheatley, Numerical Analysis)
			i_A_MLS_step = 2;
			for (int k=0;k<4;k++)
			{
				for (int j=0;j<4;j++)
					b_input_vector[j] = 0.0;
				b_input_vector[k] = 1.0;
				LU_decomposition(A_input_matrix,b_input_vector,length,A_LU_equivalent_matrix,i_A_MLS_step,i_Singular,x_solution,i);
				for (int j=0;j<4;j++)
					A_MLS_inverse[j][k] = x_solution[j];
			}
			store->beta0_MLS[i] = A_MLS_inverse[0][0];
			store->beta1x_MLS[i]= A_MLS_inverse[1][0];
			store->beta1y_MLS[i]= A_MLS_inverse[2][0];
			store->beta1z_MLS[i]= A_MLS_inverse[3][0];

		}

	}


	//c      print*,'itime, MLS i_singularMatrix_count ',
	//c     &        itime, i_singularMatrix_count
	//	sum_wab.assign(np,0);
	//	rho_sum.assign(np,0);
	for (size_t i=0;i<store->np;i++)
	{
		store->sum_wab[i] = 0;
		store->rho_sum[i] = 0;
	}
	ac->i_MLS_part = 2;

	ac->doit();  // !Perform the corrected density summation
	for (i = nstart_MLS;i<store->np;i++)
	{
		store->rho_sum_MLS[i] = store->rho_sum_MLS[i] + store->beta0_MLS[i]*store->adh*store->pm[i];  //!self contribution
		store->rho[i] = store->rho_sum_MLS[i];
		equation_of_state->eos(store->rho[i],store->TE[i],store->p[i],store->cs[i]);
	}
	delete ac;
	return;
}

void SOLVE_BASE::LU_decomposition(double A_input_matrix[4][4],double b_input_vector[4],
		int length,double A_LU_equivalent_matrix[4][4],int i_A_MLS_step,
		int &i_Singular,double x_solution[4],int i)
{
	//    !Solves Ax = b using LU Decomposition
	//    !A = LU  ==> LUx = b
	//    !Ux = y  ==> Ly  = b
	//    !Solve Ly = b using forward subsitution
	//    !Solve Ux = y using backward substitution
	double A_temp[4][4];
	double L_matrix[4][4];
	double U_matrix[4][4];
	double y_solution[4];
	double al_min,al_max;
	double tol=pow(10.0,-16);

	if(i_A_MLS_step == 1)//  !Construct LU Matrix
	{

		//     !- Check singularity -
		double amax_Aij = 0.0;
		double row_sum_min = pow(10,8);
		i_Singular = 0;
		for (int iii=0;iii<4;iii++)
		{
			double row_sum = 0.0;
			for (int jjj=0;jjj<4;jjj++)
			{
				//         !Transfer input matrix into a temporary matrix
				A_temp[iii][jjj] = A_input_matrix[iii][jjj];
				//         !Check matrix is not singular
				double A_temp_ij_abs = fabs(A_temp[iii][jjj]);
				row_sum = row_sum + A_temp_ij_abs;
				if(A_temp_ij_abs>amax_Aij)
					amax_Aij = A_temp_ij_abs;
			}
			if(row_sum<row_sum_min)
			{
				row_sum_min = row_sum;
			}

		}
		if(amax_Aij < tol) i_Singular = 1;
		if(row_sum_min < tol) i_Singular = 1;

		if(i_Singular<1)
		{
			for (int iii=0;iii<4;iii++)
				L_matrix[iii][0]=A_input_matrix[iii][0];
			for (int jjj=0;jjj<4;jjj++)
				U_matrix[0][jjj]=A_input_matrix[0][jjj] / L_matrix[0][0];

			//       !- Forward elimination using Gaussian Elimination for matrices L and U -
			for (int jjj=1;jjj<4;jjj++)
			{
				for (int iii=jjj;iii<4;iii++)
				{
					double sum_temp = 0.0;
					for(int kkk = 0;kkk<=jjj-1;kkk++)
						sum_temp = sum_temp + L_matrix[iii][kkk]*U_matrix[kkk][jjj];
					L_matrix[iii][jjj] = A_input_matrix[iii][jjj] - sum_temp;
				}
				U_matrix[jjj][jjj] = 1.0;
				for(int iii = jjj+1;iii<4;iii++)
				{
					double sum_temp= 0.0;
					for (int kkk = 0;kkk<=jjj-1;kkk++)
						sum_temp = sum_temp + L_matrix[iii][kkk]*U_matrix[kkk][jjj];
					U_matrix[jjj][iii] = (A_input_matrix[jjj][iii] - sum_temp)/L_matrix[jjj][jjj];
				}
			}
			//       !-- Suggestion by Arno Mayrhofer for singular ill-conditioned L matrices --
			al_min = pow(10.0,8);
			al_max = -pow(10.0,8);
			for (int iii=0;iii<4;iii++)
			{
				if(L_matrix[iii][iii]>al_max) al_max = L_matrix[iii][iii];
				if(L_matrix[iii][iii]<al_min) al_min = L_matrix[iii][iii];
			}
			if(fabs(al_min) < pow(100,-14)*fabs(al_max)) i_Singular = 1;
			//       !--  End of Arno Suggestion --


			//       !- Combine into one matrix -
			for(int iii=0;iii<4;iii++)  //!Row
			{
				for(int jjj = iii+1;jjj<4;jjj++) //!Column
				{
					A_LU_equivalent_matrix[iii][jjj] = U_matrix[iii][jjj];
				}
				for (int jjj = 0;jjj<=iii;jjj++) //!Column
				{
					A_LU_equivalent_matrix[iii][jjj] = L_matrix[iii][jjj];
				}
			}

		}
	}
	else  //!if(i_MLS_step.eq.2)then  !Solver equation
	{
		for (int iii=0;iii<4;iii++)
		{
			x_solution[iii] = 0.0;//  !zero initialise solution vector
			y_solution[iii] = 0.0;//  !zero initialise temp     vector
		}
		//    !- Forward Substitution using Lower Matrix -
		y_solution[0]=b_input_vector[0]/A_LU_equivalent_matrix[0][0]; //!First equation has only one unknown
		for (int iii = 1;iii<4;iii++)
		{
			double sum_Lijyj = 0.0;
			//      !- Generating the summation term -
			for (int jjj = 0;jjj<=iii-1;jjj++)
				sum_Lijyj = sum_Lijyj+ A_LU_equivalent_matrix[iii][jjj]*y_solution[jjj];
			y_solution[iii]=(b_input_vector[iii]- sum_Lijyj) /A_LU_equivalent_matrix[iii][iii];;
		}

		//    !- Back Substitution -
		x_solution[3] = y_solution[3]; //!Nth equation has only one unknown
		for (int iii = 2;iii>=0;iii--)
		{
			double sum_Uijxj = 0.0;
			//      !- Generating the summation term -
			for (int jjj=iii+1;jjj<4;jjj++)
				sum_Uijxj = sum_Uijxj+ A_LU_equivalent_matrix[iii][jjj]*x_solution[jjj];
			x_solution[iii] = (y_solution[iii] - sum_Uijxj);
		}

	}
}

void SOLVE_BASE::computeElectricFieldPotential()
{
//	//initialization
//    PETSc petsc;
//	store->m_pSolver = &petsc;
//	store->m_pSolver->Create(store->nbb, store->np-1, 24, 0);
//	for (size_t i=0;i<store->np;i++)
//		store->Phi[i] = store->Jx[i] = store->Jy[i] = store->Jz[i] = store->GPx[i] = store->GPy[i] = store->GPz[i] = store->dUxB[i] = 0;
//	AC_BASE *ac;
//	ac = new AC_ELECTRIC(store,this);
//	double sigma=1e13,lightSpeed=3e7;
//	//step 1: calculate divergence UxB
//	ac->i_ELE_part = 1;
//	ac->doit();
//	for (size_t i=0;i<store->np;i++)
//		store->m_pSolver->Add_b(i,store->dUxB[i+store->nbb]/lightSpeed);
//	//step 2: Calculate PHi
//
//	store->m_pSolver->SetTol(1e-15);
//	ac->i_ELE_part = 2;
//	ac->doit();
//	store->m_pSolver->Print_A("A.txt");
//	store->m_pSolver->Print_b("b.txt");
//
//	//	m_pSolver->Solve();
//	store->m_pSolver->Solve_withPureNeumann();
//
//	int iter;
//	double residual;
//	store->m_pSolver->GetNumIterations(&iter);
//	store->m_pSolver->GetFinalRelativeResidualNorm(&residual);
//
//	double *x = new double[store->np-store->nbb];
//	store->m_pSolver->Get_x(x);
//
//	//debug_saveCompList("comp0");
//	for (size_t i = store->nbb;i<store->np;i++)
//		store->Phi[i]= x[i-store->nbb];
//	delete [] x;
//	//step 3: Calculate J
//	ac->i_ELE_part = 3;
//	ac->doit();
//	for (size_t i=store->nbb;i<store->np;i++)
//	{
//		double UxB[3];
//		store->GetUxB(UxB,i);
//		store->Jx[i] = sigma*(-store->GPx[i]+1/lightSpeed*UxB[0]);
//		store->Jy[i] = sigma*(-store->GPy[i]+1/lightSpeed*UxB[1]);
//		store->Jz[i] = sigma*(-store->GPz[i]+1/lightSpeed*UxB[2]);
//	}
}

void SOLVE_BASE::TimeSync()
{
	  double * dt_array;
	  int gsize =0;
	  int root = 0;
	  MPI_Barrier(MPI_COMM_WORLD);
	  if (store->id == root)
	    {
	      MPI_Comm_size(MPI_COMM_WORLD,&gsize);
	      dt_array = new double[gsize];
	    }
	  MPI_Gather(&dt, 1,MPI_DOUBLE, dt_array, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
	  MPI_Barrier(MPI_COMM_WORLD);
	  if (store->id == root)
	    {
	      std::vector<double> myvector(dt_array,dt_array+gsize);
	      std::sort(myvector.begin(),myvector.end());
	      dt = myvector[0];
	    }
	  MPI_Barrier(MPI_COMM_WORLD);
	  if (store->id == root)
	    {
	      MPI_Bcast(&dt, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
	    }
	  else
	    {
	      MPI_Bcast(&dt, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
	    }
	  MPI_Barrier(MPI_COMM_WORLD);
}
