/*
 * SOLVEBASE.h
 *
 *  Created on: Jul 6, 2011
 *      Author: tony
 */

#ifndef SOLVEBASE_H_
#define SOLVEBASE_H_
#include "STORAGE.h"
#include "EOSBASE.h"
#include "KERNEL.h"

class SOLVE_BASE {
public:
	SOLVE_BASE(STORAGE*,EOS_BASE*);
	virtual ~SOLVE_BASE();
	//methods
	virtual void step()=0;
	virtual void initialization()=0;
	virtual void parallelinitialization()=0;
	void ModifyDtForPrint();
	void densityFilter();
	void LU_decomposition(double A_input_matrix[4][4],double b_input_vector[4],int length,double A_LU_equivalent_matrix[4][4],int i_A_MLS_step,int &i_Singular,double x_solution[4],int i);
	void computeElectricFieldPotential();
	virtual void allocation()=0;
	virtual void dellocation()=0;
	void TimeSync();
	//data
	STORAGE *store;
	EOS_BASE *equation_of_state;
	KERNEL_BASE *kernel_func;
	KERNEL_CORRECTION_BASE *kernel_correction_func;
	int itime;
	double time;
	double dt;
	double visc_dt;


};

#endif /* SOLVEBASE_H_ */
