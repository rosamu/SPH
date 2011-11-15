//============================================================================
// Name        : SPH.cpp
// Author      : Tongfei Guo
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>
#include "SPH.h"
#include <string>
#include "STORAGE.h"
#include "SOLVE.h"
#include "EOS.h"
#include "KERNEL.h"
#include "mpi.h"
using namespace std;

int main(int argc, char **argv) {
	STORAGE store;
	SOLVE_BASE *solve;
	EOS_BASE *equation_of_state;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&store.NumOfCpu);
	MPI_Comm_rank(MPI_COMM_WORLD,&store.id);
	//	MPI::Init ( argc, argv );
	//
	//	store.id = MPI::COMM_WORLD.Get_rank ( );
	//
	//	store.NumOfCpu = MPI::COMM_WORLD.Get_size ( );
	istringstream line;
	if (argv[1][1] == 'i')
	{
		string input(argv[2]);
		store.PN = 1;
		store.PM = 1;
		store.PL = 1;
		store.Initialization(input);
	}
	else
	{
		if (argv[2][0] == '3')
		{
			//store.PN = argv[4][0] - 48;
			//store.PM = argv[5][0] - 48;
			//store.PL = argv[6][0] - 48;
			string input = argv[4];
			line.clear();
			line.str(input);
			line>>store.PN;
			line.clear();
			input = argv[5];
			line.clear();
			line.str(input);
			line>>store.PM;
			line.clear();
			store.PL = argv[6][0] - 48;
			input=argv[8];
			store.Initialization(input);
		}
		else if (argv[2][0] == '2')
		{
			store.PN = argv[4][0] - 48;
			store.PM = 1;
			store.PL = argv[5][0] - 48;
			string input(argv[7]);
			store.Initialization(input);
		}
		else
		{
			store.PN = argv[4][0] - 48;
			store.PM = 1;
			store.PL = 1;
			string input(argv[6]);
			store.Initialization(input);
		}
	}
	//PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
	MPI_Barrier(MPI_COMM_WORLD);
	if (store.i_EoS == 4)
		equation_of_state = new EOS_SPOLY(&store);
	else if (store.i_EoS == 1)
		equation_of_state = new EOS_TAIT(&store);
	else
	{
		assert(0);
	}
	if (store.i_algorithm == 1)
		solve = new SOLVE_PREDICTOR_CORRECTOR(&store,equation_of_state);
	if (store.i_kernel == 3)
		solve->kernel_func=new KERNEL_CUBIC(&store);
	else if (store.i_kernel == 5)
		solve->kernel_func=new KERNEL_WENDLAND5(&store);
	else
	{solve->kernel_func = NULL;}
	if (store.i_kernelcorrection == 0)
		solve->kernel_correction_func = new KERNEL_CORRECTION_NC(&store);
	else
	{assert(0);}
	solve->initialization();
	string output;
	if (argv[3][1] == 'o')
		output = argv[4];
	else
		output = argv[10];
	char filename[200];
	sprintf(filename,"%s",output.c_str());
	if (store.PN * store.PM * store.PL != 1)
		sprintf(filename,"%s-nd%s",output.c_str(),store.right_flush(store.id,4));
	ofstream myfile(filename);
	double wtime = 0;
	double wtime_init = 0;
	if (store.id == 0)
	{
		wtime_init = MPI_Wtime();
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (store.i_restartRun == 1)
	{
		store.interior_read(solve->itime,solve->dt,solve->time);
		MPI_Barrier(MPI_COMM_WORLD);
	}
	store.divide(0,store.nbb1,1);
	store.i_SR_step = 0;
	store.SendRecvBuffer();
	solve->parallelinitialization();
	store.detail = 0;
	if (store.i_restartRun != 1)
	{
		store.print_out(0,0,output.c_str());
		store.StatesPrint(0,0,"Density",output.c_str());
		store.StatesPrint(0,0,"Pressure",output.c_str());
		store.StatesPrint(0,0,"Interface",output.c_str());
	}
	MPI_Barrier(MPI_COMM_WORLD);
	assert(0);
	while (solve->time < store.tmax && solve->itime < 40000)
	{
		solve->ModifyDtForPrint();
		solve->time+=solve->dt;
		solve->itime++;
		store.t = solve->time;
		solve->step();
		wtime = MPI_Wtime ( ) - wtime_init;
		myfile<<" Time is "<<solve->time<<" Step is "<<solve->itime<<" Next dt is "<< solve->dt<<" Actual Time is"<<wtime<<endl;
		if (wtime > store.interior_print_out_time_interval)
		{
			store.interior_print_out_time_interval += 1800;
			store.interior_print_out(solve->itime,solve->dt,solve->time);
		}
		if (store.i_output == 1)
		{
			store.i_output = 0;
			store.print_out(store.nprint,solve->time,output.c_str());
			store.StatesPrint(store.nprint,solve->time,"Density",output.c_str());
			store.StatesPrint(store.nprint,solve->time,"Pressure",output.c_str());
			store.StatesPrint(store.nprint,solve->time,"Interface",output.c_str());
		}
	}
	if (solve->time < store.tmax)
	{
		store.print_out(++store.nprint,solve->time,output.c_str());
		store.StatesPrint(store.nprint,solve->time,"Density",output.c_str());
		store.StatesPrint(store.nprint,solve->time,"Pressure",output.c_str());
		store.StatesPrint(store.nprint,solve->time,"Interface",output.c_str());
	}
	if (store.id == 0)
	{
		wtime = MPI_Wtime ( ) - wtime_init;
		myfile<< "  Wall clock elapsed seconds = " << wtime << endl;
		cout << "\n";
		cout << "  Wall clock elapsed seconds = " << wtime << endl;
	}
	myfile.close();
	MPI_Finalize();
	return 0;
}
