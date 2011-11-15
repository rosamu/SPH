/**********************************************************************
 * 		solver.h
 **********************************************************************/

#ifndef _FT_IFLUID_SOLVER_H_
#define _FT_IFLUID_SOLVER_H_

#include <vector>
#include <petscksp.h>
//#include <petscmg.h>
#include <petscpcmg.h>
#include <petscpc.h>
#include <assert.h>

class SOLVER
{
public:
	SOLVER(){};
	SOLVER(int ilower, int iupper, int d_nz, int o_nz){};
	virtual ~SOLVER(){};
	virtual void Create(int ilower, int iupper, int d_nz, int o_nz){};
	
	virtual void Set_A(PetscInt i, PetscInt j, double val){};	// A[i][j]=val;
	virtual void Add_A(PetscInt i, PetscInt j, double val){};	// A[i][j]=A[i][j]+val;
	virtual void Set_x(PetscInt i, double val){};	// x[i]=val;
	virtual void Set_x(double *p){};		// x[i]=p[i];
	virtual void Add_x(PetscInt i, double val){};	// x[i]=x[i]+val;
	virtual void Get_x(double *p){};	// get the x from ij_x to p.		
	virtual void Get_x(double *p, int n, int *global_index){};
	virtual void Set_b(PetscInt i, double val){};	// b[i]=val;
	virtual void Set_b(double *b){};	
	virtual void Add_b(PetscInt i, double val){};	// b[i]=b[i]+val;

	virtual void SetMaxIter(int val){};	
	virtual void GetFinalRelativeResidualNorm(double *rel_resid_norm){};
	virtual void GetNumIterations(int *num_iterations){};

	virtual void Solve(void){};	
	virtual void Solve_withPureNeumann(void){};	
	virtual void Read_A(char *filename){};
	virtual void Print_A(char *filename){};
	virtual void Read_b(char *filename){};
	virtual void Print_b(char *filename){};
	virtual void Read_x(char *filename){};
	virtual void Print_x(char *filename){};
	virtual void test(void){};
};

class PETSc: public SOLVER
{
public:	
	MPI_Comm  comm;			// set to be MPI_COMM_WORLD.
	int iLower;
	int iUpper;			// global row range
	
	Vec x;      			/* approx solution, RHS*/
	Vec b;
  	Mat A;          		/* linear system matrix */
  	
  	KSP   ksp;          		/* Krylov subspace method context */
	PC    pc;
	MatNullSpace	nullsp;
		
	PetscErrorCode ierr;
	int its;			// numer of iterations;

public:
	PETSc();
	PETSc(int ilower, int iupper, int d_nz, int o_nz);		
		// global row range of A, x, b on this processor
	~PETSc();
	void Create(int ilower, int iupper, int d_nz, int o_nz);	
		// same as Hypre(int, int)
	void Create(MPI_Comm Comm, int ilower, int iupper, int d_nz, int o_nz);
		// same as Hypre(int, int)
	
	void Reset_A();				// Set A[i][j]=0.0;
	void Reset_b();
	void Reset_x();
	void Set_A(PetscInt i, PetscInt j, double val);	// A[i][j]=val;
	void Add_A(PetscInt i, PetscInt j, double val);	// A[i][j]=A[i][j]+val;
	void Get_row_of_A(PetscInt i, PetscInt *ncol, PetscInt **cols, double **row);
	void Set_x(PetscInt i, double val);		// x[i]=val;
	void Add_x(PetscInt i, double val);		// x[i]=x[i]+val;
	void Set_b(PetscInt i, double val);		// b[i]=val;
	void Add_b(PetscInt i, double val);		// b[i]=b[i]+val;
	void Get_x(double *p);		// get the x from ij_x to p.	
	void Get_b(double *p);		// get the b from ij_x to p.
	void Get_x(double *p, int n, int *global_index);
	
	void SetMaxIter(int val); 	// Set maximum number of iterations 
	void SetTol(double val);	// Set the convergence tolerance 
	void SetKDim(int k_dim);	
			// Set the maximum size of the Krylov space 
	void GetNumIterations(PetscInt *num_iterations);	
			// Return the number of iterations taken 
	void GetFinalRelativeResidualNorm(double *rel_resid_norm);
	void Solve(void);
	void Solve_GMRES(void);
	void Solve_withPureNeumann(void);
	void Solve_withPureNeumann_GMRES(void);
	virtual void Print_A(const char *filename);
        virtual void Print_b(const char *filename);
};
#endif
