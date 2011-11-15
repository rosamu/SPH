/*
 * STORAGE.h
 *
 *  Created on: Jul 6, 2011
 *      Author: tony
 */

#ifndef STORAGE_H_
#define STORAGE_H_
#include <string>
#include <fstream>
#include <sstream>
#include "SPH.h"
#include "assert.h"
#include <list>
//#include "solver.h"
#include "mpi.h"
#include <algorithm>
#include <stdio.h>
using namespace std;
#define Index3D(i,j,k,n,m,l) ((k)*(m)*(n)+(j)*(n)+(i))

class STORAGE {
public:
	STORAGE();
	virtual ~STORAGE();

	//methods
	void check_limits();
	void Boundaries_FrontBack(int,int,int,double,double,double);
	void Boundaries_LeftRight(int,int,int,double,double,double);
	void Boundaries_Bottom(int,int,int,double,double,double);
	void Fill_Part(int,int,int,double,double,double,double,double,double);
	void Drop(double,double,double,double,double,double,double,double,double,double,double);
	void pos_veloc(double,double,double,double,double,double);
	void pressure(double,double,double);
	void pressure(double,double,double,double);
	void pressure_b(double,double,double);
	void Initialization(string);
	void divide(int,int,int);
	void print_divide();
	const char *right_flush(int,int);
	void print_out(int,double,const char*);
	void print_out_cell(int,double);
	void gradients_calc(int i,int j,double dux,double duy,double duz);
	void gradients_calc_conservative(int i,int j,double dux,double duy,double duz){assert(false);}
	void viscosity(double dot,double drx,double dry,double drz,double dux,double duy,double duz,double rr2,double cbar,double robar,double one_over_rhobar,int i,int j,int j1,int j2);
	void allocation();
	void GetB(double,double,double,double,double []);
	void GetUxB(double [],size_t);
	void CrossProduct(double a[],double b[],double c[]);//calculate axb=c
	void ParallelInitialization();
	void ParallelAdjustment();
	void SendRecvBuffer();
	void ParallelAllocation(size_t);
	void ParticleNumberSync();
	void InterfacePrint(int,double);
	void StatesPrint(int,double,char*,const char*);
	int exists(const char *fname);
	void particlegathering(int step);
	void interior_print_out(int,double,double);
	void interior_read(int&,double&,double&);
	//data

	///initial data
	int dim;
	int i_kernel;
	int i_algorithm;
	int i_densityFilter;  // Using a density filter
	int i_viscos;
	int IBC; //BC
    int i_periodicOBs[3];
    int lattice;
    int i_EoS;
    double h_SWL;
    double B;
    double gamma;
    double pinf;
    double einf;
    double p0;
    double coef;
    double eps;   // epsilon XSPH
    double rho0;
    double viscos_val;
    double visc_wall;
    double vlx[2];    //Dimensions of container box
    double vly[2];
    double vlz[2];
    double dx;
    double dy;
    double dz;
    double h;
    size_t np;
    size_t nb;
    size_t nbb1;
    size_t nbb;
    unsigned nx;
    int ivar_dt;
    double dt;
    double tmax;     // End of the run(seconds)
    double out;      // output every out seconds
    double trec_ini; // Initial output
    double dtrec_det; //Detailed recording
    double t_sta_det;
    double t_end_det;
    double interior_print_out_time_interval;
    int i_restartRun;
    double CFL_number;
    double coefficient;
    double TE0;
    int i_kernelcorrection;
    int iRiemannSolver;
    int iTVD;
    double beta_lim;
    int i_vort;
    int ndt_VerletPerform;
    int ndt_FilterPerform;
    int ndt_DBCPerform;
    int mv_xmin;
    int mv_xmax;
    int mv_ymin;
    int mv_ymax;
    int mv_zmin;
    int mv_zmax;
    double g;
    double vnorm_mass;
    double XXmin,XXmax,YYmin,YYmax,ZZmin,ZZmax;
    int N,M,L;
    int ncn,ncm,ncl;
    int i_correct_pb;
    double zmax;
    double grx,gry,grz;
    int i_output;
    double next_print_time;
    int nprint;
    vector<double>::size_type nn;
    ///operating data
    FIELD fld;
    SWEEP swp_old, swp_new;
    SWEEPDOT swpdot;
//    vector<BOX> cells;
    vector<size_t> nc_k1;
    vector<size_t> nc_k2;
    vector< vector<size_t> > ibox_k1;
    vector< vector<size_t> > ibox_k2;
    double t;
    ///for ac_post
    double bigUdot;
    double bigVdot;
    double bigWdot;
    double X_Friction;
    double Y_Friction;
    double Z_Friction;
    int nb_inFriction;
//	vector<double>	         dudx_CSPH;
//	vector<double>	         dudy_CSPH;
//	vector<double>	         dudz_CSPH;
//	vector<double>	         dvdx_CSPH;
//	vector<double>	         dvdy_CSPH;
//	vector<double>	         dvdz_CSPH;
//	vector<double>	         dwdx_CSPH;
//	vector<double>	         dwdy_CSPH;
//	vector<double>	         dwdz_CSPH;
//	vector<double>	         dTEdx_CSPH;
//	vector<double>	         dTEdy_CSPH;
//	vector<double>	         dTEdz_CSPH;
//	vector<double>	         drhodx_CSPH;
//	vector<double>	         drhody_CSPH;
//	vector<double>	         drhodz_CSPH;
//	vector<double> dudx_CSPHo;
//	vector<double> dudy_CSPHo;
//	vector<double> dudz_CSPHo;
//	vector<double> dvdx_CSPHo;
//	vector<double> dvdy_CSPHo;
//	vector<double> dvdz_CSPHo;
//	vector<double> dwdx_CSPHo;
//	vector<double> dwdy_CSPHo;
//	vector<double> dwdz_CSPHo;
//	vector<double> dTEdx_CSPHo;
//	vector<double> dTEdy_CSPHo;
//	vector<double> dTEdz_CSPHo;
//	vector<double> drhodx_CSPHo;
//	vector<double> drhody_CSPHo;
//	vector<double> drhodz_CSPHo;
	double * ax,* ay,* az,* ar,* ux,* vx,* wx,* aTE;
	double *xcor,*ycor,*zcor;
	double term2i;
	double term2j;
	//for kernel and kernel_correction
	double frx,fry,frz;
	double frxi,fryi,frzi,frxj,fryj,frzj;
	double index_tensile;
	///densityfilter
	double adh;
	double ***A_MLS;
    double i_max;
    int i_print_divide;
    double * xp,* yp,* zp,* up,* vp,* wp,* TE,* rho,* pm,* cs,* p, * meow, * meow2;
    double * xpo,* ypo,* zpo,* upo,* vpo,* wpo,* TEo,* rhoo,* po;
    double * xpdot,* ypdot,* zpdot,* updot,* vpdot,* wpdot,* TEdot,* rhodot,* pmdot,* csdot,* pdot;
    int * iflag;
    int * bflag;
	double * rho_sum_MLS,* sum_wab,* rho_sum,* beta0_MLS,* beta1x_MLS,* beta1y_MLS,* beta1z_MLS;
	int detail;
	//ELECTRONIC
//	PETSc *m_pSolver;
	double * Phi;
	double * Jx, * Jy, * Jz;
	double * GPx, * GPy, * GPz;
	double * dUxB;
	int i_Print_J;
	//Parallelization Information
	int PN, PM, PL; //Partition Dimension, read from command line;
	int local_ncn, local_ncm, local_ncl;
	int local_ncn_interior[2], local_ncm_interior[2], local_ncl_interior[2];
	int PartitionIndex[3];
    double localvlx[2];
    double localvly[2];
    double localvlz[2];
    int NumOfCpu;
    int id;
    double localXXmax, localXXmin;
    double localYYmax, localYYmin;
    double localZZmax, localZZmin;
    double xmin,ymin,zmin;
    vector<size_t> DisplacedAndBufferParticle;
    vector<double> DisplacedAndBufferParticlexp;
    vector<double> DisplacedAndBufferParticleyp;
    vector<double> DisplacedAndBufferParticlezp;
    vector<double> DisplacedAndBufferParticleup;
    vector<double> DisplacedAndBufferParticlevp;
    vector<double> DisplacedAndBufferParticlewp;
    vector<double> DisplacedAndBufferParticleTE;
    vector<double> DisplacedAndBufferParticlerho;
    vector<double> DisplacedAndBufferParticlepm;
    double * SendBuffer;
    double * RecvBuffer;
    size_t TotalNumberOfParticleSend;
    size_t TotalNumberOfParticleRecv;
    int i_SR_step;
    size_t np_max;
    double *scalar;
    int i_extendz;
};


#endif /* STORAGE_H_ */
