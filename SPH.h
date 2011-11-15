/*
 * SPH.h
 *
 *  Created on: Jul 6, 2011
 *      Author: tony
 */

#ifndef SPH_H_
#define SPH_H_
#include <vector>
using namespace std;
struct _FIELD
{
	vector<double> xp;
	vector<double> yp;
	vector<double> zp;
	vector<double> up;
	vector<double> vp;
	vector<double> wp;
	vector<double> p;
	vector<double> rho;
	vector<double> TE;
	vector<double> pm;
	vector<double> cs;
	vector<int> iflag;
};
typedef struct _FIELD FIELD;
struct _SWEEP
{
		vector<double> xp;
		vector<double> yp;
		vector<double> zp;
		vector<double> up;
		vector<double> vp;
		vector<double> wp;
		vector<double> rho;
		vector<double> TE;
		vector<double> p;
};
typedef struct _SWEEP SWEEP;
struct _SWEEPDOT
{
		vector<double> xpdot;
		vector<double> ypdot;
		vector<double> zpdot;
		vector<double> updot;
		vector<double> vpdot;
		vector<double> wpdot;
		vector<double> rhodot;
		vector<double> TEdot;
};
typedef struct _SWEEPDOT SWEEPDOT;

struct _BOX
{
	int nc[2];
	int nc_bk[2];
	vector<vector<double>::size_type> ibox_k1;
	vector<vector<double>::size_type> ibox_k2;
};

#define index3d(i,j,k,imax,jmax,kmax)  (i + (j-1)*imax +(k-1)*jmax*imax) //i,j,k start from 1
#define index2d(i,j,imax,jmax) (i + (j-1)*imax) //i,j start from 1
typedef _BOX BOX;
#endif /* SPH_H_ */
