/*
 * STORAGE.cpp
 *
 *  Created on: Jul 6, 2011
 *      Author: tony
 */

#include "STORAGE.h"
#include "assert.h"
#include "math.h"
#include <iostream>
using namespace std;

STORAGE::STORAGE():i_print_divide(0),i_Print_J(0),i_extendz(0),interior_print_out_time_interval(1800) {
}

STORAGE::~STORAGE() {
	// TODO Auto-generated destructor stub
}

void STORAGE::check_limits()
{
	int ncases=0;
	vector<double>::size_type i;
	for (i = nb;i<np;i++)
	{
		if(zp[i] > zmax && iflag[i]==1)
		{
			zmax=zp[i];
			i_max = i;
		}
		if(1)//iflag[i]==1) should be modified later a little bit more and more
		{
			if (dim == 3)
			{
				if(xp[i] > vlx[1] || xp[i] < vlx[0] || yp[i] > vly[1] || yp[i]<vly[0] || zp[i]<vlz[0])
				{
					if (xp[i] > vlx[1])
					{
						xp[i] = vlx[1]-(xp[i]-vlx[1]);
						up[i] = -up[i];
					}
					else if (xp[i] < vlx[0])
					{
						xp[i] = vlx[0] + (vlx[0]- xp[i]);
						up[i] = - up[i];
					}
					else if (yp[i] > vly[1])
					{
						yp[i] = vly[1]-(yp[i]-vly[1]);
						up[i] = -vp[i];
					}
					else if (yp[i]<vly[0])
					{
						yp[i] = vly[0] + (vly[0]- yp[i]);
						up[i] = -vp[i];
					}
					else
					{
						zp[i] = vlz[0] + (vlz[0]- zp[i]);
						up[i] = -wp[i];
					}
					ncases=ncases+1;
				}
				if (zp[i] > zmax)
					i_extendz = 1;
			}
			else
			{
				if(xp[i] > N*dx || xp[i] < 0 || zp[i]<0)
				{
					if (xp[i] > N*dx)
					{
						xp[i] = N*dx-(xp[i]-N*dx);
						up[i] = -up[i];
					}
					else if (xp[i] < 0)
					{
						xp[i] = - xp[i];
						up[i] = -up[i];
					}
					else
					{
						zp[i] = -zp[i];
						up[i] = -wp[i];
					}
					ncases=ncases+1;
				}
			}
		}
	}
	return;
}

void STORAGE::Boundaries_FrontBack(int N, int M, int L, double dx, double dy, double dz)
{
	if (PartitionIndex[1] != 0 && PartitionIndex[1] != PM - 1)
		return;
	//	- Inner Layer -
	int n0,n1,l0,l1;
	if (localvlx[0] - vlx[0]< 0.5*dx)
		n0 = 0;
	else
		n0 = (int)((localvlx[0]-vlx[0]-0.499999999*dx)/dx)+1;
	n1 = (int)((min(vlx[1],localvlx[1])-vlx[0]-0.499999999*dx)/dx)+1;
	if (localvlz[0] - vlz[0]< 0.5*dz)
		l0 = 0;
	else
		l0 = (int)((localvlz[0]-vlz[0]-0.499999999*dz)/dz)+1;
	l1 = (int)((min(vlz[1],localvlz[1])-vlz[0]-0.499999999*dz)/dz)+1;
	for (int i = n0;i<n1;i++)
		for (int k = l0;k<l1;k++)
		{
			if (PartitionIndex[1] == 0)
			{
				nn++;
				pos_veloc(vlx[0]+(i+0.5)*dx, vly[0],vlz[0]+(k+0.5)*dz,0,0,0);
				pressure_b(dx,dy,dz);
			}
			if (PartitionIndex[1] == PM - 1)
			{
				nn++;
				pos_veloc(vlx[0]+(i+0.5)*dx, vly[1],vlz[0]+(k+0.5)*dz,0,0,0);
				pressure_b(dx,dy,dz);
			}
		}

	if (localvlx[0] - vlx[0]<= 0)
		n0 = 0;
	else
		n0 = (int)((localvlx[0]-vlx[0])/dx) + 1;
	n1 = (int)((min(vlx[1],localvlx[1])+(1e-10)*dx-vlx[0])/dx);
	if (localvlz[0] - vlz[0]<= 0)
		l0 = 0;
	else
		l0 = (int)((localvlz[0]-vlz[0])/dz) + 1;
	l1 = (int)((min(vlz[1],localvlz[1])+(1e-10)*dz-vlz[0])/dz);
	//  - Outer Layer -
	for (int i = n0;i<=n1;i++)
		for (int k = l0;k<=l1;k++)
		{
			if (PartitionIndex[1] == 0)
			{
				nn++;
				pos_veloc(vlx[0]+i*dx,   vly[0]-0.5*dy,vlz[0]+k*dz,0,0,0);
				pressure_b(dx,dy,dz);
			}
			if (PartitionIndex[1] == PM - 1)
			{
				nn++;
				pos_veloc(vlx[0]+i*dx,   vly[1]+0.5*dy,vlz[0]+k*dz,0,0,0);
				pressure_b(dx,dy,dz);
			}
		}
}

void STORAGE::Boundaries_LeftRight(int N, int M, int L, double dx, double dy, double dz)
{
	if (PartitionIndex[0] != 0 && PartitionIndex[0] != PN - 1)
		return;
	//	- Inner Layer -
	int m0,m1,l0,l1;
	if (localvly[0] - vly[0]< 0.5*dy)
		m0 = 0;
	else
		m0 = (int)((localvly[0]-vly[0]-0.499999999*dy)/dy)+1;
	m1 = (int)((min(vly[1],localvly[1])-vly[0]-0.499999999*dy)/dy) + 1;
	if (localvlz[0] - vlz[0]< 0.5*dz)
		l0 = 0;
	else
		l0 = (int)((localvlz[0]-vlz[0]-0.499999999*dz)/dz)+1;
	l1 = (int)((min(vlz[1],localvlz[1])-vlz[0]-0.499999999*dz)/dz) + 1;
	if (dim == 3)
	{
		for (int j = m0;j<m1;j++)
			for (int k = l0;k<l1;k++)
			{
				if (PartitionIndex[0] == 0)
				{
					nn++;
					pos_veloc(vlx[0],vly[0]+(j+0.5)*dy,vlz[0]+(k+0.5)*dz,0,0,0);
					pressure_b(dx,dy,dz);
				}
				if (PartitionIndex[0] == PN - 1)
				{
					nn++;
					pos_veloc(vlx[1],vly[0]+(j+0.5)*dy,vlz[0]+(k+0.5)*dz,0,0,0);
					pressure_b(dx,dy,dz);
				}
			}
	}
	else
	{
		for (int k = 0;k<L;k++)
		{
			nn++;
			pos_veloc(   0,0,(k+0.5)*dz,0,0,0);
			pressure_b(dx,dy,dz);
			nn++;
			pos_veloc(N*dx,0,(k+0.5)*dz,0,0,0);
			pressure_b(dx,dy,dz);
		}

	}
	//  - Outer Layer -
	if (localvly[0] - vly[0]<= 0)
		m0 = 0;
	else
		m0 = (int)((localvly[0]-vly[0])/dy) + 1;
	m1 = (int)((min(vly[1],localvly[1])+(1e-10)*dy-vly[0])/dy);
	if (localvlz[0] - vlz[0]<= 0)
		l0 = 0;
	else
		l0 = (int)((localvlz[0]-vlz[0])/dz) + 1;
	l1 = (int)((min(vlz[1],localvlz[1])+(1e-10)*dz-vlz[0])/dz);
	if (dim == 3)
	{
		for (int j = m0;j<=m1;j++)
			for (int k = l0;k<=l1;k++)
			{
				if (PartitionIndex[0] == 0)
				{
					nn++;
					pos_veloc(vlx[0] - 0.5*dx,vly[0]+j*dy,vlz[0]+k*dz,0,0,0);
					pressure_b(dx,dy,dz);
				}
				if (PartitionIndex[0] == PN - 1)
				{
					nn++;
					pos_veloc(vlx[1] + 0.5*dx,vly[0]+j*dy,vlz[0]+k*dz,0,0,0);
					pressure_b(dx,dy,dz);
				}
			}
	}
	else
	{
		for (int k=0;k<=L;k++)
		{
			nn++;
			pos_veloc(      -0.5*dx,0,k*dz,0,0,0);
			pressure_b(dx,dy,dz);
			nn++;
			pos_veloc((N+0.5)*dx,0,k*dz,0,0,0);
			pressure_b(dx,dy,dz);
		}

	}

}

void STORAGE::Boundaries_Bottom(int N, int M, int L, double dx, double dy, double dz)
{
	if (PartitionIndex[2] != 0)
		return;
	//	- Inner Layer -
	int n0,n1,m0,m1;
	if (localvlx[0] - vlx[0]< 0.5*dx)
		n0 = 0;
	else
		n0 = (int)((localvlx[0]-vlx[0]-0.4999999*dx)/dx)+1;
	n1 = (int)((min(vlx[1],localvlx[1])-vlx[0]-0.4999999*dx)/dx)+1;
	if (localvly[0] - vly[0]< 0.5*dy)
		m0 = 0;
	else
		m0 = (int)((localvly[0]-vly[0]-0.4999999*dy)/dy)+1;
	m1 = (int)((min(vly[1],localvly[1])-vly[0]-0.4999999*dy)/dy)+1;
	if (dim == 3)
	{
		for (int i = n0;i<n1;i++)
			for (int j = m0;j<m1;j++)
			{
				if (PartitionIndex[2] == 0)
				{
					nn++;
					pos_veloc(vlx[0]+(i+0.5)*dx,vly[0]+(j+0.5)*dy,vlz[0],0,0,0);
					pressure_b(dx,dy,dz);
				}
			}
	}
	else
	{
		for (int i =0;i<N;i++)
		{
			nn++;
			pos_veloc((i+0.5)*dx,0,0,0,0,0);
			pressure_b(dx,dy,dz);
		}
	}
	//  - Outer Layer -

	if (localvlx[0] - vlx[0]<= 0)
		n0 = 0;
	else
		n0 = (int)((localvlx[0]-vlx[0])/dx) + 1;
	n1 = (int)((min(vlx[1],localvlx[1])+(1e-10)*dx-vlx[0])/dx);
	if (localvly[0] - vly[0]<= 0)
		m0 = 0;
	else
		m0 = (int)((localvly[0]-vly[0])/dy) + 1;
	m1 = (int)((min(vly[1],localvly[1])+(1e-10)*dy-vly[0])/dy);
	if (dim == 3)
	{
		for (int i = n0;i<=n1;i++)
			for (int j = m0;j<=m1;j++)
			{
				if (PartitionIndex[2] == 0)
				{
					nn++;
					pos_veloc(vlx[0]+i*dx,vly[0]+j*dy,vlz[0]-0.5*dz,0,0,0);
					pressure_b(dx,dy,dz);
				}
			}
	}
	else
	{
		for (int i =0;i<=N;i++)
		{
			nn++;
			pos_veloc(i*dx,0,-0.5*dz,0,0,0);
			pressure_b(dx,dy,dz);
		}
	}
}

void STORAGE::Fill_Part(int N,int M,int L,double dx,double dy,double dz, double u, double v, double w)
{
	int n0,n1,m0,m1,l0,l1;
	int i_intersectx = 0;
	int i_intersecty = 0;
	int i_intersectz = 0;
	// check intersection in X direction
	if ((XXmin >= localvlx[0] && XXmin <=localvlx[1] && XXmax >= localvlx[1]) ||
			(XXmin <= localvlx[0] && XXmax >=localvlx[0] && XXmax <= localvlx[1]) ||
			(XXmin >= localvlx[0] && XXmax <=localvlx[1]) ||
			(XXmin <= localvlx[0] && XXmax >=localvlx[1]))
		i_intersectx = 1;
	// check intersection in Y direction
	if ((YYmin >= localvly[0] && YYmin <=localvly[1] && YYmax >= localvly[1]) ||
			(YYmin <= localvly[0] && YYmax >=localvly[0] && YYmax <= localvly[1]) ||
			(YYmin >= localvly[0] && YYmax <=localvly[1]) ||
			(YYmin <= localvly[0] && YYmax >=localvly[1]))
		i_intersecty = 1;
	// check intersection in Z direction
	if ((ZZmin >= localvlz[0] && ZZmin <=localvlz[1] && ZZmax >= localvlz[1]) ||
			(ZZmin <= localvlz[0] && ZZmax >=localvlz[0] && ZZmax <= localvlz[1]) ||
			(ZZmin >= localvlz[0] && ZZmax <=localvlz[1]) ||
			(ZZmin <= localvlz[0] && ZZmax >=localvlz[1]))
		i_intersectz = 1;
	// Intersect in computational domain
	if (i_intersectx * i_intersecty * i_intersectz == 1)
	{
		localXXmin = std::max(XXmin,localvlx[0]);
		localXXmax = std::min(XXmax,localvlx[1]);
		localYYmin = std::max(YYmin,localvly[0]);
		localYYmax = std::min(YYmax,localvly[1]);
		localZZmin = std::max(ZZmin,localvlz[0]);
		localZZmax = std::min(ZZmax,localvlz[1]);
	}
	else
		return;
	//	n0 = (int)(localXXmin/dx);
	//	n1 = (int)(localXXmax/dx);
	//	m0 = (int)(localYYmin/dy);
	//	m1 = (int)(localYYmax/dy);
	//	l0 = (int)(localZZmin/dz);
	//	l1 = (int)(localZZmax/dz);
	if (lattice == 1 || lattice == 2)
	{
		if (localXXmin - vlx[0]< 0.5*dx)
			n0 = 0;
		else
			n0 = (int)((localXXmin-vlx[0]-0.4999999*dx)/dx)+1;
		n1 = (int)((localXXmax-vlx[0]-0.4999999*dx)/dx)+1;
		if (localYYmin - vly[0]< 0.5*dy)
			m0 = 0;
		else
			m0 = (int)((localYYmin-vly[0]-0.4999999*dy)/dy)+1;
		m1 = (int)((localYYmax-vly[0]-0.4999999*dy)/dy)+1;
		if (localZZmin - vlz[0]< 0.5*dz)
			l0 = 0;
		else
			l0 = (int)((localZZmin-vlz[0]-0.4999999*dz)/dz)+1;
		l1 = (int)((localZZmax-vlz[0]-0.4999999*dz)/dz)+1;
		if (dim == 3)
		{
			for (int i = n0;i<n1;i++)
				for (int j = m0;j<m1;j++)
					for (int k = l0;k<l1;k++)
					{

						if (lattice == 2)
						{
							nn++;
							pos_veloc(vlx[0]+(i+0.25)*dx,vly[0]+(j+0.25)*dx,vlz[0]+(k+0.25)*dx,0,0,0);
							pressure(dx,dy,dz);
							nn++;
							pos_veloc(vlx[0]+(i+0.75)*dx,vly[0]+(j+0.75)*dx,vlz[0]+(k+0.75)*dx,0,0,0);
							pressure(dx,dy,dz);
						}
						else
						{
							//for debugging
							double value;
							value = p0 * (vlx[0]+(i+0.5));
							//end for debugging
							nn++;
							pos_veloc(vlx[0]+(i+0.5)*dx,vly[0]+(j+0.5)*dx,vlz[0]+(k+0.5)*dx,0,0,0);
							pressure(dx,dy,dz,value);
						}
					}
		}
		else
		{
			for (int i = n0;i<n1;i++)
				for (int k = l0;k<l1;k++)
				{
					if (lattice == 2)
					{
						nn++;
						pos_veloc((i+0.25)*dx,0,(k+0.25)*dz,0,0,0);
						pressure(dx,dy,dz);
						nn++;
						pos_veloc((i+0.75)*dx,0,(k+0.75)*dz,0,0,0);
						pressure(dx,dy,dz);
					}
					else assert(0);
				}
		}
	}
	else
	{
		double h_r = dx;
		if (dim == 2)
		{
			assert(0);
			//			int h_m, h_n_odd, h_n_even;
			//			double h_l, h_w;
			//			h_l = localXXmax - localXXmin; //length of domain on x-axis
			//			h_w = localZZmax - localZZmin; //width of domain on z-axis
			//
			//			// compute # of rows (h_m)
			//			int h_g = (int)((h_w-h_r)/(sqrt(3)*h_r));
			//			if (h_r+(h_g+1)*sqrt(3)*h_r<=h_w)
			//				h_m = h_g+2;
			//			else
			//				h_m = h_g+1;
			//
			//			//compute # of columns (h_n_odd) for odd-numbered rows
			//			int h_f1 = (int) ((h_l-h_r)/h_r);
			//			if (h_f1%2 != 0)
			//				h_n_odd = (h_f1+1)/2;
			//			else
			//				h_n_odd = h_f1/2;
			//
			//			//compute # of columns (h_n_even) for even-numbered rows
			//			int h_f2 = (int) (h_l/h_r);
			//			if (h_f2%2 != 0)
			//				h_n_even = (h_f2+1)/2;
			//			else
			//				h_n_even = h_f2/2;
			//
			//			//compute total number of particles inside the cube
			//			if (h_m%2 != 0)
			//				nn = (h_m+1)/2*h_n_odd+(h_m-1)/2*h_n_even;
			//			else
			//				nn = h_m/2*(h_n_odd+h_n_even);
			//			cout<<"# of particles inside whole 2D domain is\t"<<nn<<endl;
			//
			//			//compute the location of particles inside the cube
			//			for (int i=1; i<=h_m; i++){ //rows
			//				if (i%2 != 0){ //odd-numbered rows
			//					for (int j=1; j<=h_n_odd; j++){ //columns
			//						pos_veloc(localXXmin+2*h_r+(j-1)*2*h_r,0,localZZmin+h_r+(i-1)*sqrt(3)*h_r,0,0,0);
			//						pressure(h_r,h_r,h_r);
			//					}
			//				}
			//				else{ //even-numbered rows
			//					for (int j=1; j<=h_n_even; j++){ //columns
			//						pos_veloc(localXXmin+h_r+(j-1)*2*h_r,0,localZZmin+h_r+(i-1)*sqrt(3)*h_r,0,0,0);
			//						pressure(h_r,h_r,h_r);
			//					}
			//				}
			//			}
		}
		if(dim == 3){
			int hh_m, hh_n_odd, hh_n_even, hh_u, hh_v_odd, hh_v_even, hh_d;
			double hh_l, hh_w, hh_h;
			hh_l = localXXmax-localXXmin; //length of domain on x-axis
			hh_w = localYYmax-localYYmin; //width of domain on y-axis
			hh_h = localZZmax-localZZmin; //height of domain on z-axis

			//compute # of layers (hh_d)
			int hh_e = (int)( (hh_h-h_r)/2/sqrt(6)/h_r*3 );
			if (h_r+(hh_e+1)*2*sqrt(6)*h_r/3 <= hh_h)
				hh_d = hh_e+2;
			else
				hh_d = hh_e+1;

			//compute # of rows (hh_m) for odd-numbered layers
			int hh_g = (int)((hh_w-h_r)/(sqrt(3)*h_r));
			if (h_r+(hh_g+1)*sqrt(3)*h_r<=hh_w)
				hh_m = hh_g+2;
			else
				hh_m = hh_g+1;

			//compute # of columns (hh_n_odd) for odd-numbered rows & odd-numbered layers
			int hh_f1 = (int) ((hh_l-h_r)/h_r);
			if (hh_f1%2 != 0)
				hh_n_odd = (hh_f1+1)/2;
			else
				hh_n_odd = hh_f1/2;

			//compute # of columns (hh_n_even) for even-numbered rows & odd-numbered layers
			int hh_f2 = (int) (hh_l/h_r);
			if (hh_f2%2 != 0)
				hh_n_even = (hh_f2+1)/2;
			else
				hh_n_even = hh_f2/2;

			//compute # of rows (hh_u) for even-numbered layers
			int hh_a = (int)( (hh_w+(1-sqrt(3)/3)*h_r)/sqrt(3)/h_r );
			if (sqrt(3)/3*h_r+(hh_a+1)*sqrt(3)*h_r<=hh_w)
				hh_u = hh_a+2;
			else
				hh_u = hh_a+1;

			//compute # of columns (hh_v_odd) for odd-numbered rows & even-numbered layers
			if (hh_f2%2 != 0)
				hh_v_odd = (hh_f2+1)/2;
			else
				hh_v_odd = hh_f2/2;

			//compute # of columns (hh_v_even) for even-numbered rows & even-numbered layers
			if (hh_f1%2 != 0)
				hh_v_even = (hh_f1+1)/2;
			else
				hh_v_even = hh_f1/2;

			//compute the total number of particles inside the cube
			int tempA, tempB;
			if (hh_d%2 != 0){ //there is an odd number of layers
				if (hh_m%2 != 0) //there is an odd number of rows
					tempA = ( (hh_m+1)/2*hh_n_odd+(hh_m-1)/2*hh_n_even ) * (hh_d+1)/2;
				else
					tempA = hh_m/2*(hh_n_odd+hh_n_even) * (hh_d+1)/2;

				if (hh_u%2 != 0) //there is an odd number of rows
					tempB = ( (hh_u+1)/2*hh_v_odd+(hh_u-1)/2*hh_v_even ) * (hh_d-1)/2;
				else
					tempB = hh_u/2*(hh_v_odd+hh_v_even) * (hh_d-1)/2;
			}
			else{ //there is an even number of layers
				if (hh_m%2 != 0) //there is an odd number of rows
					tempA = ( (hh_m+1)/2*hh_n_odd+(hh_m-1)/2*hh_n_even ) * hh_d/2;
				else
					tempA = hh_m/2*(hh_n_odd+hh_n_even) * hh_d/2;

				if (hh_u%2 != 0) //there is an odd number of rows
					tempB = ( (hh_u+1)/2*hh_v_odd+(hh_u-1)/2*hh_v_even ) * hh_d/2;
				else
					tempB = hh_u/2*(hh_v_odd+hh_v_even) * hh_d/2;
			}
			size_t nn_increase = tempA + tempB;
			//			cout<<"# of particles inside whole 3D domain is\t"<<nn_increase<<endl;

			//compute the location of particles inside the cube
			for (int i=1; i<=hh_d; i++){ //layers
				if (i%2 != 0){ //odd-numbered layers
					for (int j=1; j<=hh_m; j++){ //rows
						if (j%2 != 0){ //odd-numbered rows
							for (int k=1; k<=hh_n_odd; k++){ //columns
								++nn;
								pos_veloc(localXXmin+2*h_r+(k-1)*2*h_r,localYYmin+h_r+(j-1)*sqrt(3)*h_r,localZZmin+h_r+(i-1)*2*sqrt(6)*h_r/3,0,0,0);
								pressure(h_r,h_r,h_r);
							}
						}
						else{
							for (int k=1; k<=hh_n_even; k++){
								++nn;
								pos_veloc(localXXmin+h_r+(k-1)*2*h_r,localYYmin+h_r+(j-1)*sqrt(3)*h_r,localZZmin+h_r+(i-1)*2*sqrt(6)*h_r/3,0,0,0);
								pressure(h_r,h_r,h_r);
							}
						}
					}
				}
				else{ //even-numbered layers
					for (int j=1; j<=hh_u; j++){ //rows
						if (j%2 != 0){ //odd-numbered rows
							for (int k=1; k<=hh_v_odd; k++){ //columns
								++nn;;
								pos_veloc(localXXmin+h_r+(k-1)*2*h_r,localYYmin+sqrt(3)*h_r/3+(j-1)*sqrt(3)*h_r,localZZmin+h_r+(i-1)*2*sqrt(6)*h_r/3,0,0,0);
								pressure(h_r,h_r,h_r);
							}
						}
						else{
							for (int k=1; k<=hh_v_even; k++){
								++nn;;
								pos_veloc(localXXmin+2*h_r+(k-1)*2*h_r,localYYmin+sqrt(3)*h_r/3+(j-1)*sqrt(3)*h_r,localZZmin+h_r+(i-1)*2*sqrt(6)*h_r/3,0,0,0);
								pressure(h_r,h_r,h_r);
							}
						}
					}
				}
			}
		}
	}
}

void STORAGE::Drop(double dx,double dy,double dz,double xcen,double ycen,double zcen,double radius,double length,double u,double v,double w)
{
	double angle = atan(1);
	double L = length/2;
	u = 25 * sin(angle);
	v = 0;
	w = -25 * cos(angle);
	XXmin=xcen-(radius+L);
	XXmax=xcen+(radius+L);
	YYmin=ycen-(radius);
	YYmax=ycen+(radius);
	ZZmin=zcen-(radius+L);
	ZZmax=zcen+(radius+L);
	int n0,n1,m0,m1,l0,l1;
	int i_intersectx = 0;
	int i_intersecty = 0;
	int i_intersectz = 0;
	if ((XXmin >= localvlx[0] && XXmin <=localvlx[1] && XXmax >= localvlx[1]) ||
			(XXmin <= localvlx[0] && XXmax >=localvlx[0] && XXmax <= localvlx[1]) ||
			(XXmin >= localvlx[0] && XXmax <=localvlx[1]) ||
			(XXmin <= localvlx[0] && XXmax >=localvlx[1]))
		i_intersectx = 1;
	// check intersection in Y direction
	if ((YYmin >= localvly[0] && YYmin <=localvly[1] && YYmax >= localvly[1]) ||
			(YYmin <= localvly[0] && YYmax >=localvly[0] && YYmax <= localvly[1]) ||
			(YYmin >= localvly[0] && YYmax <=localvly[1]) ||
			(YYmin <= localvly[0] && YYmax >=localvly[1]))
		i_intersecty = 1;
	// check intersection in Z direction
	if ((ZZmin >= localvlz[0] && ZZmin <=localvlz[1] && ZZmax >= localvlz[1]) ||
			(ZZmin <= localvlz[0] && ZZmax >=localvlz[0] && ZZmax <= localvlz[1]) ||
			(ZZmin >= localvlz[0] && ZZmax <=localvlz[1]) ||
			(ZZmin <= localvlz[0] && ZZmax >=localvlz[1]))
		i_intersectz = 1;
	// Intersect in computational domain
	if (i_intersectx * i_intersecty * i_intersectz == 1)
	{
		localXXmin = std::max(XXmin,localvlx[0]);
		localXXmax = std::min(XXmax,localvlx[1]);
		localYYmin = std::max(YYmin,localvly[0]);
		localYYmax = std::min(YYmax,localvly[1]);
		localXXmin = std::max(XXmin,localvlx[0]);
		localXXmax = std::min(XXmax,localvlx[1]);
		localZZmin = std::max(ZZmin,localvlz[0]);
		localZZmax = std::min(ZZmax,localvlz[1]);
	}
	// no intersection
	// no intersection lionel.guo
	else
		return;
	if (lattice == 1 || lattice == 2)
	{
		if (localXXmin - vlx[0]< 0.5*dx)
			n0 = 0;
		else
			n0 = (int)((localXXmin-vlx[0]-0.4999999*dx)/dx)+1;
		n1 = (int)((localXXmax-vlx[0]-0.4999999*dx)/dx)+1;
		if (localYYmin - vly[0]< 0.5*dy)
			m0 = 0;
		else
			m0 = (int)((localYYmin-vly[0]-0.4999999*dy)/dy)+1;
		m1 = (int)((localYYmax-vly[0]-0.4999999*dy)/dy)+1;
		if (localZZmin - vlz[0]< 0.5*dz)
			l0 = 0;
		else
			l0 = (int)((localZZmin-vlz[0]-0.4999999*dz)/dz)+1;
		l1 = (int)((localZZmax-vlz[0]-0.4999999*dz)/dz)+1;
		double xct,zct,xcf,ycf,zcf,x,y,z,dist;
		if (dim == 3)
		{
			for (int i = n0;i<n1;i++)
				for (int j = m0;j<m1;j++)
					for (int k = l0;k<l1;k++)
					{
						x = vlx[0] + (i+0.5) * dx;
						y = vly[0] + (j+0.5) * dy;
						z = vlz[0] + (k+0.5) * dz;
						xct = x - xcen;
						ycf = y - ycen;
						zct = z - zcen;
						xcf = cos(-angle) * xct - sin(-angle) * zct;
						zcf = sin(-angle) * xct + cos(-angle) * zct;
						if(zcf > L)
						{
							dist = xcf * xcf + ycf * ycf + (zcf-L) * (zcf-L);
							dist = sqrt(dist) - radius;
						}
						else if(zcf > -L)
						{
							dist = xcf * xcf + ycf * ycf;
							dist = sqrt(dist) -radius;
						}
						else
						{
							dist = xcf * xcf + ycf * ycf + (zcf+L) * (zcf+L);
							dist = sqrt(dist) - radius;
						}
						if (dist < 0)
						{

							if (lattice == 2)
							{
								nn++;
								pos_veloc(vlx[0]+(i+0.25)*dx,vly[0]+(j+0.25)*dx,vlz[0]+(k+0.25)*dx,u,v,w);
								pressure(dx,dy,dz);
								nn++;
								pos_veloc(vlx[0]+(i+0.75)*dx,vly[0]+(j+0.75)*dx,vlz[0]+(k+0.75)*dx,u,v,w);
								pressure(dx,dy,dz);
							}
							else
							{
								//for debugging
								//end for debugging
								nn++;
								pos_veloc(vlx[0]+(i+0.5)*dx,vly[0]+(j+0.5)*dx,vlz[0]+(k+0.5)*dx,u,v,w);
								pressure(dx,dy,dz);
							}
						}
					}
		}
		else
		{
			for (int i = n0;i<n1;i++)
				for (int k = l0;k<l1;k++)
				{
					x = (i+0.5) * dx;
					z = (k+0.5) * dz;
					xct = x - xcen;
					zct = z - zcen;
					xcf = cos(-angle) * xct - sin(-angle) * zct;
					zcf = sin(-angle) * xct + cos(-angle) * zct;
					if(zcf > L)
					{
						dist = xcf * xcf + (zcf-L) * (zcf-L);
						dist = sqrt(dist) - radius;
					}
					else if(zcf > -L)
					{
						dist = xcf * xcf;
						dist = sqrt(dist) -radius;
					}
					else
					{
						dist = xcf * xcf + (zcf+L) * (zcf+L);
						dist = sqrt(dist) - radius;
					}
					if (dist < 0)
					{

						if (lattice == 2)
						{
							nn++;
							pos_veloc((i+0.25)*dx,0,(k+0.25)*dx,u,0,w);
							pressure(dx,dy,dz);
							nn++;
							pos_veloc((i+0.75)*dx,0,(k+0.75)*dx,u,0,w);
							pressure(dx,dy,dz);
						}
						else assert(0);
					}
				}
		}
	}
	else
	{
		double xct,zct,xcf,ycf,zcf,x,y,z,dist, h_r = dx;
		if(dim == 3){
			int hh_m, hh_n_odd, hh_n_even, hh_u, hh_v_odd, hh_v_even, hh_d;
			double hh_l, hh_w, hh_h;
			hh_l = localXXmax-localXXmin; //length of domain on x-axis
			hh_w = localYYmax-localYYmin; //width of domain on y-axis
			hh_h = localZZmax-localZZmin; //height of domain on z-axis

			//compute # of layers (hh_d)
			int hh_e = (int)( (hh_h-h_r)/2/sqrt(6)/h_r*3 );
			if (h_r+(hh_e+1)*2*sqrt(6)*h_r/3 <= hh_h)
				hh_d = hh_e+2;
			else
				hh_d = hh_e+1;

			//compute # of rows (hh_m) for odd-numbered layers
			int hh_g = (int)((hh_w-h_r)/(sqrt(3)*h_r));
			if (h_r+(hh_g+1)*sqrt(3)*h_r<=hh_w)
				hh_m = hh_g+2;
			else
				hh_m = hh_g+1;

			//compute # of columns (hh_n_odd) for odd-numbered rows & odd-numbered layers
			int hh_f1 = (int) ((hh_l-h_r)/h_r);
			if (hh_f1%2 != 0)
				hh_n_odd = (hh_f1+1)/2;
			else
				hh_n_odd = hh_f1/2;

			//compute # of columns (hh_n_even) for even-numbered rows & odd-numbered layers
			int hh_f2 = (int) (hh_l/h_r);
			if (hh_f2%2 != 0)
				hh_n_even = (hh_f2+1)/2;
			else
				hh_n_even = hh_f2/2;

			//compute # of rows (hh_u) for even-numbered layers
			int hh_a = (int)( (hh_w+(1-sqrt(3)/3)*h_r)/sqrt(3)/h_r );
			if (sqrt(3)/3*h_r+(hh_a+1)*sqrt(3)*h_r<=hh_w)
				hh_u = hh_a+2;
			else
				hh_u = hh_a+1;

			//compute # of columns (hh_v_odd) for odd-numbered rows & even-numbered layers
			if (hh_f2%2 != 0)
				hh_v_odd = (hh_f2+1)/2;
			else
				hh_v_odd = hh_f2/2;

			//compute # of columns (hh_v_even) for even-numbered rows & even-numbered layers
			if (hh_f1%2 != 0)
				hh_v_even = (hh_f1+1)/2;
			else
				hh_v_even = hh_f1/2;

			//compute the location of particles inside the cyclinder
			for (int i=1; i<=hh_d; i++){ //layers
				if (i%2 != 0){ //odd-numbered layers
					for (int j=1; j<=hh_m; j++){ //rows
						if (j%2 != 0){ //odd-numbered rows
							for (int k=1; k<=hh_n_odd; k++){ //columns

								x = localXXmin+2*h_r+(k-1)*2*h_r;
								y = localYYmin+h_r+(j-1)*sqrt(3)*h_r;
								z = localZZmin+h_r+(i-1)*2*sqrt(6)*h_r/3;
								xct = x - xcen;
								ycf = y - ycen;
								zct = z - zcen;
								xcf = cos(-angle) * xct - sin(-angle) * zct;
								zcf = sin(-angle) * xct + cos(-angle) * zct;
								if(zcf > L)
								{
									dist = xcf * xcf + ycf * ycf + (zcf-L) * (zcf-L);
									dist = sqrt(dist) - radius;
								}
								else if(zcf > -L)
								{
									dist = xcf * xcf + ycf * ycf;
									dist = sqrt(dist) -radius;
								}
								else
								{
									dist = xcf * xcf + ycf * ycf + (zcf+L) * (zcf+L);
									dist = sqrt(dist) - radius;
								}
								if (dist < 0)
								{
									nn++;
									pos_veloc(localXXmin+2*h_r+(k-1)*2*h_r,localYYmin+h_r+(j-1)*sqrt(3)*h_r,localZZmin+h_r+(i-1)*2*sqrt(6)*h_r/3,u,v,w);
									pressure(h_r,h_r,h_r);
								}
							}
						}
						else{
							for (int k=1; k<=hh_n_even; k++){
								x = localXXmin+h_r+(k-1)*2*h_r;
								y = localYYmin+h_r+(j-1)*sqrt(3)*h_r;
								z = localZZmin+h_r+(i-1)*2*sqrt(6)*h_r/3;
								xct = x - xcen;
								ycf = y - ycen;
								zct = z - zcen;
								xcf = cos(-angle) * xct - sin(-angle) * zct;
								zcf = sin(-angle) * xct + cos(-angle) * zct;
								if(zcf > L)
								{
									dist = xcf * xcf + ycf * ycf + (zcf-L) * (zcf-L);
									dist = sqrt(dist) - radius;
								}
								else if(zcf > -L)
								{
									dist = xcf * xcf + ycf * ycf;
									dist = sqrt(dist) -radius;
								}
								else
								{
									dist = xcf * xcf + ycf * ycf + (zcf+L) * (zcf+L);
									dist = sqrt(dist) - radius;
								}
								if (dist < 0)
								{
									nn++;
									pos_veloc(localXXmin+h_r+(k-1)*2*h_r,localYYmin+h_r+(j-1)*sqrt(3)*h_r,localZZmin+h_r+(i-1)*2*sqrt(6)*h_r/3,u,v,w);
									pressure(h_r,h_r,h_r);
								}
							}
						}
					}
				}
				else{ //even-numbered layers
					for (int j=1; j<=hh_u; j++){ //rows
						if (j%2 != 0){ //odd-numbered rows
							for (int k=1; k<=hh_v_odd; k++){ //columns
								x = localXXmin+h_r+(k-1)*2*h_r;
								y = localYYmin+sqrt(3)*h_r/3+(j-1)*sqrt(3)*h_r;
								z = localZZmin+h_r+(i-1)*2*sqrt(6)*h_r/3;
								xct = x - xcen;
								ycf = y - ycen;
								zct = z - zcen;
								xcf = cos(-angle) * xct - sin(-angle) * zct;
								zcf = sin(-angle) * xct + cos(-angle) * zct;
								if(zcf > L)
								{
									dist = xcf * xcf + ycf * ycf + (zcf-L) * (zcf-L);
									dist = sqrt(dist) - radius;
								}
								else if(zcf > -L)
								{
									dist = xcf * xcf + ycf * ycf;
									dist = sqrt(dist) -radius;
								}
								else
								{
									dist = xcf * xcf + ycf * ycf + (zcf+L) * (zcf+L);
									dist = sqrt(dist) - radius;
								}
								if (dist < 0)
								{
									nn++;
									pos_veloc(localXXmin+h_r+(k-1)*2*h_r,localYYmin+sqrt(3)*h_r/3+(j-1)*sqrt(3)*h_r,localZZmin+h_r+(i-1)*2*sqrt(6)*h_r/3,u,v,w);
									pressure(h_r,h_r,h_r);
								}
							}
						}
						else{
							for (int k=1; k<=hh_v_even; k++){
								x = localXXmin+2*h_r+(k-1)*2*h_r;
								y = localYYmin+sqrt(3)*h_r/3+(j-1)*sqrt(3)*h_r;
								z = localZZmin+h_r+(i-1)*2*sqrt(6)*h_r/3;
								xct = x - xcen;
								ycf = y - ycen;
								zct = z - zcen;
								xcf = cos(-angle) * xct - sin(-angle) * zct;
								zcf = sin(-angle) * xct + cos(-angle) * zct;
								if(zcf > L)
								{
									dist = xcf * xcf + ycf * ycf + (zcf-L) * (zcf-L);
									dist = sqrt(dist) - radius;
								}
								else if(zcf > -L)
								{
									dist = xcf * xcf + ycf * ycf;
									dist = sqrt(dist) -radius;
								}
								else
								{
									dist = xcf * xcf + ycf * ycf + (zcf+L) * (zcf+L);
									dist = sqrt(dist) - radius;
								}
								if (dist < 0)
								{
									nn++;
									pos_veloc(localXXmin+2*h_r+(k-1)*2*h_r,localYYmin+sqrt(3)*h_r/3+(j-1)*sqrt(3)*h_r,localZZmin+h_r+(i-1)*2*sqrt(6)*h_r/3,u,v,w);
									pressure(h_r,h_r,h_r);
								}
							}
						}
					}
				}
			}
		}
	}
}

void STORAGE::pos_veloc(double x,double y,double z,double u,double v,double w)
{
	fld.xp.push_back(x);
	fld.yp.push_back(y);
	fld.zp.push_back(z);
	fld.up.push_back(u);
	fld.vp.push_back(v);
	fld.wp.push_back(w);
	fld.iflag.push_back(1);

}

void STORAGE::pressure(double dx,double dy,double dz)
{
	fld.rho.push_back(rho0);
	fld.p.push_back(p0);
	fld.pm.push_back(vnorm_mass*fld.rho.back()*dx*dy*dz);
}

void STORAGE::pressure(double dx,double dy,double dz,double value)
{
	fld.rho.push_back(rho0);
	fld.p.push_back(value);
	fld.pm.push_back(vnorm_mass*fld.rho.back()*dx*dy*dz);
}

void STORAGE::pressure_b(double dx,double dy,double dz)
{
	fld.rho.push_back(rho0);
	fld.p.push_back(0);
	fld.pm.push_back(0.5*fld.rho.back()*dx*dy*dz);
}

void STORAGE::Initialization(string input)
{
	ifstream is(input.c_str());
	string textline;
	vector<string> lines_of_text;
	while (getline(is, textline))
		lines_of_text.push_back(textline);
	vector<string>::size_type line_num=0;
	istringstream line(lines_of_text[line_num]);
	line>>i_restartRun;
	line_num++;
	line.str(lines_of_text[line_num]);
	line>>i_kernel;
	line_num++;
	line.str(lines_of_text[line_num]);
	line>>i_algorithm;
	line_num++;
	line.str(lines_of_text[line_num]);
	line>>i_densityFilter;
	line_num++;
	if (i_densityFilter > 0)
	{
		line.str(lines_of_text[line_num]);
		line>>ndt_FilterPerform;
		line_num++;
	}
	else
		ndt_FilterPerform = 0;
	line.str(lines_of_text[line_num]);
	line>>i_kernelcorrection;
	line_num++;
	i_periodicOBs[0] = 0;
	i_periodicOBs[1] = 0;
	i_periodicOBs[2] = 0;
	line.str(lines_of_text[line_num]);
	line>>i_viscos;
	line_num++;
	if (i_viscos == 1)
	{
		line.str(lines_of_text[line_num]);
		line>>viscos_val;
		line_num++;
	}
	else
	{
		line.str(lines_of_text[line_num]);
		line>>viscos_val;
		line_num++;
	}
	line.str(lines_of_text[line_num]);
	line>>i_vort;
	line_num++;
	line.clear();
	line.str(lines_of_text[line_num]);
	line>>g;
	line_num++;
	line.clear();
	line.str(lines_of_text[line_num]);
	line>>dim;
	line_num++;
	line.str(lines_of_text[line_num]);
	line>>i_EoS;
	line_num++;
	line.str(lines_of_text[line_num]);
	line>>h_SWL;
	line_num++;
	if (i_EoS == 1 || i_EoS == 3)
	{
		line.str(lines_of_text[line_num]);
		line>>coef;
		line_num++;
		rho0 = 13;
		gamma = 7;
		TE0 = 0;
		p0 = 0;
		B = coef * coef * 9.81 * h_SWL * rho0/ gamma;
	}
	else if (i_EoS == 2)
	{
		coef = 0;
		rho0 = 1;
		gamma = 1.4;
		line.str(lines_of_text[line_num]);
		line>>TE0;
		line_num++;
	}
	else
	{
		coef = 0;
		line.clear();
		line.str(lines_of_text[line_num]);
		line>>rho0;
		line_num++;
		line.clear();
		line.str(lines_of_text[line_num]);
		line>>p0;
		line_num++;
		line.clear();
		line.str(lines_of_text[line_num]);
		line>>pinf;
		line_num++;
		line.clear();
		line.str(lines_of_text[line_num]);
		line>>einf;
		line_num++;
		line.clear();
		line.str(lines_of_text[line_num]);
		line>>gamma;
		line_num++;

	}
	line.clear();
	line.str(lines_of_text[line_num]);
	line>>IBC;
	line_num++;
	if (IBC == 1)
	{
		line.str(lines_of_text[line_num]);
		line>>visc_wall;
		line_num++;
		ndt_DBCPerform = 0;
	}
	else
	{
		visc_wall = 0;
		line.str(lines_of_text[line_num]);
		line>>ndt_DBCPerform;
		line_num++;
		if (ndt_DBCPerform == 0) ndt_DBCPerform = 1;
		if (i_algorithm == 2 && ndt_DBCPerform  != 1) assert(false);
		else if (i_algorithm == 2) {}
	}
	line.str(lines_of_text[line_num]);
	line>>lattice;
	line_num++;
	if (lattice == 1)
		vnorm_mass = 1;
	else if (lattice == 2)
		vnorm_mass = 0.5;
	else
		vnorm_mass = sqrt(2)/12;
	line.str(lines_of_text[line_num]);
	if (dim == 3)
		line>>vlx[0]>>vly[0]>>vlz[0];
	else if (dim == 2)
	{
		line>>vlx[0]>>vlz[0];
		vly[0] = 0;
	}
	else
		line>>vlx[0];
	line_num++;
	line.clear();
	line.str(lines_of_text[line_num]);
	if (dim == 3)
		line>>vlx[1]>>vly[1]>>vlz[1];
	else if (dim == 2)
	{
		line>>vlx[1]>>vlz[1];
		vly[1] = 0;
	}
	else
	{
		line>>vlx[1];
		vlz[1] = vly[1] = 0;

	}
	line_num++;
	line.str(lines_of_text[line_num]);
	if (dim == 3)
		line>>dx>>dy>>dz;
	else
	{
		line>>dx>>dz;
		dy = 0;
	}
	line_num++;
	line.str(lines_of_text[line_num]);
	line>>coefficient;
	h=coefficient*sqrt(dx*dx+dy*dy+dz*dz);
	line_num++;
	N = (int)((vlx[1] - vlx[0] + 0.000001*dx)/dx);
	M = (int)((vly[1] - vly[0] + 0.000001*dy)/dy);
	L = (int)((vlz[1] - vlz[0] + 0.000001*dz)/dz);
	zmax = L*dz;
	line.str(lines_of_text[line_num]);
	line>>i_periodicOBs[0]>>i_periodicOBs[1]>>i_periodicOBs[2];
	line_num++;
	ParallelInitialization();
	nn = 0;  //Initial number of particles
	if (dim == 3 && i_periodicOBs[1] !=1 && IBC != 3)
		Boundaries_FrontBack(N,M,L,dx,dy,dz);
	if (i_periodicOBs[0] != 1 && IBC != 3)
		Boundaries_LeftRight(N,M,L,dx,dy,dz);
	if (i_periodicOBs[2] != 1 && IBC != 3)
		Boundaries_Bottom(N,M,L,dx,dy,dz);
	nbb1 = nn;
	nbb = nbb1;
	nb = nbb;
	int fill_flag_part = 0;
	line.clear();
	line.str(lines_of_text[line_num]);
	line>>fill_flag_part;
	line_num++;
	if (fill_flag_part == 1 && dim == 3)
	{
		line.clear();
		line.str(lines_of_text[line_num]);
		line>>i_correct_pb;
		line_num++;
		line.clear();
		line.str(lines_of_text[line_num]);
		line>>XXmin>>XXmax;
		line_num++;
		line.clear();
		line.str(lines_of_text[line_num]);
		line>>YYmin>>YYmax;
		line_num++;
		line.clear();
		line.str(lines_of_text[line_num]);
		line>>ZZmin>>ZZmax;
		line_num++;
		line.clear();
		line.str(lines_of_text[line_num]);
		double u,v,w;
		line>>u>>v>>w;
		line_num++;
		Fill_Part(N,M,L,dx,dy,dz,u,v,w);
	}
	else if (fill_flag_part == 1 && dim == 2)
	{
		line.clear();
		line.str(lines_of_text[line_num]);
		line>>i_correct_pb;
		line_num++;
		line.clear();
		line.str(lines_of_text[line_num]);
		line>>XXmin>>XXmax;
		line_num++;
		line.clear();
		line.str(lines_of_text[line_num]);
		line>>ZZmin>>ZZmax;
		line_num++;
		double u,w;
		line>>u>>w;
		Fill_Part(N,0,L,dx,0,dz,u,0,w);
	}
	int fill_flag_drop =0;
	line.clear();
	line.str(lines_of_text[line_num]);
	line>>fill_flag_drop;
	line_num++;
	if (fill_flag_drop == 1 && dim == 3)
	{
		double xcen,ycen,zcen,radius,length,u,v,w;
		line.clear();
		line.str(lines_of_text[line_num]);
		line>>xcen>>ycen>>zcen;
		line_num++;
		line.clear();
		line.str(lines_of_text[line_num]);
		line>>radius;
		line_num++;
		line.clear();
		line.str(lines_of_text[line_num]);
		line>>length;
		line_num++;
		line.clear();
		line.str(lines_of_text[line_num]);
		line>>u>>v>>w;
		line_num++;
		Drop(dx,dy,dz,xcen,ycen,zcen,radius,length,u,v,w);
	}
	if (fill_flag_drop == 1 && dim == 2)
	{
		double xcen,zcen,radius,length,u,w;
		line.clear();
		line.str(lines_of_text[line_num]);
		line>>xcen>>zcen;
		line_num++;
		line.clear();
		line.str(lines_of_text[line_num]);
		line>>radius;
		line_num++;
		line.clear();
		line.str(lines_of_text[line_num]);
		line>>length;
		line_num++;
		line.clear();
		line.str(lines_of_text[line_num]);
		line>>u>>w;
		line_num++;
		Drop(dx,0,dz,xcen,0,zcen,radius,length,u,0,w);
	}
	np = nn;
	line.clear();
	line.str(lines_of_text[line_num]);
	line>>tmax>>out;
	line_num++;
	line.str(lines_of_text[line_num]);
	line>>trec_ini;
	line_num++;
	line.str(lines_of_text[line_num]);
	line>>dtrec_det>>t_sta_det>>t_end_det;
	line_num++;
	line.str(lines_of_text[line_num]);
	line>>dt>>ivar_dt;
	line_num++;
	line.str(lines_of_text[line_num]);
	line>>CFL_number;
	line_num++;
	line.str(lines_of_text[line_num]);
	line>>iRiemannSolver;
	line_num++;
	//for Tait
	double cmax = sqrt(B*gamma/rho0);
	double dtmax = CFL_number * dx / (1.1*cmax);
	if (dt > dtmax) dt = dtmax;
	dtmax = CFL_number * dx *dx / (101*pow(10,-6));
	if (dt > dtmax) dt = dtmax;
	//for Tait end
	if (PN != 0 || PL != 0 || PL !=0)
		ParallelAdjustment();
	ParticleNumberSync();
	allocation();
	for (size_t i=0;i<np;i++)
	{
		xp[i] = fld.xp[i];
		yp[i] = fld.yp[i];
		zp[i] = fld.zp[i];
		up[i] = fld.up[i];
		vp[i] = fld.vp[i];
		wp[i] = fld.wp[i];
		rho[i] = fld.rho[i];
		p[i] = fld.p[i];
		pm[i] = fld.pm[i];
	}
	fld.xp.clear();
	fld.yp.clear();
	fld.zp.clear();
	fld.up.clear();
	fld.vp.clear();
	fld.wp.clear();
	fld.TE.clear();
	fld.rho.clear();
	fld.pm.clear();
	fld.p.clear();
	fld.cs.clear();
}

void STORAGE::divide(int start, int end, int kind)
{
	//	ncn = int(((N+1)*dx+0.2*h)/(2*h))+1;
	//	if (dim == 3)
	//		ncm = int(((M+1)*dx+0.2*h)/(2*h))+1;
	//	else
	//		ncm = 1;
	//	ncl = int((zmax+0.2*h)/(2*h))+1;
	//	nx = ncn*ncm*ncl;
	if (i_extendz == 1 && (PartitionIndex[2] + 1) == PL)
	{
		ncl++;
		local_ncl++;
		local_ncl_interior[1]++;
		nx = local_ncn * local_ncm * local_ncl;
		zmax = ncl * 2 * h + vlz[0] - 1.5*dz;
		i_extendz = 0;
	}
	if (kind == 1)
	{
		nc_k1.clear();
		vector< vector<size_t> >::iterator it= ibox_k1.begin();
		for (;it!=ibox_k1.end();it++)
			(*it).clear();
		ibox_k1.clear();
		nc_k1.assign(nx,0);
		ibox_k1.resize(nx);
		it= ibox_k1.begin();
		for (;it!=ibox_k1.end();it++)
			(*it).clear();
	}
	else
	{
		nc_k2.clear();
		vector< vector<size_t> >::iterator it= ibox_k2.begin();
		for (;it!=ibox_k2.end();it++)
			(*it).clear();
		ibox_k2.clear();
		nc_k2.assign(nx,0);
		ibox_k2.resize(nx);
		it= ibox_k2.begin();
		for (;it!=ibox_k2.end();it++)
			(*it).clear();
	}
	for (int i = start;i<end;i++)
	{
		if (dim == 3)
		{
			double distx,disty,distz;
			distx = xp[i] - xmin;
			disty = yp[i] - ymin;
			distz = zp[i] - zmin;
			int icell,jcell,kcell;
			icell = (int)( distx / (2*h) ) + 1;
			jcell = (int)( disty / (2*h) ) + 1;
			kcell = (int)( distz / (2*h) ) + 1;
			int ii;
			ii    = icell - 1 + (jcell - 1)*local_ncn +(kcell - 1)*local_ncn*local_ncm;
			(kind == 1)?((nc_k1[ii])++):((nc_k2[ii])++);
			(kind == 1)?(ibox_k1[ii].push_back(i)):(ibox_k2[ii].push_back(i));
		}
		else
		{
			double distx,distz;
			distx = xp[i] - xmin;
			distz = zp[i] - zmin;
			int icell,kcell;
			icell = (int)( distx / (2*h) ) + 1;
			kcell = (int)( distz / (2*h) ) + 1;
			int ii;
			ii    = icell - 1 + (kcell - 1)*local_ncn;
			(kind == 1)?((nc_k1[ii])++):((nc_k2[ii])++);
			(kind == 1)?(ibox_k1[ii].push_back(i)):(ibox_k2[ii].push_back(i));
		}
	}
}

void STORAGE::print_divide()
{
	//	vector< vector<size_t> >::iterator itboxk1 = ibox_k1.begin();
	//	vector< vector<size_t> >::iterator itboxk2 = ibox_k2.begin();
	//	FILE *outfile;
	//	outfile = fopen("divide.txt","w");
	//	fprintf(outfile,"DIVIDE\n");
	//	fprintf(outfile,"KIND_1\n");
	//	int num=0;
	//	for(;itboxk1!=ibox_k1.end();++itboxk1)
	//	{
	//		vector<size_t>::iterator it=(*itboxk1).begin();
	//		++num;
	//		fprintf(outfile,"Box %d\n",num);
	//		if((*itboxk1).size()>0)
	//		{
	//			for(;it!=(*itboxk1).end();it++)
	//			{
	//				fprintf(outfile,"%lu ",*it);
	//			}
	//			fprintf(outfile,"\n");
	//		}
	//		else
	//			fprintf(outfile,"No kind_1 particle in cell %d\n",num);
	//	}
	//	fprintf(outfile,"KIND_2\n");
	//	num = 0;
	//	for(;itboxk2!=ibox_k2.end();++itboxk2)
	//	{
	//		vector<size_t>::iterator it=(*itboxk2).begin();
	//		++num;
	//		fprintf(outfile,"Box %d\n",num);
	//		if((*itboxk2).size()>0)
	//		{
	//			for(;it!=(*itboxk2).end();it++)
	//			{
	//				fprintf(outfile,"%lu ",*it);
	//			}
	//			fprintf(outfile,"\n");
	//		}
	//		else
	//			fprintf(outfile,"No kind_2 particle in cell %d\n",num);
	//	}
	//	fclose(outfile);
}
const char* STORAGE::right_flush(
		int	   n,
		int	   ndigits)
{
	static	char	s[20];
	int		i;

	if (n == 0)
		ndigits--;
	for (i = n; i > 0; i /= 10) ndigits--;

	s[0] = '\0';
	for (;ndigits > 0; ndigits--)
		(void) sprintf(s,"%s0",s);
	(void) sprintf(s,"%s%d",s,n);
	return s;
}

void STORAGE::print_out(int step,double time,const char* outputname)
{
	size_t i;
	char filename[200];
	FILE *outfile;
	sprintf(filename,"vtkoutput_%s",right_flush(step,7));
	sprintf(filename,"%s_%s",filename,outputname);
	if (PN * PM * PL != 1)
		sprintf(filename,"%s-nd%s",filename,right_flush(id,4));
	sprintf(filename,"%s.vtk",filename);
	outfile = fopen(filename,"w");
	double d = 0;
	fprintf(outfile,"# vtk DataFile Version 3.0\n");
	fprintf(outfile,"The actual time is %.8f\n",time);
	fprintf(outfile,"ASCII\n");
	fprintf(outfile,"DATASET POLYDATA\n");
	fprintf(outfile,"POINTS %lu double\n",np-nb);
	for (i =nb;i<np;i++)
	{
		if (dim == 3)
			fprintf(outfile,"%.16g %.16g %.16g\n",xp[i],yp[i],zp[i]);
		else
			fprintf(outfile,"%.16g %.16g %.16g\n",xp[i],d,zp[i]);
	}
	fprintf(outfile,"POINT_DATA %lu\n",np-nb);
	fprintf(outfile,"SCALARS pressure double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for (i = nb;i<np;i++)
	{
		fprintf(outfile,"%.16g\n",p[i]);
	}
	fprintf(outfile,"SCALARS density double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for (i =nb;i<np;i++)
	{
		fprintf(outfile,"%.16g\n",rho[i]);
	}
	fprintf(outfile,"VECTORS velocity double\n");
	for (i =nb;i<np;i++)
	{
		if (dim == 3)
			fprintf(outfile,"%.16g %.16g %.16g\n",up[i],vp[i],wp[i]);
		else
			fprintf(outfile,"%.16g %.16g %.16g\n",xp[i],d,zp[i]);
	}
	if (detail)
	{
		fprintf(outfile,"SCALARS PM double\n");
		fprintf(outfile,"LOOKUP_TABLE default\n");
		for (i =nb;i<np;i++)
		{
			fprintf(outfile,"%.16g\n",pm[i]);
		}
		fprintf(outfile,"SCALARS TE double\n");
		fprintf(outfile,"LOOKUP_TABLE default\n");
		for (i =nb;i<np;i++)
		{
			fprintf(outfile,"%.16g\n",TE[i]);
		}
	}
	if (i_Print_J == 1)
	{
		fprintf(outfile,"VECTORS Current double\n");
		for (i =nb;i<np;i++)
		{
			fprintf(outfile,"%.16g %.16g %.16g\n",Jx[i],Jy[i],Jz[i]);
		}
		fprintf(outfile,"VECTORS GradientPhi double\n");
		for (i =nb;i<np;i++)
		{
			fprintf(outfile,"%.16g %.16g %.16g\n",GPx[i],GPy[i],GPz[i]);
		}
		fprintf(outfile,"SCALARS DivergenceUxB double\n");
		fprintf(outfile,"LOOKUP_TABLE default\n");
		for (i =nb;i<np;i++)
		{
			fprintf(outfile,"%.16g\n",dUxB[i]);
		}
	}
	fclose(outfile);
	return;
}

void STORAGE::print_out_cell(int step,double time)
{
	char filename[200];
	FILE *outfile;
	sprintf(filename,"vtkoutput_%s",right_flush(step,7));
	if (PN * PM * PL != 1)
		sprintf(filename,"%s-nd%s",filename,right_flush(id,4));
	sprintf(filename,"%s.vtk",filename);
	outfile = fopen(filename,"w");
	fprintf(outfile,"# vtk DataFile Version 3.0\n");
	fprintf(outfile,"The actual time is %.8f\n",time);
	fprintf(outfile,"ASCII\n");
	fprintf(outfile,"DATASET POLYDATA\n");
	int number=0;
	for (int j = local_ncm_interior[0];j<local_ncm_interior[1]; j++)
		for (int k = local_ncl_interior[0]; k<local_ncl_interior[1];k++)
		{		      
			int IndexDisplacedAndBuffer = Index3D(1,j,k,local_ncn,local_ncm,local_ncl);
			number += nc_k2[IndexDisplacedAndBuffer];
		}
	fprintf(outfile,"POINTS %d double\n",number);
	//	for (i =nb;i<np;i++)
	//	{
	//		if (dim == 3)
	//			fprintf(outfile,"%.16g %.16g %.16g\n",xp[i],yp[i],zp[i]);
	//		else
	//			fprintf(outfile,"%.16g %.16g %.16g\n",xp[i],d,zp[i]);
	//	}
	for (int j = local_ncm_interior[0];j<local_ncm_interior[1]; j++)
		for (int k = local_ncl_interior[0]; k<local_ncl_interior[1];k++)
		{
			int IndexDisplacedAndBuffer = Index3D(1,j,k,local_ncn,local_ncm,local_ncl);
			vector<size_t>::iterator it = ibox_k2[IndexDisplacedAndBuffer].begin();
			if (nc_k2[IndexDisplacedAndBuffer] != 0)
			{
				for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
				{
					size_t num;
					num = *it;
					fprintf(outfile,"%.16g %.16g %.16g\n",xp[num],yp[num],zp[num]);
				}
			}
		}
	fprintf(outfile,"POINT_DATA %d\n",number);
	fprintf(outfile,"SCALARS pressure double\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	for (int j = local_ncm_interior[0];j<local_ncm_interior[1]; j++)
		for (int k = local_ncl_interior[0]; k<local_ncl_interior[1];k++)
		{		  
			int IndexDisplacedAndBuffer = Index3D(1,j,k,local_ncn,local_ncm,local_ncl);
			vector<size_t>::iterator it = ibox_k2[IndexDisplacedAndBuffer].begin();
			if (nc_k2[IndexDisplacedAndBuffer] != 0)
			{
				for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
				{
					size_t num;
					num = *it;
					fprintf(outfile,"%.16g\n",p[num]);
				}
			}
		}
	fprintf(outfile,"SCALARS index integer\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");
	int number_bk=0;
	for (int j = local_ncm_interior[0];j<local_ncm_interior[1]; j++)
		for (int k = local_ncl_interior[0]; k<local_ncl_interior[1];k++)
		{
			int IndexDisplacedAndBuffer = Index3D(1,j,k,local_ncn,local_ncm,local_ncl);
			vector<size_t>::iterator it = ibox_k2[IndexDisplacedAndBuffer].begin();
			if (nc_k2[IndexDisplacedAndBuffer] != 0)
			{
				for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
				{
					size_t num;
					num = *it;
					fprintf(outfile,"%lu\n",num);
					++number_bk;
				}
			}
		}
	if (number_bk != number)
		assert(0);
	//	for (i =nb;i<np;i++)
	//	{
	//		fprintf(outfile,"%.16g\n",p[i]);
	//	}
	//	fprintf(outfile,"SCALARS density double\n");
	//	fprintf(outfile,"LOOKUP_TABLE default\n");
	//	for (i =nb;i<np;i++)
	//	{
	//		fprintf(outfile,"%.16g\n",rho[i]);
	//	}
	//	fprintf(outfile,"VECTORS velocity double\n");
	//	for (i =nb;i<np;i++)
	//	{
	//		if (dim == 3)
	//			fprintf(outfile,"%.16g %.16g %.16g\n",up[i],vp[i],wp[i]);
	//		else
	//			fprintf(outfile,"%.16g %.16g %.16g\n",xp[i],d,zp[i]);
	//	}
	//	if (detail)
	//	{
	//		fprintf(outfile,"SCALARS PM double\n");
	//		fprintf(outfile,"LOOKUP_TABLE default\n");
	//		for (i =nb;i<np;i++)
	//		{
	//			fprintf(outfile,"%.16g\n",pm[i]);
	//		}
	//		fprintf(outfile,"SCALARS TE double\n");
	//		fprintf(outfile,"LOOKUP_TABLE default\n");
	//		for (i =nb;i<np;i++)
	//		{
	//			fprintf(outfile,"%.16g\n",TE[i]);
	//		}
	//	}
	//	if (i_Print_J == 1)
	//	{
	//		fprintf(outfile,"VECTORS Current double\n");
	//		for (i =nb;i<np;i++)
	//		{
	//			fprintf(outfile,"%.16g %.16g %.16g\n",Jx[i],Jy[i],Jz[i]);
	//		}
	//		fprintf(outfile,"VECTORS GradientPhi double\n");
	//		for (i =nb;i<np;i++)
	//		{
	//			fprintf(outfile,"%.16g %.16g %.16g\n",GPx[i],GPy[i],GPz[i]);
	//		}
	//		fprintf(outfile,"SCALARS DivergenceUxB double\n");
	//		fprintf(outfile,"LOOKUP_TABLE default\n");
	//		for (i =nb;i<np;i++)
	//		{
	//			fprintf(outfile,"%.16g\n",dUxB[i]);
	//		}
	//	}
	fclose(outfile);
	return;
}

void STORAGE::gradients_calc(int i,int j,double dux,double duy,double duz)
{
	//	double Vj = fld.pm[j] / fld.rho[j];
	//
	//	double frx_norm_i = Vj*frxi;
	//	double fry_norm_i = Vj*fryi;
	//	double frz_norm_i = Vj*frzi;
	//
	//	dudx_CSPH[i] = dudx_CSPH[i] - dux*frx_norm_i;
	//	dudy_CSPH[i] = dudy_CSPH[i] - dux*fry_norm_i;
	//	dudz_CSPH[i] = dudz_CSPH[i] - dux*frz_norm_i;
	//	dvdx_CSPH[i] = dvdx_CSPH[i] - duy*frx_norm_i;
	//	dvdy_CSPH[i] = dvdy_CSPH[i] - duy*fry_norm_i;
	//	dvdz_CSPH[i] = dvdz_CSPH[i] - duy*frz_norm_i;
	//	dwdx_CSPH[i] = dwdx_CSPH[i] - duz*frx_norm_i;
	//	dwdy_CSPH[i] = dwdy_CSPH[i] - duz*fry_norm_i;
	//	dwdz_CSPH[i] = dwdz_CSPH[i] - duz*frz_norm_i;
	//
	//	double Vi = fld.pm[i]/fld.rho[i];
	//
	//	double frx_norm_j = Vi*frxj;
	//	double fry_norm_j = Vi*fryj;
	//	double frz_norm_j = Vi*frzj;
	//
	//	dudx_CSPH[j] = dudx_CSPH[j] - dux*frx_norm_j;
	//	dudy_CSPH[j] = dudy_CSPH[j] - dux*fry_norm_j;
	//	dudz_CSPH[j] = dudz_CSPH[j] - dux*frz_norm_j;
	//	dvdx_CSPH[j] = dvdx_CSPH[j] - duy*frx_norm_j;
	//	dvdy_CSPH[j] = dvdy_CSPH[j] - duy*fry_norm_j;
	//	dvdz_CSPH[j] = dvdz_CSPH[j] - duy*frz_norm_j;
	//	dwdx_CSPH[j] = dwdx_CSPH[j] - duz*frx_norm_j;
	//	dwdy_CSPH[j] = dwdy_CSPH[j] - duz*fry_norm_j;
	//	dwdz_CSPH[j] = dwdz_CSPH[j] - duz*frz_norm_j;
	//
	//	double drho = fld.rho[i] - fld.rho[j];
	//	drhodx_CSPH[i] = drhodx_CSPH[i] - drho*frx_norm_i;
	//	drhody_CSPH[i] = drhody_CSPH[i] - drho*fry_norm_i;
	//	drhodz_CSPH[i] = drhodz_CSPH[i] - drho*frz_norm_i;
	//	drhodx_CSPH[j] = drhodx_CSPH[j] - drho*frx_norm_j;
	//	drhody_CSPH[j] = drhody_CSPH[j] - drho*fry_norm_j;
	//	drhodz_CSPH[j] = drhodz_CSPH[j] - drho*frz_norm_j;
	//
	//	double dTE = fld.TE[i]- fld.TE[j];
	//	dTEdx_CSPH[i] = dTEdx_CSPH[i] - dTE*frx_norm_i;
	//	dTEdy_CSPH[i] = dTEdy_CSPH[i] - dTE*fry_norm_i;
	//	dTEdz_CSPH[i] = dTEdz_CSPH[i] - dTE*frz_norm_i;
	//	dTEdx_CSPH[j] = dTEdx_CSPH[j] - dTE*frx_norm_j;
	//	dTEdy_CSPH[j] = dTEdy_CSPH[j] - dTE*fry_norm_j;
	//	dTEdz_CSPH[j] = dTEdz_CSPH[j] - dTE*frz_norm_j;
}

void STORAGE::viscosity(double dot,double drx,double dry,double drz,double dux,double duy,double duz,double rr2,double cbar,double robar,double one_over_rhobar,int i,int j,int j1,int j2)
{
	double pi_visc;
	if (dot<0)        // ! aproaching
	{
		double amubar = h*dot/(rr2 + 0.01*h*h);
		pi_visc = - viscos_val * cbar * amubar / robar;
	}
	else             //  ! going away
		pi_visc = 0;

	if (dim == 3)
	{
		ax[i] = ax[i] - pm[j] * pi_visc * frxi;
		ay[i] = ay[i] - pm[j] * pi_visc * fryi;
		az[i] = az[i] - pm[j] * pi_visc * frzi;

		ax[j] = ax[j] + pm[i] * pi_visc * frxj;
		ay[j] = ay[j] + pm[i] * pi_visc * fryj;
		az[j] = az[j] + pm[i] * pi_visc * frzj;

		term2i=0.5 * pi_visc *( frxi*dux+fryi*duy+frzi*duz);
		term2j=0.5 * pi_visc *( frxj*dux+fryj*duy+frzj*duz);
	}
	else
	{
		ax[i] = ax[i] - pm[j] * pi_visc * frxi;
		az[i] = az[i] - pm[j] * pi_visc * frzi;

		ax[j] = ax[j] + pm[i] * pi_visc * frxj;
		az[j] = az[j] + pm[i] * pi_visc * frzj;

		term2i=0.5 * pi_visc *( frxi*dux+frzi*duz);
		term2j=0.5 * pi_visc *( frxj*dux+frzj*duz);
	}
}

void STORAGE::allocation()
{
	//	np_max = np_max + np_max / 3;
	xp = new double[np_max];
	yp = new double[np_max];
	zp = new double[np_max];
	up = new double[np_max];
	vp = new double[np_max];
	wp = new double[np_max];
	TE = new double[np_max];
	rho = new double[np_max];
	pm = new double[np_max];
	cs = new double[np_max];
	p = new double[np_max];
	iflag = new int[np_max];
	bflag = new int[np_max];
	xpo = new double[np_max];
	ypo = new double[np_max];
	zpo = new double[np_max];
	upo = new double[np_max];
	vpo = new double[np_max];
	wpo = new double[np_max];
	TEo = new double[np_max];
	rhoo = new double[np_max];
	xpdot = new double[np_max];
	ypdot = new double[np_max];
	zpdot = new double[np_max];
	updot = new double[np_max];
	vpdot = new double[np_max];
	wpdot = new double[np_max];
	TEdot = new double[np_max];
	rhodot = new double[np_max];
	ax = new double[np_max];
	ay = new double[np_max];
	az = new double[np_max];
	ar = new double[np_max];
	ux = new double[np_max];
	vx = new double[np_max];
	wx = new double[np_max];
	aTE = new double[np_max];
	xcor = new double[np_max];
	ycor = new double[np_max];
	zcor = new double[np_max];
	A_MLS = new double**[np_max];
	for (size_t i = 0; i < np_max; ++i)
	{
		A_MLS[i] = new double*[4];
		for (int j = 0; j < 4; ++j)
			A_MLS[i][j] = new double[4];
	}
	rho_sum_MLS = new double[np_max];
	sum_wab = new double[np_max];
	rho_sum = new double[np_max];
	beta0_MLS = new double[np_max];
	beta1x_MLS = new double[np_max];
	beta1y_MLS = new double[np_max];
	beta1z_MLS = new double[np_max];
	for (size_t i=0;i<np;i++)
	{
		xp[i]=0;
		yp[i]=0;
		zp[i]=0;
		up[i]=0;
		vp[i]=0;
		wp[i]=0;
		TE[i]=0;
		rho[i]=0;
		pm[i]=0;
		cs[i]=0;
		p[i]=0;
		iflag[i] = 1;
		bflag[i] = 0;
		xpo[i]=0;
		ypo[i]=0;
		zpo[i]=0;
		upo[i]=0;
		vpo[i]=0;
		wpo[i]=0;
		TEo[i]=0;
		rhoo[i]=0;
		xpdot[i] = 0;
		ypdot[i] = 0;
		zpdot[i] = 0;
		updot[i] = 0;
		vpdot[i] = 0;
		wpdot[i] = 0;
		TEdot[i] = 0;
		rhodot[i] = 0;
		xpdot[i]=0;
		ypdot[i]=0;
		zpdot[i]=0;
		updot[i]=0;
		vpdot[i]=0;
		wpdot[i]=0;
		TEdot[i]=0;
		rhodot[i]=0;
		ax[i]=0;
		ay[i]=0;
		az[i]=0;
		ar[i]=0;
		ux[i]=0;
		vx[i]=0;
		wx[i]=0;
		aTE[i]=0;
		xcor[i]=0;
		ycor[i]=0;
		zcor[i]=0;
	}
	int n0,n1,m0,m1,l0,l1;
	if ((localvlx[0] - vlx[0]< 0) || (localvlx[0] == vlx[0]))
		n0 = 0;
	else
		n0 = (int)((localvlx[0]-vlx[0])/dx) + 1;
	//	if (local_ncn_interior[0] != 0)
	//		n0--;
	n1 = (int)((min(vlx[1],localvlx[1])+(1e-10)*dx-vlx[0])/dx);
	if (local_ncn_interior[1] != local_ncn)
		n1++;
	if (localvly[0] - vly[0]< 0)
		m0 = 0;
	else
		m0 = (int)((localvly[0]-vly[0])/dy) + 1;
	m1 = (int)((min(vly[1],localvly[1])+(1e-10)*dy-vly[0])/dy);
	if (localvlz[0] - vlz[0]< 0)
		l0 = 0;
	else
		l0 = (int)((localvlz[0]-vlz[0])/dz) + 1;
	l1 = (int)((min(vlz[1],localvlz[1])+(1e-10)*dz-vlz[0])/dz);
	scalar = new double[(n1-n0+1) * (m1-m0+1) * (l1-l0+1)];
}

void STORAGE::GetB(double x,double y,double z,double t,double B[])
{
	double zz = z + 6 * t;
	B[0] = 120*sqrt(0.5 * (1.0 - tanh((zz - 9.5)/0.62)));
	B[1] = 0;
	B[2] = 0;
}

void STORAGE::GetUxB(double UxB[],size_t i)
{
	double vel[3]={up[i],vp[i],wp[i]};
	double v_B[3] = {0,0,-6};
	vel[0] = vel[0] - v_B[0];
	vel[1] = vel[1] - v_B[1];
	vel[2] = vel[2] - v_B[2];
	double B[3];
	GetB(xp[i],yp[i],zp[i],t,B);
	CrossProduct(vel,B,UxB);
}

void STORAGE::CrossProduct(double a[],double b[],double c[])
{
	c[0] = a[1]*b[2]-a[2]*b[1];
	c[1] = -a[0]*b[2]-a[2]*b[0];
	c[2] = a[0]*b[1]-a[1]*b[0];
}

void STORAGE::ParallelInitialization()
{
	if (dim != 3)
		assert(0);
	if (dim == 3)
	{
		PartitionIndex[2] = id / (PN * PM);
		PartitionIndex[1] = (id - PartitionIndex[2] * (PN * PM)) / PN;
		PartitionIndex[0] = id - PartitionIndex[2] * (PN * PM) - PartitionIndex[1] * PN;
	}
	else if (dim == 2)
	{
		PartitionIndex[2] = id / PN;
		PartitionIndex[0] = id - PartitionIndex[2] * PN;
	}
	if (IBC != 3)
		ncn = int(((N+3)*dx)/(2*h))+1;
	else
		ncn = int((N*dx)/(2*h))+1;
	if (dim == 3)
	{
		if (IBC != 3)
			ncm = int(((M+3)*dx)/(2*h))+1;
		else
			ncm = int((M*dx)/(2*h))+1;
	}
	else
		ncm = 1;
	if (IBC != 3)
		ncl= int(((L+3)*dx)/(2*h))+1;
	else
		ncl = int((L*dx)/(2*h))+1;
	if (IBC !=3)
	{
		zmax = ncl*2*h + vlz[0] - 1.5*dz;
	}
	if (PartitionIndex[0] + 1 != PN)
	{
		if (ncn % PN == 0)
			local_ncn = ncn / PN;
		else
			local_ncn = ncn / PN + 1;
		if (IBC == 3)
		{
			localvlx[0] = PartitionIndex[0] * local_ncn * 2 * h;
			localvlx[1] = (PartitionIndex[0] + 1) * local_ncn * 2 * h;
		}
		else
		{
			localvlx[0] = vlx[0] - 1.5 * dx + PartitionIndex[0] * local_ncn * 2 * h;
			localvlx[1] = vlx[0] - 1.5 * dx + (PartitionIndex[0] + 1) * local_ncn * 2 * h;
		}
	}
	else
	{
		if (ncn % PN == 0)
		{
			local_ncn = ncn / PN;
			if (IBC == 3)
			{
				localvlx[0] = PartitionIndex[0] * local_ncn * 2 * h;
				localvlx[1] = vlx[1];
			}
			else
			{
				localvlx[0] = vlx[0] - 1.5 * dx + PartitionIndex[0] * local_ncn * 2 * h;
				localvlx[1] = vlx[1] + 1.5 * dx;
			}
		}
		else
		{
			local_ncn = ncn - PartitionIndex[0] * (ncn / PN + 1);
			if (IBC == 3)
			{
				localvlx[0] = vlx[0] + PartitionIndex[0] * (ncn / PN + 1) * 2 * h;
				localvlx[1] = vlx[1];
			}
			else
			{
				localvlx[0] = vlx[0] - 1.5 * dx + PartitionIndex[0] * (ncn / PN + 1) * 2 * h;
				localvlx[1] = vlx[1] + 1.5 * dx;
			}
		}
	}
	if (PartitionIndex[1] + 1 != PM)
	{
		if (ncm % PM == 0)
			local_ncm = ncm / PM;
		else
			local_ncm = ncm / PM + 1;
		if (IBC == 3)
		{
			localvly[0] = PartitionIndex[1] * local_ncm * 2 * h;
			localvly[1] = (PartitionIndex[1] + 1) * local_ncm * 2 * h;
		}
		else
		{
			localvly[0] = vly[0] - 1.5 * dy + PartitionIndex[1] * local_ncm * 2 * h;
			localvly[1] = vly[0] - 1.5 * dy + (PartitionIndex[1] + 1) * local_ncm * 2 * h;
		}
	}
	else
	{
		if (ncm % PM == 0)
		{
			local_ncm = ncm / PM;
			if (IBC == 3)
			{
				localvly[0] = PartitionIndex[1] * local_ncm * 2 * h;
				localvly[1] = vly[1];
			}
			else
			{
				localvly[0] = vly[0] - 1.5 * dy + PartitionIndex[1] * local_ncm * 2 * h;
				localvly[1] = vly[1] + 1.5 * dy;
			}
		}
		else
		{
			local_ncm = ncm - PartitionIndex[1] * (ncm / PM + 1);
			if (IBC == 3)
			{
				localvly[0] = PartitionIndex[1] * (ncm / PM + 1) * 2 * h;
				localvly[1] = vly[1];
			}
			else
			{
				localvly[0] = vly[0] - 1.5 * dy +  PartitionIndex[1] * (ncm / PM + 1) * 2 * h;
				localvly[1] = vly[1] + 1.5 * dy;
			}
		}
	}
	if (PartitionIndex[2] + 1 != PL)
	{
		if (ncl % PL == 0)
			local_ncl = ncl / PL;
		else
			local_ncl = ncl / PL + 1;
		if (IBC == 3)
		{
			localvlz[0] = PartitionIndex[2] * local_ncl * 2 * h;
			localvlz[1] = (PartitionIndex[2] + 1) * local_ncl * 2 * h;
		}
		else
		{
			localvlz[0] = vlz[0] - 1.5 * dz +  PartitionIndex[2] * local_ncl * 2 * h;
			localvlz[1] = vlz[0] - 1.5 * dz +  (PartitionIndex[2] + 1) * local_ncl * 2 * h;
		}
	}
	else
	{
		if (ncl % PL == 0)
		{
			local_ncl = ncl / PL;
			if (IBC == 3)
			{
				localvlz[0] = PartitionIndex[2] * local_ncl * 2 * h;
				localvlz[1] = vlz[1];
			}
			else
			{
				localvlz[0] = vlz[0] - 1.5 * dz +  PartitionIndex[2] * local_ncl * 2 * h;
				localvlz[1] = vlz[1] + 1.5 * dz;
			}
		}
		else
		{
			local_ncl = ncl - PartitionIndex[2] * (ncl / PL + 1);
			if (IBC == 3)
			{
				localvlz[0] = PartitionIndex[2] * (ncl / PL + 1) * 2 * h;
				localvlz[1] = vlz[1];
			}
			else
			{
				localvlz[0] = vlz[0] - 1.5 * dz +  PartitionIndex[2] * (ncl / PL + 1) * 2 * h;
				localvlz[1] = vlz[1] + 1.5 * dz;
			}
		}
	}
	xmin = localvlx[0];
	ymin = localvly[0];
	zmin = localvlz[0];
}

void STORAGE::ParallelAdjustment()
{
	if (dim != 3)
		assert(0);
	local_ncn_interior[0] = 0;
	local_ncn_interior[1] = local_ncn;
	local_ncm_interior[0] = 0;
	local_ncm_interior[1] = local_ncm;
	local_ncl_interior[0] = 0;
	local_ncl_interior[1] = local_ncl;
	if (PartitionIndex[0] != 0)
	{
		//localvlx[0] -= 2 *h;
		++local_ncn;
		local_ncn_interior[0] = 1;
		++local_ncn_interior[1];
	}
	if (PartitionIndex[0] != PN - 1)
	{
		//localvlx[1] += 2 *h;
		++local_ncn;
	}
	if (PartitionIndex[1] != 0)
	{
		//localvly[0] -= 2 *h;
		++local_ncm;
		local_ncm_interior[0] = 1;
		++local_ncm_interior[1];
	}
	if (PartitionIndex[1] != PM - 1)
	{
		//localvly[1] += 2 *h;
		++local_ncm;
	}
	if (PartitionIndex[2] != 0)
	{
		//localvlz[0] -= 2 *h;
		++local_ncl;
		local_ncl_interior[0] = 1;
		++local_ncl_interior[1];
	}
	if (PartitionIndex[2] != PL - 1)
	{
		//localvlz[1] += 2 *h;
		++local_ncl;
	}
	nx = local_ncn * local_ncm * local_ncl;
	zmax = vlz[1];
	if (PartitionIndex[0] != 0)
		xmin -= 2 * h;
	if (PartitionIndex[1] != 0)
		ymin -= 2 * h;
	if (PartitionIndex[2] != 0)
		zmin -= 2 * h;
}

void STORAGE::SendRecvBuffer()
{
	divide(nb,np,2);
	if (PN * PM * PL != 1)
	{
		DisplacedAndBufferParticlexp.clear();
		DisplacedAndBufferParticleyp.clear();
		DisplacedAndBufferParticlezp.clear();
		DisplacedAndBufferParticleup.clear();
		DisplacedAndBufferParticlevp.clear();
		DisplacedAndBufferParticlewp.clear();
		DisplacedAndBufferParticleTE.clear();
		DisplacedAndBufferParticlerho.clear();
		DisplacedAndBufferParticlepm.clear();
		DisplacedAndBufferParticle.clear();
		int IndexDisplacedAndBuffer = 0;
		TotalNumberOfParticleSend = TotalNumberOfParticleRecv = 0;
		size_t NumOfParticleSend, NumOfParticleRecv;
		size_t DestinationID;
		size_t SourceID;
		MPI_Status status;
		//26 Neighbours
		//East-North-Up:25
		MPI_Barrier(MPI_COMM_WORLD);
		if (PN != 1 && PM != 1 && PL != 1)
		{
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != PM -1 && PartitionIndex[2] != PL -1)
			{
				NumOfParticleSend = 0;
				IndexDisplacedAndBuffer = nx-1;
				vector<size_t>::iterator it = ibox_k2[IndexDisplacedAndBuffer].begin();
				if (nc_k2[IndexDisplacedAndBuffer] != 0)
				{
					for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
					{
						size_t num;
						num=*it;
						DisplacedAndBufferParticle.push_back(num); TotalNumberOfParticleSend++;bflag[num] = 1;
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 25, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 25, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != PM -1 && PartitionIndex[2] != PL -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 25, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 25, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//East-North-Mid:24
		if (PN != 1 && PM != 1)
		{
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != PM -1)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int k = local_ncl_interior[0]; k<local_ncl_interior[1]; k++)
				{

					IndexDisplacedAndBuffer = Index3D(local_ncn-1,local_ncm-1,k,local_ncn,local_ncm,local_ncl);
					it = ibox_k2[IndexDisplacedAndBuffer].begin();
					if (nc_k2[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num); TotalNumberOfParticleSend++;bflag[num] = 1;
						}
					}

				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 24, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0)
			{
				SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 24, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != PM -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 24, MPI_COMM_WORLD);
				}
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 24, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//East-North-Down:23
		if (PN != 1 && PM != 1 && PL != 1)
		{
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != PM -1 && PartitionIndex[2] != 0)
			{
				NumOfParticleSend = 0;
				IndexDisplacedAndBuffer = Index3D(local_ncn-1,local_ncm-1,0,local_ncn,local_ncm,local_ncl);
				vector<size_t>::iterator it = ibox_k2[IndexDisplacedAndBuffer].begin();
				if (nc_k2[IndexDisplacedAndBuffer] != 0)
				{
					for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
					{
						size_t num;
						num=*it;
						DisplacedAndBufferParticle.push_back(num); TotalNumberOfParticleSend++;bflag[num] = 1;
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 23, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0 && PartitionIndex[2] != PL - 1)
			{
				SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2] + 1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 23, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != PM -1 && PartitionIndex[2] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 23, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0 && PartitionIndex[2] != PL - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 23, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//East-Mid-UP:22
		if (PN != 1 && PL != 1)
		{
			if (PartitionIndex[0] != PN -1 && PartitionIndex[2] != PL -1)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int j = local_ncm_interior[0];j<local_ncm_interior[1]; j++)
				{

					IndexDisplacedAndBuffer = Index3D(local_ncn-1,j,local_ncl-1,local_ncn,local_ncm,local_ncl);
					it = ibox_k2[IndexDisplacedAndBuffer].begin();
					if (nc_k2[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num); TotalNumberOfParticleSend++;bflag[num] = 1;
						}
					}

				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 22, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[2] != 0)
			{
				SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 22, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != PN -1 && PartitionIndex[2] != PL -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 22, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[2] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 22, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//East-Mid-Mid:21
		if (PN != 1)
		{
			if (PartitionIndex[0] != PN -1)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int j = local_ncm_interior[0]; j<local_ncm_interior[1];j++)
					for (int k = local_ncl_interior[0]; k<local_ncl_interior[1]; k++)
					{

						IndexDisplacedAndBuffer = Index3D(local_ncn-1,j,k,local_ncn,local_ncm,local_ncl);
						it = ibox_k2[IndexDisplacedAndBuffer].begin();
						if (nc_k2[IndexDisplacedAndBuffer] != 0)
						{
							for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
							{
								size_t num;
								num=*it;
								DisplacedAndBufferParticle.push_back(num); TotalNumberOfParticleSend++;bflag[num] = 1;
							}
						}

					}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2],PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 21, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != 0)
			{
				SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2],PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 21, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != PN -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2],PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 21, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
					//			BufferParticle.clear();
				}
			}
			if (PartitionIndex[0] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2],PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 21, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//East-Mid-Down:20
		if (PN!= 1 && PL!=1)
		{
			if (PartitionIndex[0] != PN -1 && PartitionIndex[2] != 0)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int j = local_ncm_interior[0];j<local_ncm_interior[1]; j++)
				{

					IndexDisplacedAndBuffer = Index3D(local_ncn-1,j,0,local_ncn,local_ncm,local_ncl);
					it = ibox_k2[IndexDisplacedAndBuffer].begin();
					if (nc_k2[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num); TotalNumberOfParticleSend++;bflag[num] = 1;
						}
					}

				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 20, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[2] != PL - 1)
			{
				SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 20, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != PN -1 && PartitionIndex[2] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 20, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[2] != PL - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 20, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//East-South-Up:19
		if (PN != 1 && PM != 1 && PL != 1)
		{
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != 0 && PartitionIndex[2] != PL -1)
			{
				NumOfParticleSend = 0;
				IndexDisplacedAndBuffer = Index3D(local_ncn-1,0,local_ncl-1,local_ncn,local_ncm,local_ncl);
				vector<size_t>::iterator it = ibox_k2[IndexDisplacedAndBuffer].begin();
				if (nc_k2[IndexDisplacedAndBuffer] != 0)
				{
					for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
					{
						size_t num;
						num=*it;
						DisplacedAndBufferParticle.push_back(num); TotalNumberOfParticleSend++;bflag[num] = 1;
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2]+1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 19, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != 0)
			{
				SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 19, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != 0 && PartitionIndex[2] != PL -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 19, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
					//			BufferParticle.clear();
				}
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 19, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//East-South-Mid:18
		if (PN != 1 && PM != 1)
		{
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != 0)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int k = local_ncl_interior[0]; k<local_ncl_interior[1]; k++)
				{

					IndexDisplacedAndBuffer = Index3D(local_ncn-1,0 ,k,local_ncn,local_ncm,local_ncl);
					it = ibox_k2[IndexDisplacedAndBuffer].begin();
					if (nc_k2[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num); TotalNumberOfParticleSend++;bflag[num] = 1;
						}
					}

				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 18, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1)
			{
				SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 18, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 18, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
					//			BufferParticle.clear();
				}
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 18, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//East-South-Down:17
		if (PN != 1 && PM != 1 && PL != 1)
		{
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				NumOfParticleSend = 0;
				IndexDisplacedAndBuffer = Index3D(local_ncn-1,0,0,local_ncn,local_ncm,local_ncl);
				vector<size_t>::iterator it = ibox_k2[IndexDisplacedAndBuffer].begin();
				if (nc_k2[IndexDisplacedAndBuffer] != 0)
				{
					for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
					{
						size_t num;
						num=*it;
						DisplacedAndBufferParticle.push_back(num); TotalNumberOfParticleSend++;bflag[num] = 1;
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 17, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != PL - 1)
			{
				SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2] + 1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 17, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 17, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
					//			BufferParticle.clear();
				}
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != PL - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 117, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
		}
		//Mid-North-Up:16
		if (PM != 1 && PL != 1)
		{
			if (PartitionIndex[1] != PM - 1 && PartitionIndex[2] != PL -1)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int i = local_ncn_interior[0]; i<local_ncn_interior[1]; i++)
				{

					IndexDisplacedAndBuffer = Index3D(i,local_ncm-1 ,local_ncl-1,local_ncn,local_ncm,local_ncl);
					it = ibox_k2[IndexDisplacedAndBuffer].begin();
					if (nc_k2[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num); TotalNumberOfParticleSend++;bflag[num] = 1;
						}
					}

				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 16, MPI_COMM_WORLD);
			}
			if (PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 16, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[1] != PM - 1 && PartitionIndex[2] != PL -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 16, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 16, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//Mid-North-Mid:15
		if (PM !=1)
		{
			if (PartitionIndex[1] != PM - 1)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int i = local_ncn_interior[0]; i<local_ncn_interior[1];i++)
					for (int k = local_ncl_interior[0]; k<local_ncl_interior[1]; k++)
					{

						IndexDisplacedAndBuffer = Index3D(i,local_ncm-1,k,local_ncn,local_ncm,local_ncl);
						it = ibox_k2[IndexDisplacedAndBuffer].begin();
						if (nc_k2[IndexDisplacedAndBuffer] != 0)
						{
							for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
							{
								size_t num;
								num=*it;
								DisplacedAndBufferParticle.push_back(num); TotalNumberOfParticleSend++;bflag[num] = 1;
							}
						}

					}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 15, MPI_COMM_WORLD);
			}
			if (PartitionIndex[1] != 0)
			{
				SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 15, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[1] != PM - 1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 15, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[1] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 15, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//Mid-North-Down:14
		if (PM != 1 && PL != 1)
		{
			if (PartitionIndex[1] != PM - 1 && PartitionIndex[2] != 0)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int i = local_ncn_interior[0]; i<local_ncn_interior[1]; i++)
				{

					IndexDisplacedAndBuffer = Index3D(i,local_ncm-1,0,local_ncn,local_ncm,local_ncl);
					it = ibox_k2[IndexDisplacedAndBuffer].begin();
					if (nc_k2[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num); TotalNumberOfParticleSend++;bflag[num] = 1;
						}
					}

				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 14, MPI_COMM_WORLD);
			}
			if (PartitionIndex[1] != 0 && PartitionIndex[2] != PL - 1)
			{
				SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2]+1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 14, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[1] != PM - 1 && PartitionIndex[2] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 14, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[1] != 0 && PartitionIndex[2] != PL - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 14, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//Mid-Mid-UP:13
		if (PL != 1)
		{
			if (PartitionIndex[2] != PL -1)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int i = local_ncn_interior[0]; i<local_ncn_interior[1];i++)
					for (int j = local_ncm_interior[0];j<local_ncm_interior[1]; j++)
					{

						IndexDisplacedAndBuffer = Index3D(i,j,local_ncl-1,local_ncn,local_ncm,local_ncl);
						it = ibox_k2[IndexDisplacedAndBuffer].begin();
						if (nc_k2[IndexDisplacedAndBuffer] != 0)
						{
							for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
							{
								size_t num;
								num=*it;
								DisplacedAndBufferParticle.push_back(num); TotalNumberOfParticleSend++;bflag[num] = 1;
							}
						}

					}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 13, MPI_COMM_WORLD);
			}
			if (PartitionIndex[2] != 0)
			{
				SourceID = Index3D(PartitionIndex[0],PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 13, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[2] != PL -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 13, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
					//			BufferParticle.clear();
				}
			}
			if (PartitionIndex[2] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0],PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 13, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//Mid-Mid-Down:12
		if (PL != 1)
		{
			if (PartitionIndex[2] != 0)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int i = local_ncn_interior[0]; i<local_ncn_interior[1];i++)
					for (int j = local_ncm_interior[0];j<local_ncm_interior[1]; j++)
					{

						IndexDisplacedAndBuffer = Index3D(i,j,0,local_ncn,local_ncm,local_ncl);
						it = ibox_k2[IndexDisplacedAndBuffer].begin();
						if (nc_k2[IndexDisplacedAndBuffer] != 0)
						{
							for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
							{
								size_t num;
								num=*it;
								DisplacedAndBufferParticle.push_back(num); TotalNumberOfParticleSend++;bflag[num] = 1;
							}
						}

					}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 12, MPI_COMM_WORLD);
			}
			if (PartitionIndex[2] != PL - 1)
			{
				SourceID = Index3D(PartitionIndex[0],PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 12, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[2] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 12, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[2] != PL - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0],PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 12, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//Mid-South-Up:11
		if (PM != 1 && PL != 1)
		{
			if (PartitionIndex[1] != 0 && PartitionIndex[2] != PL -1)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int i = local_ncn_interior[0]; i<local_ncn_interior[1]; i++)
				{

					IndexDisplacedAndBuffer = Index3D(i,0,local_ncl-1,local_ncn,local_ncm,local_ncl);
					it = ibox_k2[IndexDisplacedAndBuffer].begin();
					if (nc_k2[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num); TotalNumberOfParticleSend++;bflag[num] = 1;
						}
					}

				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2]+1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 11, MPI_COMM_WORLD);
			}
			if (PartitionIndex[1] != PM - 1 && PartitionIndex[2] != 0)
			{
				SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 11, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[1] != 0 && PartitionIndex[2] != PL -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 11, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[1] != PM - 1 && PartitionIndex[2] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 11, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//Mid-South-Mid:10
		if (PM != 1)
		{
			if (PartitionIndex[1] != 0)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int i = local_ncn_interior[0]; i<local_ncn_interior[1];i++)
					for (int k = local_ncl_interior[0]; k<local_ncl_interior[1]; k++)
					{

						IndexDisplacedAndBuffer = Index3D(i,0,k,local_ncn,local_ncm,local_ncl);
						it = ibox_k2[IndexDisplacedAndBuffer].begin();
						if (nc_k2[IndexDisplacedAndBuffer] != 0)
						{
							for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
							{
								size_t num;
								num=*it;
								DisplacedAndBufferParticle.push_back(num); TotalNumberOfParticleSend++;bflag[num] = 1;
							}
						}

					}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 10, MPI_COMM_WORLD);
			}
			if (PartitionIndex[1] != PM - 1)
			{
				SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 10, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[1] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 10, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[1] != PM - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 10, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//Mid-South-Down:9
		if (PM != 1 && PL != 1)
		{
			if (PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int i = local_ncn_interior[0]; i<local_ncn_interior[1]; i++)
				{

					IndexDisplacedAndBuffer = Index3D(i,0,0,local_ncn,local_ncm,local_ncl);
					it = ibox_k2[IndexDisplacedAndBuffer].begin();
					if (nc_k2[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num); TotalNumberOfParticleSend++;bflag[num] = 1;
						}
					}

				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 19, MPI_COMM_WORLD);
			}
			if (PartitionIndex[1] != PM - 1 && PartitionIndex[2] != PL - 1)
			{
				SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 9, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 9, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[1] != PM - 1 && PartitionIndex[2] != PL - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 9, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//West-North-Up:8
		if (PN != 1 && PM != 1 && PL != 1)
		{
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != PL -1)
			{
				NumOfParticleSend = 0;
				IndexDisplacedAndBuffer = Index3D(0,local_ncm-1,local_ncl-1,local_ncn,local_ncm,local_ncl);
				vector<size_t>::iterator it = ibox_k2[IndexDisplacedAndBuffer].begin();
				if (nc_k2[IndexDisplacedAndBuffer] != 0)
				{
					for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
					{
						size_t num;
						num=*it;
						DisplacedAndBufferParticle.push_back(num); TotalNumberOfParticleSend++;bflag[num] = 1;
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 8, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 8, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != PL -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 8, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 8, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//West-North-Mid:7
		if (PN != 1 && PM != 1)
		{
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int k = local_ncl_interior[0]; k<local_ncl_interior[1]; k++)
				{

					IndexDisplacedAndBuffer = Index3D(0,local_ncm-1,k,local_ncn,local_ncm,local_ncl);
					it = ibox_k2[IndexDisplacedAndBuffer].begin();
					if (nc_k2[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num); TotalNumberOfParticleSend++;bflag[num] = 1;
						}
					}

				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 7, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != 0)
			{
				SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 7, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 7, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 7, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//West-North-Down:6
		if (PN != 1 && PM != 1 && PL != 1)
		{
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != 0)
			{
				NumOfParticleSend = 0;
				IndexDisplacedAndBuffer = Index3D(0,local_ncm-1,0,local_ncn,local_ncm,local_ncl);
				vector<size_t>::iterator it = ibox_k2[IndexDisplacedAndBuffer].begin();
				if (nc_k2[IndexDisplacedAndBuffer] != 0)
				{
					for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
					{
						size_t num;
						num=*it;
						DisplacedAndBufferParticle.push_back(num); TotalNumberOfParticleSend++;bflag[num] = 1;
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 6, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != 0 && PartitionIndex[2] != PL - 1)
			{
				SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2] + 1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 6, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 6, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != 0 && PartitionIndex[2] != PL - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 6, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//West-Mid-UP:5
		if (PN != 1 && PL != 1)
		{
			if (PartitionIndex[0] != 0 && PartitionIndex[2] != PL -1)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int j = local_ncm_interior[0];j<local_ncm_interior[1]; j++)
				{

					IndexDisplacedAndBuffer = Index3D(0,j,local_ncl-1,local_ncn,local_ncm,local_ncl);
					it = ibox_k2[IndexDisplacedAndBuffer].begin();
					if (nc_k2[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num); TotalNumberOfParticleSend++;bflag[num] = 1;
						}
					}

				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 5, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[2] != 0)
			{
				SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 5, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != 0 && PartitionIndex[2] != PL -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 5, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[2] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 5, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//West-Mid-Mid:4
		if (PN != 1)
		{
			if (PartitionIndex[0] != 0)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int j = local_ncm_interior[0]; j<local_ncm_interior[1];j++)
					for (int k = local_ncl_interior[0]; k<local_ncl_interior[1]; k++)
					{

						IndexDisplacedAndBuffer = Index3D(0,j,k,local_ncn,local_ncm,local_ncl);
						it = ibox_k2[IndexDisplacedAndBuffer].begin();
						if (nc_k2[IndexDisplacedAndBuffer] != 0)
						{
							for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
							{
								size_t num;
								num=*it;
								DisplacedAndBufferParticle.push_back(num); TotalNumberOfParticleSend++;bflag[num] = 1;
							}
						}

					}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2],PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 4, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != PN - 1)
			{
				SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2],PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 4, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2],PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 4, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != PN - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2],PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 4, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//West-Mid-Down:3
		if (PN!= 1 && PL!=1)
		{
			if (PartitionIndex[0] != 0 && PartitionIndex[2] != 0)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int j = local_ncm_interior[0];j<local_ncm_interior[1]; j++)
				{

					IndexDisplacedAndBuffer = Index3D(0,j,0,local_ncn,local_ncm,local_ncl);
					it = ibox_k2[IndexDisplacedAndBuffer].begin();
					if (nc_k2[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num); TotalNumberOfParticleSend++;bflag[num] = 1;
						}
					}

				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 3, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[2] != PL - 1)
			{
				SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 3, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != 0 && PartitionIndex[2] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 3, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[2] != PL - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 3, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//West-South-Up:2
		if (PN != 1 && PM != 1 && PL != 1)
		{
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0 && PartitionIndex[2] != PL -1)
			{
				NumOfParticleSend = 0;
				IndexDisplacedAndBuffer = Index3D(0,0,local_ncl-1,local_ncn,local_ncm,local_ncl);
				vector<size_t>::iterator it = ibox_k2[IndexDisplacedAndBuffer].begin();
				if (nc_k2[IndexDisplacedAndBuffer] != 0)
				{
					for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
					{
						size_t num;
						num=*it;
						DisplacedAndBufferParticle.push_back(num); TotalNumberOfParticleSend++;bflag[num] = 1;
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2]+1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 2, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != 0)
			{
				SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 2, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0 && PartitionIndex[2] != PL -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 2, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 2, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//West-South-Mid:1
		if (PN != 1 && PM != 1)
		{
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int k = local_ncl_interior[0]; k<local_ncl_interior[1]; k++)
				{

					IndexDisplacedAndBuffer = Index3D(0,0,k,local_ncn,local_ncm,local_ncl);
					it = ibox_k2[IndexDisplacedAndBuffer].begin();
					if (nc_k2[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num); TotalNumberOfParticleSend++;bflag[num] = 1;
						}
					}

				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 1, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != PM - 1)
			{
				SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 1, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 1, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != PM - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 1, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//West-South-Down:0
		if (PN != 1 && PM != 1 && PL != 1)
		{
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				NumOfParticleSend = 0;
				IndexDisplacedAndBuffer = Index3D(0,0,0,local_ncn,local_ncm,local_ncl);
				vector<size_t>::iterator it = ibox_k2[IndexDisplacedAndBuffer].begin();
				if (nc_k2[IndexDisplacedAndBuffer] != 0)
				{
					for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
					{
						size_t num;
						num=*it;
						DisplacedAndBufferParticle.push_back(num); TotalNumberOfParticleSend++;bflag[num] = 1;
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 0, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != PL - 1)
			{
				SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 0, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 0, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != PL - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 0, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//Finishing DisplacedAndBuffer particle transfer
		if (TotalNumberOfParticleRecv != 0 || TotalNumberOfParticleSend != 0)
		{
			if (((i_SR_step == 2) && (TotalNumberOfParticleRecv != DisplacedAndBufferParticlerho.size())) ||
					((i_SR_step == 1) && (2 * TotalNumberOfParticleRecv != DisplacedAndBufferParticlerho.size())) )
				assert(0);
			size_t nptemp = np + TotalNumberOfParticleRecv - TotalNumberOfParticleSend;
			size_t np_max_temp = np_max;
			if (nptemp >= np_max_temp)
			{
				delete[] p;
				delete[] xpdot;
				delete[] ypdot;
				delete[] zpdot;
				delete[] updot;
				delete[] vpdot;
				delete[] wpdot;
				delete[] rhodot;
				delete[] ax;
				delete[] ay;
				delete[] az;
				delete[] ar;
				delete[] ux;
				delete[] vx;
				delete[] wx;
				delete[] aTE;
				delete[] xcor;
				delete[] ycor;
				delete[] zcor;
				delete [] cs;
				for (size_t i = 0; i < np_max; ++i)
				{
					for (int j = 0; j < 4; ++j)
						delete [] A_MLS[i][j];
					delete [] A_MLS[i];
				}
				delete [] A_MLS;
				delete[] rho_sum_MLS;
				delete[] sum_wab;
				delete[] rho_sum;
				delete[] beta0_MLS;
				delete[] beta1x_MLS;
				delete[] beta1y_MLS;
				delete[] beta1z_MLS;
				np_max = nptemp + nptemp / 5; //new np_max
				p = new double[np_max];
				xpdot = new double[np_max];
				ypdot = new double[np_max];
				zpdot = new double[np_max];
				updot = new double[np_max];
				vpdot = new double[np_max];
				wpdot = new double[np_max];
				TEdot = new double[np_max];
				rhodot = new double[np_max];
				ax = new double[np_max];
				ay = new double[np_max];
				az = new double[np_max];
				ar = new double[np_max];
				ux = new double[np_max];
				vx = new double[np_max];
				wx = new double[np_max];
				aTE = new double[np_max];
				xcor = new double[np_max];
				ycor = new double[np_max];
				zcor = new double[np_max];
				cs = new double[np_max];
				A_MLS = new double**[np_max];
				for (size_t i = 0; i < np_max; ++i)
				{
					A_MLS[i] = new double*[4];
					for (int j = 0; j < 4; ++j)
						A_MLS[i][j] = new double[4];
				}
				rho_sum_MLS = new double[np_max];
				sum_wab = new double[np_max];
				rho_sum = new double[np_max];
				beta0_MLS = new double[np_max];
				beta1x_MLS = new double[np_max];
				beta1y_MLS = new double[np_max];
				beta1z_MLS = new double[np_max];

			}
			//			double * xpt = new double[nptemp];
			//			double * ypt = new double[nptemp];
			//			double * zpt = new double[nptemp];
			//			double * upt = new double[nptemp];
			//			double * vpt = new double[nptemp];
			//			double * wpt = new double[nptemp];
			//			double * pt = new double[nptemp];
			//			double * rhot = new double[nptemp];
			//			double * pmt = new double[nptemp];
			//			double * xpot=0;
			//			double * ypot=0;
			//			double * zpot=0;
			//			double * upot=0;
			//			double * vpot=0;
			//			double * wpot=0;
			//			double * pot=0;
			//			double * rhoot=0;
			//			if (i_SR_step == 1)
			//			{
			//				xpot = new double[nptemp];
			//				ypot = new double[nptemp];
			//				zpot = new double[nptemp];
			//				upot = new double[nptemp];
			//				vpot = new double[nptemp];
			//				wpot = new double[nptemp];
			//				pot = new double[nptemp];
			//				rhoot = new double[nptemp];
			//			}
			size_t count = nb;
			for (size_t i = 0; i<nb; i++)
			{
				xpdot[i] = xp[i];
				ypdot[i] = yp[i];
				zpdot[i] = zp[i];
				updot[i] = up[i];
				vpdot[i] = vp[i];
				wpdot[i] = wp[i];
				TEdot[i] = TE[i];
				rhodot[i] = rho[i];
				xcor[i] = pm[i];
				if (i_SR_step == 1)
				{
					ax[i] = xpo[i];
					ay[i] = ypo[i];
					az[i] = zpo[i];
					ux[i] = upo[i];
					vx[i] = vpo[i];
					wx[i] = wpo[i];
					aTE[i] = TEo[i];
					ar[i] = rhoo[i];
				}
			}
			for (size_t i = nb; i<np; i++)
			{
				if (bflag[i] != 1)
				{
					xpdot[count] = xp[i];
					ypdot[count] = yp[i];
					zpdot[count] = zp[i];
					updot[count] = up[i];
					vpdot[count] = vp[i];
					wpdot[count] = wp[i];
					TEdot[count] = TE[i];
					rhodot[count] = rho[i];
					xcor[count] = pm[i];
					if (i_SR_step == 1)
					{
						ax[count] = xpo[i];
						ay[count] = ypo[i];
						az[count] = zpo[i];
						ux[count] = upo[i];
						vx[count] = vpo[i];
						wx[count] = wpo[i];
						aTE[count] = TEo[i];
						ar[count] = rhoo[i];
					}
					++count;
				}
			}
			for (size_t i = 0; i < TotalNumberOfParticleRecv; i++)
			{

				xpdot[count] = (i_SR_step == 1)?DisplacedAndBufferParticlexp[i*2]:DisplacedAndBufferParticlexp[i];
				ypdot[count] = (i_SR_step == 1)?DisplacedAndBufferParticleyp[i*2]:DisplacedAndBufferParticleyp[i];
				zpdot[count] = (i_SR_step == 1)?DisplacedAndBufferParticlezp[i*2]:DisplacedAndBufferParticlezp[i];
				updot[count] = (i_SR_step == 1)?DisplacedAndBufferParticleup[i*2]:DisplacedAndBufferParticleup[i];
				vpdot[count] = (i_SR_step == 1)?DisplacedAndBufferParticlevp[i*2]:DisplacedAndBufferParticlevp[i];
				wpdot[count] = (i_SR_step == 1)?DisplacedAndBufferParticlewp[i*2]:DisplacedAndBufferParticlewp[i];
				TEdot[count] = (i_SR_step == 1)?DisplacedAndBufferParticleTE[i*2]:DisplacedAndBufferParticleTE[i];
				rhodot[count] = (i_SR_step == 1)?DisplacedAndBufferParticlerho[i*2]:DisplacedAndBufferParticlerho[i];
				xcor[count] = DisplacedAndBufferParticlepm[i];
				if (i_SR_step == 1)
				{
					ax[count] = DisplacedAndBufferParticlexp[i*2+1];
					ay[count] = DisplacedAndBufferParticleyp[i*2+1];
					az[count] = DisplacedAndBufferParticlezp[i*2+1];
					ux[count] = DisplacedAndBufferParticleup[i*2+1];
					vx[count] = DisplacedAndBufferParticlevp[i*2+1];
					wx[count] = DisplacedAndBufferParticlewp[i*2+1];
					aTE[count] = DisplacedAndBufferParticleTE[i*2+1];
					ar[count] = DisplacedAndBufferParticlerho[i*2+1];
				}
				++count;
			}
			if (nptemp != count) assert(0);
			np = nptemp;
			if (nptemp >= np_max_temp)
			{
				delete [] xp;
				delete [] yp;
				delete [] zp;
				delete [] up;
				delete [] vp;
				delete [] wp;
				delete [] TE;
				delete [] rho;
				delete [] pm;
				delete [] xpo;
				delete [] ypo;
				delete [] zpo;
				delete [] upo;
				delete [] vpo;
				delete [] wpo;
				delete [] TEo;
				delete [] rhoo;
				delete [] bflag;
				xp = new double[np_max];
				yp = new double[np_max];
				zp = new double[np_max];
				up = new double[np_max];
				vp = new double[np_max];
				wp = new double[np_max];
				TE = new double[np_max];
				rho = new double[np_max];
				pm = new double[np_max];
				xpo = new double[np_max];
				ypo = new double[np_max];
				zpo = new double[np_max];
				upo = new double[np_max];
				vpo = new double[np_max];
				wpo = new double[np_max];
				TEo = new double[np_max];
				rhoo = new double[np_max];
				bflag = new int[np_max];
			}
			for (size_t i = 0; i < np; i++ )
			{
				xp[i] = xpdot[i];
				yp[i] = ypdot[i];
				zp[i] = zpdot[i];
				up[i] = updot[i];
				vp[i] = vpdot[i];
				wp[i] = wpdot[i];
				TE[i] = TEdot[i];
				rho[i] = rhodot[i];
				pm[i] = xcor[i];
				if (i_SR_step == 1)
				{
					xpo[i] = ax[i];
					ypo[i] = ay[i];
					zpo[i] = az[i];
					upo[i] = ux[i];
					vpo[i] = vx[i];
					wpo[i] = wx[i];
					TEo[i] = aTE[i];
					rhoo[i] = ar[i];
				}
				else
				{
					xpo[i] = 0;
					ypo[i] = 0;
					zpo[i] = 0;
					upo[i] = 0;
					vpo[i] = 0;
					wpo[i] = 0;
					TEo[i] = 0;
					rhoo[i] = 0;
				}
			}
			divide(nb,np,2);
			//			delete [] xpt;
			//			delete [] ypt;
			//			delete [] zpt;
			//			delete [] upt;
			//			delete [] vpt;
			//			delete [] wpt;
			//			delete [] pt;
			//			delete [] rhot;
			//			delete [] pmt;
			//			if (i_SR_step == 1)
			//			{
			//				delete [] xpot;
			//				delete [] ypot;
			//				delete [] zpot;
			//				delete [] upot;
			//				delete [] vpot;
			//				delete [] wpot;
			//				delete [] pot;
			//				delete [] rhoot;
			//			}
			DisplacedAndBufferParticlexp.clear();
			DisplacedAndBufferParticleyp.clear();
			DisplacedAndBufferParticlezp.clear();
			DisplacedAndBufferParticleup.clear();
			DisplacedAndBufferParticlevp.clear();
			DisplacedAndBufferParticlewp.clear();
			DisplacedAndBufferParticleTE.clear();
			DisplacedAndBufferParticlerho.clear();
			DisplacedAndBufferParticlepm.clear();
			DisplacedAndBufferParticle.clear();
			TotalNumberOfParticleSend = TotalNumberOfParticleRecv = 0;
		}


		if (PN != 1 && PM != 1 && PL != 1)
		{
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != PM -1 && PartitionIndex[2] != PL -1)
			{
				NumOfParticleSend = 0;
				IndexDisplacedAndBuffer = Index3D(local_ncn-2,local_ncm-2,local_ncl-2,local_ncn,local_ncm,local_ncl);
				vector<size_t>::iterator it = ibox_k1[IndexDisplacedAndBuffer].begin();
				if (nc_k1[IndexDisplacedAndBuffer] != 0)
				{
					for(;it!=ibox_k1[IndexDisplacedAndBuffer].end();++it)
					{
						size_t num;
						num=*it;
						DisplacedAndBufferParticle.push_back(num);
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 25, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 25, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != PM -1 && PartitionIndex[2] != PL -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 25, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 25, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//East-North-Mid:24
		if (PN != 1 && PM != 1)
		{
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != PM -1)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int k = local_ncl_interior[0]; k<local_ncl_interior[1]; k++)
				{

					IndexDisplacedAndBuffer = Index3D(local_ncn-2,local_ncm-2,k,local_ncn,local_ncm,local_ncl);
					it = ibox_k1[IndexDisplacedAndBuffer].begin();
					if (nc_k1[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k1[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num);
						}
					}
				}

				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 24, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0)
			{
				SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 24, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != PM -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 24, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 24, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//East-North-Down:23
		if (PN != 1 && PM != 1 && PL != 1)
		{
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != PM -1 && PartitionIndex[2] != 0)
			{
				NumOfParticleSend = 0;
				IndexDisplacedAndBuffer = Index3D(local_ncn-2,local_ncm-2,1,local_ncn,local_ncm,local_ncl);
				vector<size_t>::iterator it = ibox_k1[IndexDisplacedAndBuffer].begin();
				if (nc_k1[IndexDisplacedAndBuffer] != 0)
				{
					for(;it!=ibox_k1[IndexDisplacedAndBuffer].end();++it)
					{
						size_t num;
						num=*it;
						DisplacedAndBufferParticle.push_back(num);
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 23, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0 && PartitionIndex[2] != PL - 1)
			{
				SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2] + 1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 23, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != PM -1 && PartitionIndex[2] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 23, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0 && PartitionIndex[2] != PL - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 23, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//East-Mid-UP:22
		if (PN != 1 && PL != 1)
		{
			if (PartitionIndex[0] != PN -1 && PartitionIndex[2] != PL -1)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int j = local_ncm_interior[0];j<local_ncm_interior[1]; j++)
				{

					IndexDisplacedAndBuffer = Index3D(local_ncn-2,j,local_ncl-2,local_ncn,local_ncm,local_ncl);
					it = ibox_k1[IndexDisplacedAndBuffer].begin();
					if (nc_k1[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k1[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num);
						}
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 22, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[2] != 0)
			{
				SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 22, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != PN -1 && PartitionIndex[2] != PL -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 22, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
					//			BufferParticle.clear();
				}
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[2] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 22, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//East-Mid-Mid:21
		if (PN != 1)
		{
			if (PartitionIndex[0] != PN -1)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int j = local_ncm_interior[0]; j<local_ncm_interior[1];j++)
					for (int k = local_ncl_interior[0]; k<local_ncl_interior[1]; k++)
					{
						IndexDisplacedAndBuffer = Index3D(local_ncn-2,j,k,local_ncn,local_ncm,local_ncl);
						it = ibox_k1[IndexDisplacedAndBuffer].begin();
						if (nc_k1[IndexDisplacedAndBuffer] != 0)
						{
							for(;it!=ibox_k1[IndexDisplacedAndBuffer].end();++it)
							{
								size_t num;
								num=*it;
								DisplacedAndBufferParticle.push_back(num);
							}
						}
					}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2],PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 21, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != 0)
			{
				SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2],PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 21, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != PN -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2],PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 21, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2],PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 21, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//East-Mid-Down:20
		if (PN!= 1 && PL!=1)
		{
			if (PartitionIndex[0] != PN -1 && PartitionIndex[2] != 0)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int j = local_ncm_interior[0];j<local_ncm_interior[1]; j++)
				{

					IndexDisplacedAndBuffer = Index3D(local_ncn-2,j,1,local_ncn,local_ncm,local_ncl);
					it = ibox_k1[IndexDisplacedAndBuffer].begin();
					if (nc_k1[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k1[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num);
						}
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 20, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[2] != PL - 1)
			{
				SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 20, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != PN -1 && PartitionIndex[2] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 20, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[2] != PL - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 20, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//East-South-Up:19
		if (PN != 1 && PM != 1 && PL != 1)
		{
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != 0 && PartitionIndex[2] != PL -1)
			{
				NumOfParticleSend = 0;
				IndexDisplacedAndBuffer = Index3D(local_ncn-2,1,local_ncl-2,local_ncn,local_ncm,local_ncl);
				vector<size_t>::iterator it = ibox_k1[IndexDisplacedAndBuffer].begin();
				if (nc_k1[IndexDisplacedAndBuffer] != 0)
				{
					for(;it!=ibox_k1[IndexDisplacedAndBuffer].end();++it)
					{
						size_t num;
						num=*it;
						DisplacedAndBufferParticle.push_back(num);
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2]+1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 19, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != 0)
			{
				SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 19, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != 0 && PartitionIndex[2] != PL -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 19, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
					//			BufferParticle.clear();
				}
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 19, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//East-South-Mid:18
		if (PN != 1 && PM != 1)
		{
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != 0)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int k = local_ncl_interior[0]; k<local_ncl_interior[1]; k++)
				{

					IndexDisplacedAndBuffer = Index3D(local_ncn-2, 1, k,local_ncn,local_ncm,local_ncl);
					it = ibox_k1[IndexDisplacedAndBuffer].begin();
					if (nc_k1[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k1[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num);
						}
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 18, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1)
			{
				SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 18, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 18, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 18, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//East-South-Down:17
		if (PN != 1 && PM != 1 && PL != 1)
		{
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				NumOfParticleSend = 0;
				IndexDisplacedAndBuffer = Index3D(local_ncn-2,1,1,local_ncn,local_ncm,local_ncl);
				vector<size_t>::iterator it = ibox_k1[IndexDisplacedAndBuffer].begin();
				if (nc_k1[IndexDisplacedAndBuffer] != 0)
				{
					for(;it!=ibox_k1[IndexDisplacedAndBuffer].end();++it)
					{
						size_t num;
						num=*it;
						DisplacedAndBufferParticle.push_back(num);
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 17, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != PL - 1)
			{
				SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2] + 1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 17, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 17, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
					//			BufferParticle.clear();
				}
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != PL - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 117, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
		}
		//Mid-North-Up:16
		if (PM != 1 && PL != 1)
		{
			if (PartitionIndex[1] != PM - 1 && PartitionIndex[2] != PL -1)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int i = local_ncn_interior[0]; i<local_ncn_interior[1]; i++)
				{

					IndexDisplacedAndBuffer = Index3D(i,local_ncm-2 ,local_ncl-2,local_ncn,local_ncm,local_ncl);
					it = ibox_k1[IndexDisplacedAndBuffer].begin();
					if (nc_k1[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k1[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num);
						}
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 16, MPI_COMM_WORLD);
			}
			if (PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 16, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[1] != PM - 1 && PartitionIndex[2] != PL -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 16, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 16, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//Mid-North-Mid:15
		if (PM !=1)
		{
			if (PartitionIndex[1] != PM - 1)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int i = local_ncn_interior[0]; i<local_ncn_interior[1];i++)
					for (int k = local_ncl_interior[0]; k<local_ncl_interior[1]; k++)
					{

						IndexDisplacedAndBuffer = Index3D(i,local_ncm-2,k,local_ncn,local_ncm,local_ncl);
						it = ibox_k1[IndexDisplacedAndBuffer].begin();
						if (nc_k1[IndexDisplacedAndBuffer] != 0)
						{
							for(;it!=ibox_k1[IndexDisplacedAndBuffer].end();++it)
							{
								size_t num;
								num=*it;
								DisplacedAndBufferParticle.push_back(num);
							}
						}
					}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 15, MPI_COMM_WORLD);
			}
			if (PartitionIndex[1] != 0)
			{
				SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 15, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[1] != PM - 1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 15, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[1] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 15, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//Mid-North-Down:14
		if (PM != 1 && PL != 1)
		{
			if (PartitionIndex[1] != PM - 1 && PartitionIndex[2] != 0)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int i = local_ncn_interior[0]; i<local_ncn_interior[1]; i++)
				{

					IndexDisplacedAndBuffer = Index3D(i,local_ncm-2,1,local_ncn,local_ncm,local_ncl);
					it = ibox_k1[IndexDisplacedAndBuffer].begin();
					if (nc_k1[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k1[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num);
						}
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 14, MPI_COMM_WORLD);
			}
			if (PartitionIndex[1] != 0 && PartitionIndex[2] != PL - 1)
			{
				SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2]+1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 14, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[1] != PM - 1 && PartitionIndex[2] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 14, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[1] != 0 && PartitionIndex[2] != PL - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 14, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//Mid-Mid-UP:13
		if (PL != 1)
		{
			if (PartitionIndex[2] != PL -1)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int i = local_ncn_interior[0]; i<local_ncn_interior[1];i++)
					for (int j = local_ncm_interior[0];j<local_ncm_interior[1]; j++)
					{

						IndexDisplacedAndBuffer = Index3D(i,j,local_ncl-2,local_ncn,local_ncm,local_ncl);
						it = ibox_k1[IndexDisplacedAndBuffer].begin();
						if (nc_k1[IndexDisplacedAndBuffer] != 0)
						{
							for(;it!=ibox_k1[IndexDisplacedAndBuffer].end();++it)
							{
								size_t num;
								num=*it;
								DisplacedAndBufferParticle.push_back(num);
							}
						}
					}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 13, MPI_COMM_WORLD);
			}
			if (PartitionIndex[2] != 0)
			{
				SourceID = Index3D(PartitionIndex[0],PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 13, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[2] != PL -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 13, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
					//			BufferParticle.clear();
				}
			}
			if (PartitionIndex[2] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0],PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 13, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//Mid-Mid-Down:12
		if (PL != 1)
		{
			if (PartitionIndex[2] != 0)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int i = local_ncn_interior[0]; i<local_ncn_interior[1];i++)
					for (int j = local_ncm_interior[0];j<local_ncm_interior[1]; j++)
					{

						IndexDisplacedAndBuffer = Index3D(i,j,1,local_ncn,local_ncm,local_ncl);
						it = ibox_k1[IndexDisplacedAndBuffer].begin();
						if (nc_k1[IndexDisplacedAndBuffer] != 0)
						{
							for(;it!=ibox_k1[IndexDisplacedAndBuffer].end();++it)
							{
								size_t num;
								num=*it;
								DisplacedAndBufferParticle.push_back(num);
							}
						}
					}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 12, MPI_COMM_WORLD);
			}
			if (PartitionIndex[2] != PL - 1)
			{
				SourceID = Index3D(PartitionIndex[0],PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 12, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[2] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 12, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[2] != PL - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0],PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 12, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//Mid-South-Up:11
		if (PM != 1 && PL != 1)
		{
			if (PartitionIndex[1] != 0 && PartitionIndex[2] != PL -1)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int i = local_ncn_interior[0]; i<local_ncn_interior[1]; i++)
				{

					IndexDisplacedAndBuffer = Index3D(i,1,local_ncl-2,local_ncn,local_ncm,local_ncl);
					it = ibox_k1[IndexDisplacedAndBuffer].begin();
					if (nc_k1[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k1[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num);
						}
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2]+1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 11, MPI_COMM_WORLD);
			}
			if (PartitionIndex[1] != PM - 1 && PartitionIndex[2] != 0)
			{
				SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 11, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[1] != 0 && PartitionIndex[2] != PL -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 11, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[1] != PM - 1 && PartitionIndex[2] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 11, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//Mid-South-Mid:10
		if (PM != 1)
		{
			if (PartitionIndex[1] != 0)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int i = local_ncn_interior[0]; i<local_ncn_interior[1];i++)
					for (int k = local_ncl_interior[0]; k<local_ncl_interior[1]; k++)
					{

						IndexDisplacedAndBuffer = Index3D(i,1,k,local_ncn,local_ncm,local_ncl);
						it = ibox_k1[IndexDisplacedAndBuffer].begin();
						if (nc_k1[IndexDisplacedAndBuffer] != 0)
						{
							for(;it!=ibox_k1[IndexDisplacedAndBuffer].end();++it)
							{
								size_t num;
								num=*it;
								DisplacedAndBufferParticle.push_back(num);
							}
						}
					}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 10, MPI_COMM_WORLD);
			}
			if (PartitionIndex[1] != PM - 1)
			{
				SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 10, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[1] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 10, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[1] != PM - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 10, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//Mid-South-Down:9
		if (PM != 1 && PL != 1)
		{
			if (PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int i = local_ncn_interior[0]; i<local_ncn_interior[1]; i++)
				{

					IndexDisplacedAndBuffer = Index3D(i,1,1,local_ncn,local_ncm,local_ncl);
					it = ibox_k1[IndexDisplacedAndBuffer].begin();
					if (nc_k1[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k1[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num);
						}
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 19, MPI_COMM_WORLD);
			}
			if (PartitionIndex[1] != PM - 1 && PartitionIndex[2] != PL - 1)
			{
				SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 9, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 9, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[1] != PM - 1 && PartitionIndex[2] != PL - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 9, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//West-North-Up:8
		if (PN != 1 && PM != 1 && PL != 1)
		{
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != PL -1)
			{
				NumOfParticleSend = 0;
				IndexDisplacedAndBuffer = Index3D(1,local_ncm-1,local_ncl-2,local_ncn,local_ncm,local_ncl);
				vector<size_t>::iterator it = ibox_k1[IndexDisplacedAndBuffer].begin();
				if (nc_k1[IndexDisplacedAndBuffer] != 0)
				{
					for(;it!=ibox_k1[IndexDisplacedAndBuffer].end();++it)
					{
						size_t num;
						num=*it;
						DisplacedAndBufferParticle.push_back(num);
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 8, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 8, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != PL -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 8, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 8, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//West-North-Mid:7
		if (PN != 1 && PM != 1)
		{
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int k = local_ncl_interior[0]; k<local_ncl_interior[1]; k++)
				{

					IndexDisplacedAndBuffer = Index3D(1,local_ncm-2,k,local_ncn,local_ncm,local_ncl);
					it = ibox_k1[IndexDisplacedAndBuffer].begin();
					if (nc_k1[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k1[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num);
						}
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 7, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != 0)
			{
				SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 7, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 7, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 7, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//West-North-Down:6
		if (PN != 1 && PM != 1 && PL != 1)
		{
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != 0)
			{
				NumOfParticleSend = 0;
				IndexDisplacedAndBuffer = Index3D(1,local_ncm-2,1,local_ncn,local_ncm,local_ncl);
				vector<size_t>::iterator it = ibox_k1[IndexDisplacedAndBuffer].begin();
				if (nc_k1[IndexDisplacedAndBuffer] != 0)
				{
					for(;it!=ibox_k1[IndexDisplacedAndBuffer].end();++it)
					{
						size_t num;
						num=*it;
						DisplacedAndBufferParticle.push_back(num);
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 6, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != 0 && PartitionIndex[2] != PL - 1)
			{
				SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2] + 1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 6, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 6, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != 0 && PartitionIndex[2] != PL - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 6, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//West-Mid-UP:5
		if (PN != 1 && PL != 1)
		{
			if (PartitionIndex[0] != 0 && PartitionIndex[2] != PL -1)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int j = local_ncm_interior[0];j<local_ncm_interior[1]; j++)
				{

					IndexDisplacedAndBuffer = Index3D(1,j,local_ncl-2,local_ncn,local_ncm,local_ncl);
					it = ibox_k1[IndexDisplacedAndBuffer].begin();
					if (nc_k1[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k1[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num);
						}
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 5, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[2] != 0)
			{
				SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 5, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != 0 && PartitionIndex[2] != PL -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 5, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[2] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 5, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//West-Mid-Mid:4
		if (PN != 1)
		{
			if (PartitionIndex[0] != 0)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int j = local_ncm_interior[0]; j<local_ncm_interior[1];j++)
					for (int k = local_ncl_interior[0]; k<local_ncl_interior[1]; k++)
					{

						IndexDisplacedAndBuffer = Index3D(1,j,k,local_ncn,local_ncm,local_ncl);
						it = ibox_k1[IndexDisplacedAndBuffer].begin();
						if (nc_k1[IndexDisplacedAndBuffer] != 0)
						{
							for(;it!=ibox_k1[IndexDisplacedAndBuffer].end();++it)
							{
								size_t num;
								num=*it;
								DisplacedAndBufferParticle.push_back(num);
							}
						}
					}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2],PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 4, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != PN - 1)
			{
				SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2],PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 4, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2],PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 4, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != PN - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2],PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 4, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//West-Mid-Down:3
		if (PN!= 1 && PL!=1)
		{
			if (PartitionIndex[0] != 0 && PartitionIndex[2] != 0)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int j = local_ncm_interior[0];j<local_ncm_interior[1]; j++)
				{

					IndexDisplacedAndBuffer = Index3D(1,j,1,local_ncn,local_ncm,local_ncl);
					it = ibox_k1[IndexDisplacedAndBuffer].begin();
					if (nc_k1[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k1[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num);
						}
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 3, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[2] != PL - 1)
			{
				SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 3, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != 0 && PartitionIndex[2] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 3, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[2] != PL - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 3, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//West-South-Up:2
		if (PN != 1 && PM != 1 && PL != 1)
		{
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0 && PartitionIndex[2] != PL -1)
			{
				NumOfParticleSend = 0;
				IndexDisplacedAndBuffer = Index3D(1,1,local_ncl-2,local_ncn,local_ncm,local_ncl);
				vector<size_t>::iterator it = ibox_k1[IndexDisplacedAndBuffer].begin();
				if (nc_k1[IndexDisplacedAndBuffer] != 0)
				{
					for(;it!=ibox_k1[IndexDisplacedAndBuffer].end();++it)
					{
						size_t num;
						num=*it;
						DisplacedAndBufferParticle.push_back(num);
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2]+1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 2, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != 0)
			{
				SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 2, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0 && PartitionIndex[2] != PL -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 2, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 2, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//West-South-Mid:1
		if (PN != 1 && PM != 1)
		{
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int k = local_ncl_interior[0]; k<local_ncl_interior[1]; k++)
				{

					IndexDisplacedAndBuffer = Index3D(1,1,k,local_ncn,local_ncm,local_ncl);
					it = ibox_k1[IndexDisplacedAndBuffer].begin();
					if (nc_k1[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k1[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num);
						}
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 1, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != PM - 1)
			{
				SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 1, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 1, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != PM - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 1, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//West-South-Down:0
		if (PN != 1 && PM != 1 && PL != 1)
		{
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				NumOfParticleSend = 0;
				IndexDisplacedAndBuffer = Index3D(1,1,1,local_ncn,local_ncm,local_ncl);
				vector<size_t>::iterator it = ibox_k1[IndexDisplacedAndBuffer].begin();
				if (nc_k1[IndexDisplacedAndBuffer] != 0)
				{
					for(;it!=ibox_k1[IndexDisplacedAndBuffer].end();++it)
					{
						size_t num;
						num=*it;
						DisplacedAndBufferParticle.push_back(num);
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 0, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != PL - 1)
			{
				SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 0, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 0, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != PL - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 0, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		if (TotalNumberOfParticleRecv !=0)
		{
			if (TotalNumberOfParticleRecv != DisplacedAndBufferParticlepm.size())
				assert(0);
			size_t nptemp = nbb1 + TotalNumberOfParticleRecv + np - nbb;
			size_t nbbtemp = nbb1 + TotalNumberOfParticleRecv;
			size_t np_max_temp = np_max;
			if (nptemp >= np_max_temp)
			{
				delete[] p;
				delete[] xpdot;
				delete[] ypdot;
				delete[] zpdot;
				delete[] updot;
				delete[] vpdot;
				delete[] wpdot;
				delete[] rhodot;
				delete[] ax;
				delete[] ay;
				delete[] az;
				delete[] ar;
				delete[] ux;
				delete[] vx;
				delete[] wx;
				delete[] aTE;
				delete[] xcor;
				delete[] ycor;
				delete[] zcor;
				delete [] cs;
				for (size_t i = 0; i < np_max; ++i)
				{
					for (int j = 0; j < 4; ++j)
						delete [] A_MLS[i][j];
					delete [] A_MLS[i];
				}
				delete [] A_MLS;
				delete[] rho_sum_MLS;
				delete[] sum_wab;
				delete[] rho_sum;
				delete[] beta0_MLS;
				delete[] beta1x_MLS;
				delete[] beta1y_MLS;
				delete[] beta1z_MLS;
				np_max = nptemp + nptemp / 5; //new np_max
				p = new double[np_max];
				xpdot = new double[np_max];
				ypdot = new double[np_max];
				zpdot = new double[np_max];
				updot = new double[np_max];
				vpdot = new double[np_max];
				wpdot = new double[np_max];
				TEdot = new double[np_max];
				rhodot = new double[np_max];
				ax = new double[np_max];
				ay = new double[np_max];
				az = new double[np_max];
				ar = new double[np_max];
				ux = new double[np_max];
				vx = new double[np_max];
				wx = new double[np_max];
				aTE = new double[np_max];
				xcor = new double[np_max];
				ycor = new double[np_max];
				zcor = new double[np_max];
				cs = new double[np_max];
				A_MLS = new double**[np_max];
				for (size_t i = 0; i < np_max; ++i)
				{
					A_MLS[i] = new double*[4];
					for (int j = 0; j < 4; ++j)
						A_MLS[i][j] = new double[4];
				}
				rho_sum_MLS = new double[np_max];
				sum_wab = new double[np_max];
				rho_sum = new double[np_max];
				beta0_MLS = new double[np_max];
				beta1x_MLS = new double[np_max];
				beta1y_MLS = new double[np_max];
				beta1z_MLS = new double[np_max];

			}
			//			double * xpt = new double[nptemp];
			//			double * ypt = new double[nptemp];
			//			double * zpt = new double[nptemp];
			//			double * upt = new double[nptemp];
			//			double * vpt = new double[nptemp];
			//			double * wpt = new double[nptemp];
			//			double * pt = new double[nptemp];
			//			double * rhot = new double[nptemp];
			//			double * pmt = new double[nptemp];
			//			double * xpot=0;
			//			double * ypot=0;
			//			double * zpot=0;
			//			double * upot=0;
			//			double * vpot=0;
			//			double * wpot=0;
			//			double * pot=0;
			//			double * rhoot=0;
			//			if (i_SR_step == 1)
			//			{
			//				xpot = new double[nptemp];
			//				ypot = new double[nptemp];
			//				zpot = new double[nptemp];
			//				upot = new double[nptemp];
			//				vpot = new double[nptemp];
			//				wpot = new double[nptemp];
			//				pot = new double[nptemp];
			//				rhoot = new double[nptemp];
			//			}
			size_t count = 0;
			for (size_t i = 0; i<nbb1; i++)
			{
				xpdot[i] = xp[i];
				ypdot[i] = yp[i];
				zpdot[i] = zp[i];
				updot[i] = up[i];
				vpdot[i] = vp[i];
				wpdot[i] = wp[i];
				TEdot[i] = TE[i];
				rhodot[i] = rho[i];
				xcor[i] = pm[i];
				if (i_SR_step == 1)
				{
					ax[i] = xpo[i];
					ay[i] = ypo[i];
					az[i] = zpo[i];
					ux[i] = upo[i];
					vx[i] = vpo[i];
					wx[i] = wpo[i];
					aTE[i] = TEo[i];
					ar[i] = rhoo[i];
				}
				++count;
			}
			for (size_t i = 0; i < TotalNumberOfParticleRecv; i++)
			{
				xpdot[count] = (i_SR_step == 1)?DisplacedAndBufferParticlexp[i*2]:DisplacedAndBufferParticlexp[i];
				ypdot[count] = (i_SR_step == 1)?DisplacedAndBufferParticleyp[i*2]:DisplacedAndBufferParticleyp[i];
				zpdot[count] = (i_SR_step == 1)?DisplacedAndBufferParticlezp[i*2]:DisplacedAndBufferParticlezp[i];
				updot[count] = (i_SR_step == 1)?DisplacedAndBufferParticleup[i*2]:DisplacedAndBufferParticleup[i];
				vpdot[count] = (i_SR_step == 1)?DisplacedAndBufferParticlevp[i*2]:DisplacedAndBufferParticlevp[i];
				wpdot[count] = (i_SR_step == 1)?DisplacedAndBufferParticlewp[i*2]:DisplacedAndBufferParticlewp[i];
				TEdot[count] = (i_SR_step == 1)?DisplacedAndBufferParticleTE[i*2]:DisplacedAndBufferParticleTE[i];
				rhodot[count] = (i_SR_step == 1)?DisplacedAndBufferParticlerho[i*2]:DisplacedAndBufferParticlerho[i];
				xcor[count] = DisplacedAndBufferParticlepm[i];
				if (i_SR_step == 1)
				{
					ax[count] = DisplacedAndBufferParticlexp[i*2+1];
					ay[count] = DisplacedAndBufferParticleyp[i*2+1];
					az[count] = DisplacedAndBufferParticlezp[i*2+1];
					ux[count] = DisplacedAndBufferParticleup[i*2+1];
					vx[count] = DisplacedAndBufferParticlevp[i*2+1];
					wx[count] = DisplacedAndBufferParticlewp[i*2+1];
					aTE[count] = DisplacedAndBufferParticleTE[i*2+1];
					ar[count] = DisplacedAndBufferParticlerho[i*2+1];
				}
				++count;
			}
			for (size_t i = nbb; i<np; i++)
			{
				xpdot[count] = xp[i];
				ypdot[count] = yp[i];
				zpdot[count] = zp[i];
				updot[count] = up[i];
				vpdot[count] = vp[i];
				wpdot[count] = wp[i];
				TEdot[count] = TE[i];
				rhodot[count] = rho[i];
				xcor[count] = pm[i];
				if (i_SR_step == 1)
				{
					ax[count] = xpo[i];
					ay[count] = ypo[i];
					az[count] = zpo[i];
					ux[count] = upo[i];
					vx[count] = vpo[i];
					wx[count] = wpo[i];
					aTE[count] = TEo[i];
					ar[count] = rhoo[i];
				}
				++count;
			}
			nb = nbbtemp + nb - nbb;
			nbb = nbbtemp;
			if (nptemp != count) assert(0);
			np = nptemp;
			if (nptemp >= np_max_temp)
			{
				delete [] xp;
				delete [] yp;
				delete [] zp;
				delete [] up;
				delete [] vp;
				delete [] wp;
				delete [] TE;
				delete [] rho;
				delete [] pm;
				delete [] xpo;
				delete [] ypo;
				delete [] zpo;
				delete [] upo;
				delete [] vpo;
				delete [] wpo;
				delete [] TEo;
				delete [] rhoo;
				delete [] bflag;
				xp = new double[np_max];
				yp = new double[np_max];
				zp = new double[np_max];
				up = new double[np_max];
				vp = new double[np_max];
				wp = new double[np_max];
				TE = new double[np_max];
				rho = new double[np_max];
				pm = new double[np_max];
				xpo = new double[np_max];
				ypo = new double[np_max];
				zpo = new double[np_max];
				upo = new double[np_max];
				vpo = new double[np_max];
				wpo = new double[np_max];
				TEo = new double[np_max];
				rhoo = new double[np_max];
				bflag = new int[np_max];
			}
			for (size_t i = 0; i< np; i++ )
			{
				xp[i] = xpdot[i];
				yp[i] = ypdot[i];
				zp[i] = zpdot[i];
				up[i] = updot[i];
				vp[i] = vpdot[i];
				wp[i] = wpdot[i];
				TE[i] = TEdot[i];
				rho[i] = rhodot[i];
				pm[i] = xcor[i];
				if (i_SR_step == 1)
				{
					xpo[i] = ax[i];
					ypo[i] = ay[i];
					zpo[i] = az[i];
					upo[i] = ux[i];
					vpo[i] = vx[i];
					wpo[i] = wx[i];
					TEo[i] = aTE[i];
					rhoo[i] = ar[i];
				}
				else
				{
					xpo[i] = 0;
					ypo[i] = 0;
					zpo[i] = 0;
					upo[i] = 0;
					vpo[i] = 0;
					wpo[i] = 0;
					TEo[i] = 0;
					rhoo[i] = 0;
				}
			}
			//			delete [] xpt;
			//			delete [] ypt;
			//			delete [] zpt;
			//			delete [] upt;
			//			delete [] vpt;
			//			delete [] wpt;
			//			delete [] pt;
			//			delete [] rhot;
			//			delete [] pmt;
			//			if (i_SR_step == 1)
			//			{
			//				delete [] xpot;
			//				delete [] ypot;
			//				delete [] zpot;
			//				delete [] upot;
			//				delete [] vpot;
			//				delete [] wpot;
			//				delete [] pot;
			//				delete [] rhoot;
			//			}
			DisplacedAndBufferParticlexp.clear();
			DisplacedAndBufferParticleyp.clear();
			DisplacedAndBufferParticlezp.clear();
			DisplacedAndBufferParticleup.clear();
			DisplacedAndBufferParticlevp.clear();
			DisplacedAndBufferParticlewp.clear();
			DisplacedAndBufferParticleTE.clear();
			DisplacedAndBufferParticlerho.clear();
			DisplacedAndBufferParticlepm.clear();
			DisplacedAndBufferParticle.clear();
			divide(0,nbb,1);
			divide(nb,np,2);
			IndexDisplacedAndBuffer = 0;
			TotalNumberOfParticleSend = TotalNumberOfParticleRecv = 0;
		}
		//26 Neighbours
		//East-North-Up:25
		if (PN != 1 && PM != 1 && PL != 1)
		{
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != PM -1 && PartitionIndex[2] != PL -1)
			{
				NumOfParticleSend = 0;
				IndexDisplacedAndBuffer = Index3D(local_ncn-2,local_ncm-2,local_ncl-2,local_ncn,local_ncm,local_ncl);
				vector<size_t>::iterator it = ibox_k2[IndexDisplacedAndBuffer].begin();
				if (nc_k2[IndexDisplacedAndBuffer] != 0)
				{
					for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
					{
						size_t num;
						num=*it;
						DisplacedAndBufferParticle.push_back(num);
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 25, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 25, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != PM -1 && PartitionIndex[2] != PL -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 25, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 25, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//East-North-Mid:24
		if (PN != 1 && PM != 1)
		{
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != PM -1)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int k = local_ncl_interior[0]; k<local_ncl_interior[1]; k++)
				{

					IndexDisplacedAndBuffer = Index3D(local_ncn-2,local_ncm-2,k,local_ncn,local_ncm,local_ncl);
					it = ibox_k2[IndexDisplacedAndBuffer].begin();
					if (nc_k2[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num);
						}
					}
				}

				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 24, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0)
			{
				SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 24, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != PM -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 24, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 24, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//East-North-Down:23
		if (PN != 1 && PM != 1 && PL != 1)
		{
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != PM -1 && PartitionIndex[2] != 0)
			{
				NumOfParticleSend = 0;
				IndexDisplacedAndBuffer = Index3D(local_ncn-2,local_ncm-2,1,local_ncn,local_ncm,local_ncl);
				vector<size_t>::iterator it = ibox_k2[IndexDisplacedAndBuffer].begin();
				if (nc_k2[IndexDisplacedAndBuffer] != 0)
				{
					for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
					{
						size_t num;
						num=*it;
						DisplacedAndBufferParticle.push_back(num);
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 23, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0 && PartitionIndex[2] != PL - 1)
			{
				SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2] + 1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 23, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != PM -1 && PartitionIndex[2] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 23, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0 && PartitionIndex[2] != PL - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 23, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//East-Mid-UP:22
		if (PN != 1 && PL != 1)
		{
			if (PartitionIndex[0] != PN -1 && PartitionIndex[2] != PL -1)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int j = local_ncm_interior[0];j<local_ncm_interior[1]; j++)
				{

					IndexDisplacedAndBuffer = Index3D(local_ncn-2,j,local_ncl-2,local_ncn,local_ncm,local_ncl);
					it = ibox_k2[IndexDisplacedAndBuffer].begin();
					if (nc_k2[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num);
						}
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 22, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[2] != 0)
			{
				SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 22, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != PN -1 && PartitionIndex[2] != PL -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 22, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
					//			BufferParticle.clear();
				}
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[2] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 22, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//East-Mid-Mid:21
		if (PN != 1)
		{
			if (PartitionIndex[0] != PN -1)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int j = local_ncm_interior[0]; j<local_ncm_interior[1];j++)
					for (int k = local_ncl_interior[0]; k<local_ncl_interior[1]; k++)
					{
						IndexDisplacedAndBuffer = Index3D(local_ncn-2,j,k,local_ncn,local_ncm,local_ncl);
						it = ibox_k2[IndexDisplacedAndBuffer].begin();
						if (nc_k2[IndexDisplacedAndBuffer] != 0)
						{
							for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
							{
								size_t num;
								num=*it;
								DisplacedAndBufferParticle.push_back(num);
							}
						}
					}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2],PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 21, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != 0)
			{
				SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2],PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 21, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != PN -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2],PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 21, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2],PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 21, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//East-Mid-Down:20
		if (PN!= 1 && PL!=1)
		{
			if (PartitionIndex[0] != PN -1 && PartitionIndex[2] != 0)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int j = local_ncm_interior[0];j<local_ncm_interior[1]; j++)
				{

					IndexDisplacedAndBuffer = Index3D(local_ncn-2,j,1,local_ncn,local_ncm,local_ncl);
					it = ibox_k2[IndexDisplacedAndBuffer].begin();
					if (nc_k2[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num);
						}
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 20, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[2] != PL - 1)
			{
				SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 20, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != PN -1 && PartitionIndex[2] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 20, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[2] != PL - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 20, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//East-South-Up:19
		if (PN != 1 && PM != 1 && PL != 1)
		{
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != 0 && PartitionIndex[2] != PL -1)
			{
				NumOfParticleSend = 0;
				IndexDisplacedAndBuffer = Index3D(local_ncn-2,1,local_ncl-2,local_ncn,local_ncm,local_ncl);
				vector<size_t>::iterator it = ibox_k2[IndexDisplacedAndBuffer].begin();
				if (nc_k2[IndexDisplacedAndBuffer] != 0)
				{
					for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
					{
						size_t num;
						num=*it;
						DisplacedAndBufferParticle.push_back(num);
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2]+1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 19, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != 0)
			{
				SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 19, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != 0 && PartitionIndex[2] != PL -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 19, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
					//			BufferParticle.clear();
				}
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 19, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//East-South-Mid:18
		if (PN != 1 && PM != 1)
		{
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != 0)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int k = local_ncl_interior[0]; k<local_ncl_interior[1]; k++)
				{

					IndexDisplacedAndBuffer = Index3D(local_ncn-2, 1, k,local_ncn,local_ncm,local_ncl);
					it = ibox_k2[IndexDisplacedAndBuffer].begin();
					if (nc_k2[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num);
						}
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 18, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1)
			{
				SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 18, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 18, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 18, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//East-South-Down:17
		if (PN != 1 && PM != 1 && PL != 1)
		{
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				NumOfParticleSend = 0;
				IndexDisplacedAndBuffer = Index3D(local_ncn-2,1,1,local_ncn,local_ncm,local_ncl);
				vector<size_t>::iterator it = ibox_k2[IndexDisplacedAndBuffer].begin();
				if (nc_k2[IndexDisplacedAndBuffer] != 0)
				{
					for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
					{
						size_t num;
						num=*it;
						DisplacedAndBufferParticle.push_back(num);
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 17, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != PL - 1)
			{
				SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2] + 1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 17, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != PN -1 && PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 17, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
					//			BufferParticle.clear();
				}
			}
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != PL - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 117, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
		}
		//Mid-North-Up:16
		if (PM != 1 && PL != 1)
		{
			if (PartitionIndex[1] != PM - 1 && PartitionIndex[2] != PL -1)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int i = local_ncn_interior[0]; i<local_ncn_interior[1]; i++)
				{

					IndexDisplacedAndBuffer = Index3D(i,local_ncm-2 ,local_ncl-2,local_ncn,local_ncm,local_ncl);
					it = ibox_k2[IndexDisplacedAndBuffer].begin();
					if (nc_k2[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num);
						}
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 16, MPI_COMM_WORLD);
			}
			if (PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 16, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[1] != PM - 1 && PartitionIndex[2] != PL -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 16, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 16, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//Mid-North-Mid:15
		if (PM !=1)
		{
			if (PartitionIndex[1] != PM - 1)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int i = local_ncn_interior[0]; i<local_ncn_interior[1];i++)
					for (int k = local_ncl_interior[0]; k<local_ncl_interior[1]; k++)
					{

						IndexDisplacedAndBuffer = Index3D(i,local_ncm-2,k,local_ncn,local_ncm,local_ncl);
						it = ibox_k2[IndexDisplacedAndBuffer].begin();
						if (nc_k2[IndexDisplacedAndBuffer] != 0)
						{
							for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
							{
								size_t num;
								num=*it;
								DisplacedAndBufferParticle.push_back(num);
							}
						}
					}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 15, MPI_COMM_WORLD);
			}
			if (PartitionIndex[1] != 0)
			{
				SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 15, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[1] != PM - 1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 15, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[1] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 15, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//Mid-North-Down:14
		if (PM != 1 && PL != 1)
		{
			if (PartitionIndex[1] != PM - 1 && PartitionIndex[2] != 0)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int i = local_ncn_interior[0]; i<local_ncn_interior[1]; i++)
				{

					IndexDisplacedAndBuffer = Index3D(i,local_ncm-2,1,local_ncn,local_ncm,local_ncl);
					it = ibox_k2[IndexDisplacedAndBuffer].begin();
					if (nc_k2[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num);
						}
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 14, MPI_COMM_WORLD);
			}
			if (PartitionIndex[1] != 0 && PartitionIndex[2] != PL - 1)
			{
				SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2]+1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 14, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[1] != PM - 1 && PartitionIndex[2] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 14, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[1] != 0 && PartitionIndex[2] != PL - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 14, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//Mid-Mid-UP:13
		if (PL != 1)
		{
			if (PartitionIndex[2] != PL -1)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int i = local_ncn_interior[0]; i<local_ncn_interior[1];i++)
					for (int j = local_ncm_interior[0];j<local_ncm_interior[1]; j++)
					{

						IndexDisplacedAndBuffer = Index3D(i,j,local_ncl-2,local_ncn,local_ncm,local_ncl);
						it = ibox_k2[IndexDisplacedAndBuffer].begin();
						if (nc_k2[IndexDisplacedAndBuffer] != 0)
						{
							for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
							{
								size_t num;
								num=*it;
								DisplacedAndBufferParticle.push_back(num);
							}
						}
					}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 13, MPI_COMM_WORLD);
			}
			if (PartitionIndex[2] != 0)
			{
				SourceID = Index3D(PartitionIndex[0],PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 13, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[2] != PL -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 13, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
					//			BufferParticle.clear();
				}
			}
			if (PartitionIndex[2] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0],PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 13, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//Mid-Mid-Down:12
		if (PL != 1)
		{
			if (PartitionIndex[2] != 0)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int i = local_ncn_interior[0]; i<local_ncn_interior[1];i++)
					for (int j = local_ncm_interior[0];j<local_ncm_interior[1]; j++)
					{

						IndexDisplacedAndBuffer = Index3D(i,j,1,local_ncn,local_ncm,local_ncl);
						it = ibox_k2[IndexDisplacedAndBuffer].begin();
						if (nc_k2[IndexDisplacedAndBuffer] != 0)
						{
							for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
							{
								size_t num;
								num=*it;
								DisplacedAndBufferParticle.push_back(num);
							}
						}
					}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 12, MPI_COMM_WORLD);
			}
			if (PartitionIndex[2] != PL - 1)
			{
				SourceID = Index3D(PartitionIndex[0],PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 12, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[2] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 12, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[2] != PL - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0],PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 12, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//Mid-South-Up:11
		if (PM != 1 && PL != 1)
		{
			if (PartitionIndex[1] != 0 && PartitionIndex[2] != PL -1)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int i = local_ncn_interior[0]; i<local_ncn_interior[1]; i++)
				{

					IndexDisplacedAndBuffer = Index3D(i,1,local_ncl-2,local_ncn,local_ncm,local_ncl);
					it = ibox_k2[IndexDisplacedAndBuffer].begin();
					if (nc_k2[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num);
						}
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2]+1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 11, MPI_COMM_WORLD);
			}
			if (PartitionIndex[1] != PM - 1 && PartitionIndex[2] != 0)
			{
				SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 11, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[1] != 0 && PartitionIndex[2] != PL -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 11, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[1] != PM - 1 && PartitionIndex[2] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 11, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//Mid-South-Mid:10
		if (PM != 1)
		{
			if (PartitionIndex[1] != 0)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int i = local_ncn_interior[0]; i<local_ncn_interior[1];i++)
					for (int k = local_ncl_interior[0]; k<local_ncl_interior[1]; k++)
					{

						IndexDisplacedAndBuffer = Index3D(i,1,k,local_ncn,local_ncm,local_ncl);
						it = ibox_k2[IndexDisplacedAndBuffer].begin();
						if (nc_k2[IndexDisplacedAndBuffer] != 0)
						{
							for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
							{
								size_t num;
								num=*it;
								DisplacedAndBufferParticle.push_back(num);
							}
						}
					}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 10, MPI_COMM_WORLD);
			}
			if (PartitionIndex[1] != PM - 1)
			{
				SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 10, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[1] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 10, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[1] != PM - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 10, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//Mid-South-Down:9
		if (PM != 1 && PL != 1)
		{
			if (PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int i = local_ncn_interior[0]; i<local_ncn_interior[1]; i++)
				{

					IndexDisplacedAndBuffer = Index3D(i,1,1,local_ncn,local_ncm,local_ncl);
					it = ibox_k2[IndexDisplacedAndBuffer].begin();
					if (nc_k2[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num);
						}
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 19, MPI_COMM_WORLD);
			}
			if (PartitionIndex[1] != PM - 1 && PartitionIndex[2] != PL - 1)
			{
				SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 9, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0],PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 9, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[1] != PM - 1 && PartitionIndex[2] != PL - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0],PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 9, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//West-North-Up:8
		if (PN != 1 && PM != 1 && PL != 1)
		{
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != PL -1)
			{
				NumOfParticleSend = 0;
				IndexDisplacedAndBuffer = Index3D(1,local_ncm-1,local_ncl-2,local_ncn,local_ncm,local_ncl);
				vector<size_t>::iterator it = ibox_k2[IndexDisplacedAndBuffer].begin();
				if (nc_k2[IndexDisplacedAndBuffer] != 0)
				{
					for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
					{
						size_t num;
						num=*it;
						DisplacedAndBufferParticle.push_back(num);
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 8, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 8, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != PL -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 8, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 8, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//West-North-Mid:7
		if (PN != 1 && PM != 1)
		{
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int k = local_ncl_interior[0]; k<local_ncl_interior[1]; k++)
				{

					IndexDisplacedAndBuffer = Index3D(1,local_ncm-2,k,local_ncn,local_ncm,local_ncl);
					it = ibox_k2[IndexDisplacedAndBuffer].begin();
					if (nc_k2[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num);
						}
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 7, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != 0)
			{
				SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 7, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 7, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 7, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//West-North-Down:6
		if (PN != 1 && PM != 1 && PL != 1)
		{
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != 0)
			{
				NumOfParticleSend = 0;
				IndexDisplacedAndBuffer = Index3D(1,local_ncm-2,1,local_ncn,local_ncm,local_ncl);
				vector<size_t>::iterator it = ibox_k2[IndexDisplacedAndBuffer].begin();
				if (nc_k2[IndexDisplacedAndBuffer] != 0)
				{
					for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
					{
						size_t num;
						num=*it;
						DisplacedAndBufferParticle.push_back(num);
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 6, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != 0 && PartitionIndex[2] != PL - 1)
			{
				SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2] + 1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 6, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 6, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != 0 && PartitionIndex[2] != PL - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]-1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 6, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//West-Mid-UP:5
		if (PN != 1 && PL != 1)
		{
			if (PartitionIndex[0] != 0 && PartitionIndex[2] != PL -1)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int j = local_ncm_interior[0];j<local_ncm_interior[1]; j++)
				{

					IndexDisplacedAndBuffer = Index3D(1,j,local_ncl-2,local_ncn,local_ncm,local_ncl);
					it = ibox_k2[IndexDisplacedAndBuffer].begin();
					if (nc_k2[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num);
						}
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 5, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[2] != 0)
			{
				SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 5, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != 0 && PartitionIndex[2] != PL -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 5, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[2] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 5, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//West-Mid-Mid:4
		if (PN != 1)
		{
			if (PartitionIndex[0] != 0)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int j = local_ncm_interior[0]; j<local_ncm_interior[1];j++)
					for (int k = local_ncl_interior[0]; k<local_ncl_interior[1]; k++)
					{

						IndexDisplacedAndBuffer = Index3D(1,j,k,local_ncn,local_ncm,local_ncl);
						it = ibox_k2[IndexDisplacedAndBuffer].begin();
						if (nc_k2[IndexDisplacedAndBuffer] != 0)
						{
							for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
							{
								size_t num;
								num=*it;
								DisplacedAndBufferParticle.push_back(num);
							}
						}
					}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2],PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 4, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != PN - 1)
			{
				SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2],PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 4, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2],PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 4, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != PN - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2],PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 4, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//West-Mid-Down:3
		if (PN!= 1 && PL!=1)
		{
			if (PartitionIndex[0] != 0 && PartitionIndex[2] != 0)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int j = local_ncm_interior[0];j<local_ncm_interior[1]; j++)
				{

					IndexDisplacedAndBuffer = Index3D(1,j,1,local_ncn,local_ncm,local_ncl);
					it = ibox_k2[IndexDisplacedAndBuffer].begin();
					if (nc_k2[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num);
						}
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 3, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[2] != PL - 1)
			{
				SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 3, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != 0 && PartitionIndex[2] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1],PartitionIndex[2]-1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 3, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[2] != PL - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1],PartitionIndex[2]+1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 3, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//West-South-Up:2
		if (PN != 1 && PM != 1 && PL != 1)
		{
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0 && PartitionIndex[2] != PL -1)
			{
				NumOfParticleSend = 0;
				IndexDisplacedAndBuffer = Index3D(1,1,local_ncl-2,local_ncn,local_ncm,local_ncl);
				vector<size_t>::iterator it = ibox_k2[IndexDisplacedAndBuffer].begin();
				if (nc_k2[IndexDisplacedAndBuffer] != 0)
				{
					for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
					{
						size_t num;
						num=*it;
						DisplacedAndBufferParticle.push_back(num);
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2]+1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 2, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != 0)
			{
				SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 2, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0 && PartitionIndex[2] != PL -1)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 2, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != 0)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 2, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//West-South-Mid:1
		if (PN != 1 && PM != 1)
		{
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0)
			{
				NumOfParticleSend = 0;
				vector<size_t>::iterator it;
				for (int k = local_ncl_interior[0]; k<local_ncl_interior[1]; k++)
				{

					IndexDisplacedAndBuffer = Index3D(1,1,k,local_ncn,local_ncm,local_ncl);
					it = ibox_k2[IndexDisplacedAndBuffer].begin();
					if (nc_k2[IndexDisplacedAndBuffer] != 0)
					{
						for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
						{
							size_t num;
							num=*it;
							DisplacedAndBufferParticle.push_back(num);
						}
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 1, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != PM - 1)
			{
				SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 1, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2],PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 1, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != PM - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2],PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 1, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//West-South-Down:0
		if (PN != 1 && PM != 1 && PL != 1)
		{
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				NumOfParticleSend = 0;
				IndexDisplacedAndBuffer = Index3D(1,1,1,local_ncn,local_ncm,local_ncl);
				vector<size_t>::iterator it = ibox_k2[IndexDisplacedAndBuffer].begin();
				if (nc_k2[IndexDisplacedAndBuffer] != 0)
				{
					for(;it!=ibox_k2[IndexDisplacedAndBuffer].end();++it)
					{
						size_t num;
						num=*it;
						DisplacedAndBufferParticle.push_back(num);
					}
				}
				NumOfParticleSend = DisplacedAndBufferParticle.size();
				DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
				MPI_Send(&NumOfParticleSend, 1, MPI_UNSIGNED_LONG, DestinationID, 0, MPI_COMM_WORLD);
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != PL - 1)
			{
				SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
				MPI_Recv(&NumOfParticleRecv, 1, MPI_UNSIGNED_LONG, SourceID, 0, MPI_COMM_WORLD, &status);
				int temp = NumOfParticleRecv;
				TotalNumberOfParticleRecv += NumOfParticleRecv;
				if (temp != 0)
				{
					size_t size = (i_SR_step == 1)?(temp * 17):(temp * 9);
					RecvBuffer =  (double *) malloc (size * sizeof(double));
					for (size_t i = 0; i< size;i++)
						RecvBuffer[i] = 0;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (PartitionIndex[0] != 0 && PartitionIndex[1] != 0 && PartitionIndex[2] != 0)
			{
				int temp = NumOfParticleSend;
				if (temp != 0)
				{
					SendBuffer =  (double *) malloc (((i_SR_step == 1)?(temp * 17):(temp * 9) )* sizeof(double));
					int index = 0;
					vector<size_t>::iterator it;
					for(it = DisplacedAndBufferParticle.begin();it!=DisplacedAndBufferParticle.end();++it)
					{
						size_t i = *it;
						SendBuffer[index] = xp[i];
						SendBuffer[index+1] = yp[i];
						SendBuffer[index+2] = zp[i];
						SendBuffer[index+3] = up[i];
						SendBuffer[index+4] = vp[i];
						SendBuffer[index+5] = wp[i];
						SendBuffer[index+6] = TE[i];
						SendBuffer[index+7] = rho[i];
						SendBuffer[index+8] = pm[i];
						if  (i_SR_step == 1)
						{
							SendBuffer[index+9] = xpo[i];
							SendBuffer[index+10] = ypo[i];
							SendBuffer[index+11] = zpo[i];
							SendBuffer[index+12] = upo[i];
							SendBuffer[index+13] = vpo[i];
							SendBuffer[index+14] = wpo[i];
							SendBuffer[index+15] =TEo[i];
							SendBuffer[index+16] = rhoo[i];
						}
						index += (i_SR_step == 1)?17:9;
					}
					DestinationID = Index3D(PartitionIndex[0]-1,PartitionIndex[1]-1,PartitionIndex[2]-1,PN,PM,PL);
					MPI_Send(SendBuffer, (i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, DestinationID, 0, MPI_COMM_WORLD);
					DisplacedAndBufferParticle.clear();
					free (SendBuffer);
					SendBuffer = NULL;
				}
			}
			if (PartitionIndex[0] != PN - 1 && PartitionIndex[1] != PM - 1 && PartitionIndex[2] != PL - 1)
			{
				size_t temp = NumOfParticleRecv;
				if (temp != 0)
				{
					SourceID = Index3D(PartitionIndex[0]+1,PartitionIndex[1]+1,PartitionIndex[2]+1,PN,PM,PL);
					MPI_Recv(RecvBuffer,(i_SR_step == 1)?(temp * 17):(temp * 9), MPI_DOUBLE, SourceID, 0, MPI_COMM_WORLD, &status);
					size_t count = 0;
					for (size_t i=0; i<NumOfParticleRecv;i++,count+=(i_SR_step == 1)?17:9)
					{
						DisplacedAndBufferParticlexp.push_back(RecvBuffer[count]);
						DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+1]);
						DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+2]);
						DisplacedAndBufferParticleup.push_back(RecvBuffer[count+3]);
						DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+4]);
						DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+5]);
						DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+6]);
						DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+7]);
						DisplacedAndBufferParticlepm.push_back(RecvBuffer[count+8]);
						if  (i_SR_step == 1)
						{
							DisplacedAndBufferParticlexp.push_back(RecvBuffer[count+9]);
							DisplacedAndBufferParticleyp.push_back(RecvBuffer[count+10]);
							DisplacedAndBufferParticlezp.push_back(RecvBuffer[count+11]);
							DisplacedAndBufferParticleup.push_back(RecvBuffer[count+12]);
							DisplacedAndBufferParticlevp.push_back(RecvBuffer[count+13]);
							DisplacedAndBufferParticlewp.push_back(RecvBuffer[count+14]);
							DisplacedAndBufferParticleTE.push_back(RecvBuffer[count+15]);
							DisplacedAndBufferParticlerho.push_back(RecvBuffer[count+16]);
						};
					}
					free (RecvBuffer);
					RecvBuffer = NULL;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		if (TotalNumberOfParticleRecv !=0)
		{
			if (TotalNumberOfParticleRecv != DisplacedAndBufferParticlepm.size())
				assert(0);
			size_t nptemp = nbb + TotalNumberOfParticleRecv + np - nb;
			size_t nbtemp = nbb + TotalNumberOfParticleRecv;
			size_t np_max_temp = np_max;
			if (nptemp >= np_max_temp)
			{
				delete[] p;
				delete[] xpdot;
				delete[] ypdot;
				delete[] zpdot;
				delete[] updot;
				delete[] vpdot;
				delete[] wpdot;
				delete[] rhodot;
				delete[] ax;
				delete[] ay;
				delete[] az;
				delete[] ar;
				delete[] ux;
				delete[] vx;
				delete[] wx;
				delete[] aTE;
				delete[] xcor;
				delete[] ycor;
				delete[] zcor;
				delete [] cs;
				for (size_t i = 0; i < np_max; ++i)
				{
					for (int j = 0; j < 4; ++j)
						delete [] A_MLS[i][j];
					delete [] A_MLS[i];
				}
				delete [] A_MLS;
				delete[] rho_sum_MLS;
				delete[] sum_wab;
				delete[] rho_sum;
				delete[] beta0_MLS;
				delete[] beta1x_MLS;
				delete[] beta1y_MLS;
				delete[] beta1z_MLS;
				np_max = nptemp + nptemp / 5; //new np_max
				p = new double [np_max];
				xpdot = new double[np_max];
				ypdot = new double[np_max];
				zpdot = new double[np_max];
				updot = new double[np_max];
				vpdot = new double[np_max];
				wpdot = new double[np_max];
				TEdot = new double[np_max];
				rhodot = new double[np_max];
				ax = new double[np_max];
				ay = new double[np_max];
				az = new double[np_max];
				ar = new double[np_max];
				ux = new double[np_max];
				vx = new double[np_max];
				wx = new double[np_max];
				aTE = new double[np_max];
				xcor = new double[np_max];
				ycor = new double[np_max];
				zcor = new double[np_max];
				cs = new double[np_max];
				A_MLS = new double**[np_max];
				for (size_t i = 0; i < np_max; ++i)
				{
					A_MLS[i] = new double*[4];
					for (int j = 0; j < 4; ++j)
						A_MLS[i][j] = new double[4];
				}
				rho_sum_MLS = new double[np_max];
				sum_wab = new double[np_max];
				rho_sum = new double[np_max];
				beta0_MLS = new double[np_max];
				beta1x_MLS = new double[np_max];
				beta1y_MLS = new double[np_max];
				beta1z_MLS = new double[np_max];

			}
			//			double * xpt = new double[nptemp];
			//			double * ypt = new double[nptemp];
			//			double * zpt = new double[nptemp];
			//			double * upt = new double[nptemp];
			//			double * vpt = new double[nptemp];
			//			double * wpt = new double[nptemp];
			//			double * pt = new double[nptemp];
			//			double * rhot = new double[nptemp];
			//			double * pmt = new double[nptemp];
			//			double * xpot=0;
			//			double * ypot=0;
			//			double * zpot=0;
			//			double * upot=0;
			//			double * vpot=0;
			//			double * wpot=0;
			//			double * pot=0;
			//			double * rhoot=0;
			//			if (i_SR_step == 1)
			//			{
			//				xpot = new double[nptemp];
			//				ypot = new double[nptemp];
			//				zpot = new double[nptemp];
			//				upot = new double[nptemp];
			//				vpot = new double[nptemp];
			//				wpot = new double[nptemp];
			//				pot = new double[nptemp];
			//				rhoot = new double[nptemp];
			//			}
			size_t count = 0;
			for (size_t i = 0; i<nbb; i++)
			{
				xpdot[i] = xp[i];
				ypdot[i] = yp[i];
				zpdot[i] = zp[i];
				updot[i] = up[i];
				vpdot[i] = vp[i];
				wpdot[i] = wp[i];
				TEdot[i] = TE[i];
				rhodot[i] = rho[i];
				xcor[i] = pm[i];
				if (i_SR_step == 1)
				{
					ax[i] = xpo[i];
					ay[i] = ypo[i];
					az[i] = zpo[i];
					ux[i] = upo[i];
					vx[i] = vpo[i];
					wx[i] = wpo[i];
					aTE[i] = TEo[i];
					ar[i] = rhoo[i];
				}
				++count;
			}
			for (size_t i = 0; i < TotalNumberOfParticleRecv; i++)
			{
				xpdot[count] = (i_SR_step == 1)?DisplacedAndBufferParticlexp[i*2]:DisplacedAndBufferParticlexp[i];
				ypdot[count] = (i_SR_step == 1)?DisplacedAndBufferParticleyp[i*2]:DisplacedAndBufferParticleyp[i];
				zpdot[count] = (i_SR_step == 1)?DisplacedAndBufferParticlezp[i*2]:DisplacedAndBufferParticlezp[i];
				updot[count] = (i_SR_step == 1)?DisplacedAndBufferParticleup[i*2]:DisplacedAndBufferParticleup[i];
				vpdot[count] = (i_SR_step == 1)?DisplacedAndBufferParticlevp[i*2]:DisplacedAndBufferParticlevp[i];
				wpdot[count] = (i_SR_step == 1)?DisplacedAndBufferParticlewp[i*2]:DisplacedAndBufferParticlewp[i];
				TEdot[count] = (i_SR_step == 1)?DisplacedAndBufferParticleTE[i*2]:DisplacedAndBufferParticleTE[i];
				rhodot[count] = (i_SR_step == 1)?DisplacedAndBufferParticlerho[i*2]:DisplacedAndBufferParticlerho[i];
				xcor[count] = DisplacedAndBufferParticlepm[i];
				if (i_SR_step == 1)
				{
					ax[count] = DisplacedAndBufferParticlexp[i*2+1];
					ay[count] = DisplacedAndBufferParticleyp[i*2+1];
					az[count] = DisplacedAndBufferParticlezp[i*2+1];
					ux[count] = DisplacedAndBufferParticleup[i*2+1];
					vx[count] = DisplacedAndBufferParticlevp[i*2+1];
					wx[count] = DisplacedAndBufferParticlewp[i*2+1];
					aTE[count] = DisplacedAndBufferParticleTE[i*2+1];
					ar[count] = DisplacedAndBufferParticlerho[i*2+1];
				}
				++count;
			}
			for (size_t i = nb; i<np; i++)
			{
				xpdot[count] = xp[i];
				ypdot[count] = yp[i];
				zpdot[count] = zp[i];
				updot[count] = up[i];
				vpdot[count] = vp[i];
				wpdot[count] = wp[i];
				TEdot[count] = TE[i];
				rhodot[count] = rho[i];
				xcor[count] = pm[i];
				if (i_SR_step == 1)
				{
					ax[count] = xpo[i];
					ay[count] = ypo[i];
					az[count] = zpo[i];
					ux[count] = upo[i];
					vx[count] = vpo[i];
					wx[count] = wpo[i];
					aTE[count] = TEo[i];
					ar[count] = rhoo[i];
				}
				++count;
			}
			nb = nbtemp;
			if (nptemp != count) assert(0);
			np = nptemp;
			if (nptemp >= np_max_temp)
			{
				delete [] xp;
				delete [] yp;
				delete [] zp;
				delete [] up;
				delete [] vp;
				delete [] wp;
				delete [] TE;
				delete [] rho;
				delete [] pm;
				delete [] xpo;
				delete [] ypo;
				delete [] zpo;
				delete [] upo;
				delete [] vpo;
				delete [] wpo;
				delete [] TEo;
				delete [] rhoo;
				delete [] bflag;
				xp = new double[np_max];
				yp = new double[np_max];
				zp = new double[np_max];
				up = new double[np_max];
				vp = new double[np_max];
				wp = new double[np_max];
				TE = new double[np_max];
				rho = new double[np_max];
				pm = new double[np_max];
				xpo = new double[np_max];
				ypo = new double[np_max];
				zpo = new double[np_max];
				upo = new double[np_max];
				vpo = new double[np_max];
				wpo = new double[np_max];
				TEo = new double[np_max];
				rhoo = new double[np_max];
				bflag = new int[np_max];
			}
			for (size_t i = 0; i< np; i++ )
			{
				xp[i] = xpdot[i];
				yp[i] = ypdot[i];
				zp[i] = zpdot[i];
				up[i] = updot[i];
				vp[i] = vpdot[i];
				wp[i] = wpdot[i];
				TE[i] = TEdot[i];
				rho[i] = rhodot[i];
				pm[i] = xcor[i];
				if (i_SR_step == 1)
				{
					xpo[i] = ax[i];
					ypo[i] = ay[i];
					zpo[i] = az[i];
					upo[i] = ux[i];
					vpo[i] = vx[i];
					wpo[i] = wx[i];
					TEo[i] = aTE[i];
					rhoo[i] = ar[i];
				}
				else
				{
					xpo[i] = 0;
					ypo[i] = 0;
					zpo[i] = 0;
					upo[i] = 0;
					vpo[i] = 0;
					wpo[i] = 0;
					TEo[i] = 0;
					rhoo[i] = 0;
				}
			}
			//			delete [] xpt;
			//			delete [] ypt;
			//			delete [] zpt;
			//			delete [] upt;
			//			delete [] vpt;
			//			delete [] wpt;
			//			delete [] pt;
			//			delete [] rhot;
			//			delete [] pmt;
			//			if (i_SR_step == 1)
			//			{
			//				delete [] xpot;
			//				delete [] ypot;
			//				delete [] zpot;
			//				delete [] upot;
			//				delete [] vpot;
			//				delete [] wpot;
			//				delete [] pot;
			//				delete [] rhoot;
			//			}
			DisplacedAndBufferParticlexp.clear();
			DisplacedAndBufferParticleyp.clear();
			DisplacedAndBufferParticlezp.clear();
			DisplacedAndBufferParticleup.clear();
			DisplacedAndBufferParticlevp.clear();
			DisplacedAndBufferParticlewp.clear();
			DisplacedAndBufferParticleTE.clear();
			DisplacedAndBufferParticlerho.clear();
			DisplacedAndBufferParticlepm.clear();
			DisplacedAndBufferParticle.clear();
			divide(nbb,np,2);
		}
	}
}

void STORAGE::ParticleNumberSync()
{
	size_t * num_array;
	int gsize =0;
	int root = 0;
	MPI_Barrier(MPI_COMM_WORLD);
	if (id == root)
	{
		MPI_Comm_size(MPI_COMM_WORLD,&gsize);
		num_array = new size_t[gsize];
	}
	MPI_Gather(&np, 1,MPI_UNSIGNED_LONG, num_array, 1, MPI_UNSIGNED_LONG, root, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if (id == root)
	{
		std::vector<size_t> myvector(num_array,num_array+gsize);
		std::sort(myvector.begin(),myvector.end());
		np_max = myvector[gsize-1];
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (id == root)
	{
		MPI_Bcast(&np_max, 1, MPI_UNSIGNED_LONG, root, MPI_COMM_WORLD);
	}
	else
	{
		MPI_Bcast(&np_max, 1, MPI_UNSIGNED_LONG, root, MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

void STORAGE::StatesPrint(int step, double time, char* statesname, const char* outputname)
{
	//buildmesh
	int n0,n1,m0,m1,l0,l1;
	if ((localvlx[0] - vlx[0]< 0) || (localvlx[0] == vlx[0]))
		n0 = 0;
	else
		n0 = (int)((localvlx[0]-vlx[0])/dx) + 1;
	//	if (local_ncn_interior[0] != 0)
	//		n0--;
	n1 = (int)((min(vlx[1],localvlx[1])+(1e-8)*dx-vlx[0])/dx);
	if (local_ncn_interior[1] != local_ncn)
		n1++;
	if (localvly[0] - vly[0]< 0)
		m0 = 0;
	else
		m0 = (int)((localvly[0]-vly[0])/dy) + 1;
	m1 = (int)((min(vly[1],localvly[1])+(1e-8)*dy-vly[0])/dy);
	if (localvlz[0] - vlz[0]< 0)
		l0 = 0;
	else
		l0 = (int)((localvlz[0]-vlz[0])/dz) + 1;
	l1 = (int)((min(vlz[1],localvlz[1])+(1e-8)*dz-vlz[0])/dz);
	int dimension[3] = {n1-n0 +1, m1-m0+1,l1-l0+1};
	//	scalar = new double**[n1-n0+1];
	//	for (int i = 0; i<(n1-n0+1); i++)
	//	{
	//		scalar[i] = new double*[m1-m0+1];
	//		for (int j = 0;j<(m1-m0+1); j++)
	//			scalar[i][j] = new double[l1-l0+1];
	//	}
	for (int i = 0; i<(n1-n0+1); i++)
		for (int j = 0;j<(m1-m0+1); j++)
			for (int k = 0;k<(l1-l0+1);k++)
			{
				int index = Index3D(i,j,k,(n1-n0+1),(m1-m0+1),(l1-l0+1));
				if (statesname[0] == 'D' || statesname[0] == 'p')
					scalar[index] = 0;
				else if (statesname[0] == 'p' || statesname[0] == 'p')
					scalar[index] = 0;
				else
					scalar[index] = -1;
			}
	double origin[3];
	origin[0] = vlx[0] + n0 * dx;
	origin[1] = vly[0] + m0 * dy;
	origin[2] = vlz[0] + l0 * dz;
	//finish mesh building

	//locate position
	for (size_t i = nbb;i<np;i++)
	{
		if (xp[i] <= vlx[0]+0.0000001 || xp[i] >= vlx[1]-0.0000001 ||
				yp[i] <= vly[0]+0.0000001 || yp[i] >= vly[1]-0.0000001 ||
				zp[i] <= vlz[0]+0.0000001 || zp[i] >= vlz[1]-0.0000001)
			continue;
		double distx,disty,distz; //to the origin
		distx = xp[i] - origin[0];
		disty = yp[i] - origin[1];
		distz = zp[i] - origin[2];
		int icell,jcell,kcell;
		if (distx >= -dx && distx <0)
		{
			icell = -1;
		}
		else if (distx >=0 && distx < (n1-n0 + 1) * dx)
		{
			icell = (int)(distx / dx);
		}
		else
			continue;
		if (disty >= -dy && disty <0)
		{
			jcell = -1;
		}
		else if (disty >=0 && disty < (m1-m0 + 1) * dy)
		{
			jcell = (int)(disty / dy);
		}
		else
			continue;
		if (distz >= -dz && distz <0)
		{
			kcell = -1;
		}
		else if (distz >=0 && distz < (l1-l0 + 1) * dz)
		{
			kcell = (int)(distz / dz);
		}
		else
			continue;
		if (icell>=0 && icell <= (n1-n0-1) &&
				jcell>=0 && jcell <= (m1-m0-1) &&
				kcell>=0 && kcell <= (l1-l0-1))
		{
			double x0,x1,y0,y1,z0,z1;
			x0 = origin[0] + icell * dx;
			x1 = x0 + dx;
			y0 = origin[1] + jcell * dy;
			y1 = y0 + dy;
			z0 = origin[2] + kcell * dz;
			z1 = z0 + dz;
			double distancex[2],distancey[2],distancez[2]; //to faces
			distancex[0] = (xp[i] - x0)/dx;
			distancex[1] = (x1 - xp[i])/dx;
			distancey[0] = (yp[i] - y0)/dy;
			distancey[1] = (y1 - yp[i])/dy;
			distancez[0] = (zp[i] - z0)/dz;
			distancez[1] = (z1 - zp[i])/dz;
			//get volumn fraction
			double volumn[8]={0};
			//index
			//0:west-south-down; 1:east-south-down; 2:east-north-down; 3:west-north-down
			//4:west-south-up;   5:east-south-up;   6:east-north-up;   7:east-north-up
			volumn[0] = distancex[0] * distancey[0] * distancez[0];
			volumn[1] = distancex[1] * distancey[0] * distancez[0];
			volumn[2] = distancex[1] * distancey[1] * distancez[0];
			volumn[3] = distancex[0] * distancey[1] * distancez[0];
			volumn[4] = distancex[0] * distancey[0] * distancez[1];
			volumn[5] = distancex[1] * distancey[0] * distancez[1];
			volumn[6] = distancex[1] * distancey[1] * distancez[1];
			volumn[7] = distancex[0] * distancey[1] * distancez[1];
			double volumn_inverse[8];
			for (int count1 = 0; count1<8; count1++)
			{
				volumn_inverse[count1] = 1;
				for (int count2 = 0; count2<8; count2++)
				{
					if (count1 == count2) continue;
					else
						volumn_inverse[count1] = volumn_inverse[count1] * volumn[count2];
				}
			}
			double total=0;
			for (int count1 = 0; count1<8; count1++)
				total = total + volumn_inverse[count1];
			for (int count1 = 0; count1<8; count1++)
			{
				if (statesname[0] == 'D' || statesname[0] == 'p')
					volumn_inverse[count1] = rho[i] * volumn_inverse[count1] / total;
				else if (statesname[0] == 'p' || statesname[0] == 'p')
					volumn_inverse[count1] = p[i] * volumn_inverse[count1] / total;
				else
					volumn_inverse[count1] = 2 * volumn_inverse[count1] / total;
			}

			//deposit value
			int index;
			index = Index3D(icell,jcell,kcell,(n1-n0+1),(m1-m0+1),(l1-l0+1));
			scalar[index] += volumn_inverse[0];
			index = Index3D(icell+1,jcell,kcell,(n1-n0+1),(m1-m0+1),(l1-l0+1));
			scalar[index] += volumn_inverse[1];
			index = Index3D(icell+1,jcell+1,kcell,(n1-n0+1),(m1-m0+1),(l1-l0+1));
			scalar[index] += volumn_inverse[2];
			index = Index3D(icell,jcell+1,kcell,(n1-n0+1),(m1-m0+1),(l1-l0+1));
			scalar[index] += volumn_inverse[3];
			index = Index3D(icell,jcell,kcell+1,(n1-n0+1),(m1-m0+1),(l1-l0+1));
			scalar[index] += volumn_inverse[4];
			index = Index3D(icell+1,jcell,kcell+1,(n1-n0+1),(m1-m0+1),(l1-l0+1));
			scalar[index] += volumn_inverse[5];
			index = Index3D(icell+1,jcell+1,kcell+1,(n1-n0+1),(m1-m0+1),(l1-l0+1));
			scalar[index] += volumn_inverse[6];
			index = Index3D(icell,jcell+1,kcell+1,(n1-n0+1),(m1-m0+1),(l1-l0+1));
			scalar[index] += volumn_inverse[7];
		}
		//		else if(icell == -1 && jcell == -1 && kcell == -1)
		//		{
		//			double x0,x1,y0,y1,z0,z1;
		//			x0 = origin[0] - dx;
		//			x1 = x0 + dx;
		//			y0 = origin[1] - dy;
		//			y1 = y0 + dy;
		//			z0 = origin[2] - dz;
		//			z1 = z0 + dz;
		//			double distancex[2],distancey[2],distancez[2]; //to faces
		//			distancex[0] = (xp[i] - x0)/dx;
		//			distancex[1] = (x1 - xp[i])/dx;
		//			distancey[0] = (yp[i] - y0)/dy;
		//			distancey[1] = (y1 - yp[i])/dy;
		//			distancez[0] = (zp[i] - z0)/dz;
		//			distancez[1] = (z1 - zp[i])/dz;
		//			//get volumn fraction
		//			double volumn[8]={0};
		//			//index
		//			//0:west-south-down; 1:east-south-down; 2:east-north-down; 3:west-north-down
		//			//4:west-south-up;   5:east-south-up;   6:east-north-up;   7:east-north-up
		//			volumn[0] = distancex[0] * distancey[0] * distancez[0];
		//			volumn[1] = distancex[1] * distancey[0] * distancez[0];
		//			volumn[2] = distancex[1] * distancey[1] * distancez[0];
		//			volumn[3] = distancex[0] * distancey[1] * distancez[0];
		//			volumn[4] = distancex[0] * distancey[0] * distancez[1];
		//			volumn[5] = distancex[1] * distancey[0] * distancez[1];
		//			volumn[6] = distancex[1] * distancey[1] * distancez[1];
		//			volumn[7] = distancex[0] * distancey[1] * distancez[1];
		//			double volumn_inverse[8];
		//			for (int i = 0; i<8; i++)
		//			{
		//				volumn_inverse[i] = 1;
		//				for (int j = 0; j<8; j++)
		//				{
		//					if (i == j) continue;
		//					else
		//						volumn_inverse[i] = volumn_inverse[i] * volumn[j];
		//				}
		//			}
		//			double total=0;
		//			for (int i = 0; i<8; i++)
		//				total = total + volumn_inverse[i];
		//			for (int i = 0; i<8; i++)
		//				volumn_inverse[i] = 2 * volumn_inverse[i] / total;
		//
		//			//deposit value
		//			int index;
		//			index = Index3D(icell+1,jcell+1,kcell+1,(n1-n0+1),(m1-m0+1),(l1-l0+1));
		//			scalar[index] += volumn_inverse[6];
		//		}
		//		else if (icell==-1 && jcell==-1 &&
		//				kcell>=0 && kcell <= (l1-l0-1))
		//		{
		//			double x0,x1,y0,y1,z0,z1;
		//			x0 = origin[0] - dx;
		//			x1 = x0 + dx;
		//			y0 = origin[1] - dy;
		//			y1 = y0 + dy;
		//			z0 = origin[2] + kcell * dz;
		//			z1 = z0 + dz;
		//			double distancex[2],distancey[2],distancez[2]; //to faces
		//			distancex[0] = (xp[i] - x0)/dx;
		//			distancex[1] = (x1 - xp[i])/dx;
		//			distancey[0] = (yp[i] - y0)/dy;
		//			distancey[1] = (y1 - yp[i])/dy;
		//			distancez[0] = (zp[i] - z0)/dz;
		//			distancez[1] = (z1 - zp[i])/dz;
		//			//get volumn fraction
		//			double volumn[8]={0};
		//			//index
		//			//0:west-south-down; 1:east-south-down; 2:east-north-down; 3:west-north-down
		//			//4:west-south-up;   5:east-south-up;   6:east-north-up;   7:east-north-up
		//			volumn[0] = distancex[0] * distancey[0] * distancez[0];
		//			volumn[1] = distancex[1] * distancey[0] * distancez[0];
		//			volumn[2] = distancex[1] * distancey[1] * distancez[0];
		//			volumn[3] = distancex[0] * distancey[1] * distancez[0];
		//			volumn[4] = distancex[0] * distancey[0] * distancez[1];
		//			volumn[5] = distancex[1] * distancey[0] * distancez[1];
		//			volumn[6] = distancex[1] * distancey[1] * distancez[1];
		//			volumn[7] = distancex[0] * distancey[1] * distancez[1];
		//			double volumn_inverse[8];
		//			for (int i = 0; i<8; i++)
		//			{
		//				volumn_inverse[i] = 1;
		//				for (int j = 0; j<8; j++)
		//				{
		//					if (i == j) continue;
		//					else
		//						volumn_inverse[i] = volumn_inverse[i] * volumn[j];
		//				}
		//			}
		//			double total=0;
		//			for (int i = 0; i<8; i++)
		//				total = total + volumn_inverse[i];
		//			for (int i = 0; i<8; i++)
		//				volumn_inverse[i] = 2 * volumn_inverse[i] / total;
		//
		//			//deposit value
		//			int index;
		//			index = Index3D(icell+1,jcell+1,kcell,(n1-n0+1),(m1-m0+1),(l1-l0+1));
		//			scalar[index] += volumn_inverse[2];
		//			index = Index3D(icell+1,jcell+1,kcell+1,(n1-n0+1),(m1-m0+1),(l1-l0+1));
		//			scalar[index] += volumn_inverse[6];
		//		}
		else if (icell == (n1-n0) &&
				jcell>=0 && jcell <= (m1-m0-1) &&
				kcell>=0 && kcell <= (l1-l0-1))
		{
			double x0,x1,y0,y1,z0,z1;
			x0 = origin[0] + icell * dx;
			x1 = x0 + dx;
			y0 = origin[1] + jcell * dy;
			y1 = y0 + dy;
			z0 = origin[2] + kcell * dz;
			z1 = z0 + dz;
			double distancex[2],distancey[2],distancez[2]; //to faces
			distancex[0] = (xp[i] - x0)/dx;
			distancex[1] = (x1 - xp[i])/dx;
			distancey[0] = (yp[i] - y0)/dy;
			distancey[1] = (y1 - yp[i])/dy;
			distancez[0] = (zp[i] - z0)/dz;
			distancez[1] = (z1 - zp[i])/dz;
			//get volumn fraction
			double volumn[8]={0};
			//index
			//0:west-south-down; 1:east-south-down; 2:east-north-down; 3:west-north-down
			//4:west-south-up;   5:east-south-up;   6:east-north-up;   7:east-north-up
			volumn[0] = distancex[0] * distancey[0] * distancez[0];
			volumn[1] = distancex[1] * distancey[0] * distancez[0];
			volumn[2] = distancex[1] * distancey[1] * distancez[0];
			volumn[3] = distancex[0] * distancey[1] * distancez[0];
			volumn[4] = distancex[0] * distancey[0] * distancez[1];
			volumn[5] = distancex[1] * distancey[0] * distancez[1];
			volumn[6] = distancex[1] * distancey[1] * distancez[1];
			volumn[7] = distancex[0] * distancey[1] * distancez[1];
			double volumn_inverse[8];
			for (int count1 = 0; count1<8; count1++)
			{
				volumn_inverse[count1] = 1;
				for (int count2 = 0; count2<8; count2++)
				{
					if (count1 == count2) continue;
					else
						volumn_inverse[count1] = volumn_inverse[count1] * volumn[count2];
				}
			}
			double total=0;
			for (int count1 = 0; count1<8; count1++)
				total = total + volumn_inverse[count1];
			for (int count1 = 0; count1<8; count1++)
			{
				if (statesname[0] == 'D' || statesname[0] == 'p')
					volumn_inverse[count1] = rho[i] * volumn_inverse[count1] / total;
				else if (statesname[0] == 'p' || statesname[0] == 'p')
					volumn_inverse[count1] = p[i] * volumn_inverse[count1] / total;
				else
					volumn_inverse[count1] = 2 * volumn_inverse[count1] / total;
			}

			//deposit value
			int index;
			index = Index3D(icell,jcell,kcell,(n1-n0+1),(m1-m0+1),(l1-l0+1));
			scalar[index] += volumn_inverse[0];
			index = Index3D(icell,jcell+1,kcell,(n1-n0+1),(m1-m0+1),(l1-l0+1));
			scalar[index] += volumn_inverse[3];
			index = Index3D(icell,jcell,kcell+1,(n1-n0+1),(m1-m0+1),(l1-l0+1));
			scalar[index] += volumn_inverse[4];
			index = Index3D(icell,jcell+1,kcell+1,(n1-n0+1),(m1-m0+1),(l1-l0+1));
			scalar[index] += volumn_inverse[7];
		}
		else if (icell == -1 &&
				jcell>=0 && jcell <= (m1-m0-1) &&
				kcell>=0 && kcell <= (l1-l0-1))
		{
			double x0,x1,y0,y1,z0,z1;
			x0 = origin[0] + icell * dx;
			x1 = x0 + dx;
			y0 = origin[1] + jcell * dy;
			y1 = y0 + dy;
			z0 = origin[2] + kcell * dz;
			z1 = z0 + dz;
			double distancex[2],distancey[2],distancez[2]; //to faces
			distancex[0] = (xp[i] - x0)/dx;
			distancex[1] = (x1 - xp[i])/dx;
			distancey[0] = (yp[i] - y0)/dy;
			distancey[1] = (y1 - yp[i])/dy;
			distancez[0] = (zp[i] - z0)/dz;
			distancez[1] = (z1 - zp[i])/dz;
			//get volumn fraction
			double volumn[8]={0};
			//index
			//0:west-south-down; 1:east-south-down; 2:east-north-down; 3:west-north-down
			//4:west-south-up;   5:east-south-up;   6:east-north-up;   7:east-north-up
			volumn[0] = distancex[0] * distancey[0] * distancez[0];
			volumn[1] = distancex[1] * distancey[0] * distancez[0];
			volumn[2] = distancex[1] * distancey[1] * distancez[0];
			volumn[3] = distancex[0] * distancey[1] * distancez[0];
			volumn[4] = distancex[0] * distancey[0] * distancez[1];
			volumn[5] = distancex[1] * distancey[0] * distancez[1];
			volumn[6] = distancex[1] * distancey[1] * distancez[1];
			volumn[7] = distancex[0] * distancey[1] * distancez[1];
			double volumn_inverse[8];
			for (int count1 = 0; count1<8; count1++)
			{
				volumn_inverse[count1] = 1;
				for (int count2 = 0; count2<8; count2++)
				{
					if (count1 == count2) continue;
					else
						volumn_inverse[count1] = volumn_inverse[count1] * volumn[count2];
				}
			}
			double total=0;
			for (int count1 = 0; count1<8; count1++)
				total = total + volumn_inverse[count1];
			for (int count1 = 0; count1<8; count1++)
			{
				if (statesname[0] == 'D' || statesname[0] == 'p')
					volumn_inverse[count1] = rho[i] * volumn_inverse[count1] / total;
				else if (statesname[0] == 'p' || statesname[0] == 'p')
					volumn_inverse[count1] = p[i] * volumn_inverse[count1] / total;
				else
					volumn_inverse[count1] = 2 * volumn_inverse[count1] / total;
			}

			int index;
			index = Index3D(icell+1,jcell,kcell,(n1-n0+1),(m1-m0+1),(l1-l0+1));
			scalar[index] += volumn_inverse[1];
			index = Index3D(icell+1,jcell+1,kcell,(n1-n0+1),(m1-m0+1),(l1-l0+1));
			scalar[index] += volumn_inverse[2];
			index = Index3D(icell+1,jcell,kcell+1,(n1-n0+1),(m1-m0+1),(l1-l0+1));
			scalar[index] += volumn_inverse[5];
			index = Index3D(icell+1,jcell+1,kcell+1,(n1-n0+1),(m1-m0+1),(l1-l0+1));
			scalar[index] += volumn_inverse[6];
		}
		//		else if (icell>=0 && icell <= (n1-n0-1) &&
		//				jcell == (m1-m0) &&
		//				kcell>=0 && kcell <= (l1-l0-1))
		//		{
		//			double x0,x1,y0,y1,z0,z1;
		//			x0 = origin[0] + icell * dx;
		//			x1 = x0 + dx;
		//			y0 = origin[1] + jcell * dy;
		//			y1 = y0 + dy;
		//			z0 = origin[2] + kcell * dz;
		//			z1 = z0 + dz;
		//			double distancex[2],distancey[2],distancez[2]; //to faces
		//			distancex[0] = (xp[i] - x0)/dx;
		//			distancex[1] = (x1 - xp[i])/dx;
		//			distancey[0] = (yp[i] - y0)/dy;
		//			distancey[1] = (y1 - yp[i])/dy;
		//			distancez[0] = (zp[i] - z0)/dz;
		//			distancez[1] = (z1 - zp[i])/dz;
		//			//get volumn fraction
		//			double volumn[8]={0};
		//			//index
		//			//0:west-south-down; 1:east-south-down; 2:east-north-down; 3:west-north-down
		//			//4:west-south-up;   5:east-south-up;   6:east-north-up;   7:east-north-up
		//			volumn[0] = distancex[0] * distancey[0] * distancez[0];
		//			volumn[1] = distancex[1] * distancey[0] * distancez[0];
		//			volumn[2] = distancex[1] * distancey[1] * distancez[0];
		//			volumn[3] = distancex[0] * distancey[1] * distancez[0];
		//			volumn[4] = distancex[0] * distancey[0] * distancez[1];
		//			volumn[5] = distancex[1] * distancey[0] * distancez[1];
		//			volumn[6] = distancex[1] * distancey[1] * distancez[1];
		//			volumn[7] = distancex[0] * distancey[1] * distancez[1];
		//			double volumn_inverse[8];
		//			for (int i = 0; i<8; i++)
		//			{
		//				volumn_inverse[i] = 1;
		//				for (int j = 0; j<8; j++)
		//				{
		//					if (i == j) continue;
		//					else
		//						volumn_inverse[i] = volumn_inverse[i] * volumn[j];
		//				}
		//			}
		//			double total=0;
		//			for (int i = 0; i<8; i++)
		//				total = total + volumn_inverse[i];
		//			for (int i = 0; i<8; i++)
		//				volumn_inverse[i] = 2 * volumn_inverse[i] / total;
		//
		//			//deposit value
		//			int index;
		//			index = Index3D(icell,jcell,kcell,(n1-n0+1),(m1-m0+1),(l1-l0+1));
		//			scalar[index] += volumn_inverse[0];
		//			index = Index3D(icell+1,jcell,kcell,(n1-n0+1),(m1-m0+1),(l1-l0+1));
		//			scalar[index] += volumn_inverse[1];
		//			index = Index3D(icell,jcell,kcell+1,(n1-n0+1),(m1-m0+1),(l1-l0+1));
		//			scalar[index] += volumn_inverse[4];
		//			index = Index3D(icell+1,jcell,kcell+1,(n1-n0+1),(m1-m0+1),(l1-l0+1));
		//			scalar[index] += volumn_inverse[5];
		//		}
		//		else if (icell>=0 && icell <= (n1-n0-1) &&
		//				jcell == -1 &&
		//				kcell>=0 && kcell <= (l1-l0-1))
		//		{
		//			double x0,x1,y0,y1,z0,z1;
		//			x0 = origin[0] + icell * dx;
		//			x1 = x0 + dx;
		//			y0 = origin[1] + jcell * dy;
		//			y1 = y0 + dy;
		//			z0 = origin[2] + kcell * dz;
		//			z1 = z0 + dz;
		//			double distancex[2],distancey[2],distancez[2]; //to faces
		//			distancex[0] = (xp[i] - x0)/dx;
		//			distancex[1] = (x1 - xp[i])/dx;
		//			distancey[0] = (yp[i] - y0)/dy;
		//			distancey[1] = (y1 - yp[i])/dy;
		//			distancez[0] = (zp[i] - z0)/dz;
		//			distancez[1] = (z1 - zp[i])/dz;
		//			//get volumn fraction
		//			double volumn[8]={0};
		//			//index
		//			//0:west-south-down; 1:east-south-down; 2:east-north-down; 3:west-north-down
		//			//4:west-south-up;   5:east-south-up;   6:east-north-up;   7:east-north-up
		//			volumn[0] = distancex[0] * distancey[0] * distancez[0];
		//			volumn[1] = distancex[1] * distancey[0] * distancez[0];
		//			volumn[2] = distancex[1] * distancey[1] * distancez[0];
		//			volumn[3] = distancex[0] * distancey[1] * distancez[0];
		//			volumn[4] = distancex[0] * distancey[0] * distancez[1];
		//			volumn[5] = distancex[1] * distancey[0] * distancez[1];
		//			volumn[6] = distancex[1] * distancey[1] * distancez[1];
		//			volumn[7] = distancex[0] * distancey[1] * distancez[1];
		//			double volumn_inverse[8];
		//			for (int i = 0; i<8; i++)
		//			{
		//				volumn_inverse[i] = 1;
		//				for (int j = 0; j<8; j++)
		//				{
		//					if (i == j) continue;
		//					else
		//						volumn_inverse[i] = volumn_inverse[i] * volumn[j];
		//				}
		//			}
		//			double total=0;
		//			for (int i = 0; i<8; i++)
		//				total = total + volumn_inverse[i];
		//			for (int i = 0; i<8; i++)
		//				volumn_inverse[i] = 2 * volumn_inverse[i] / total;
		//
		//			//deposit value
		//			int index;
		//			index = Index3D(icell+1,jcell+1,kcell,(n1-n0+1),(m1-m0+1),(l1-l0+1));
		//			scalar[index] += volumn_inverse[2];
		//			index = Index3D(icell,jcell+1,kcell,(n1-n0+1),(m1-m0+1),(l1-l0+1));
		//			scalar[index] += volumn_inverse[3];
		//			index = Index3D(icell+1,jcell+1,kcell+1,(n1-n0+1),(m1-m0+1),(l1-l0+1));
		//			scalar[index] += volumn_inverse[6];
		//			index = Index3D(icell,jcell+1,kcell+1,(n1-n0+1),(m1-m0+1),(l1-l0+1));
		//			scalar[index] += volumn_inverse[7];
		//		}
		//		if (icell == -1 &&
		//				jcell == -1 &&
		//				kcell>=0 && kcell <= (l1-l0-1))
		//		{
		//			double x0,x1,y0,y1,z0,z1;
		//			x0 = origin[0] + icell * dx;
		//			x1 = x0 + dx;
		//			y0 = origin[1] + jcell * dy;
		//			y1 = y0 + dy;
		//			z0 = origin[2] + kcell * dz;
		//			z1 = z0 + dz;
		//			double distancex[2],distancey[2],distancez[2]; //to faces
		//			distancex[0] = (xp[i] - x0)/dx;
		//			distancex[1] = (x1 - xp[i])/dx;
		//			distancey[0] = (yp[i] - y0)/dy;
		//			distancey[1] = (y1 - yp[i])/dy;
		//			distancez[0] = (zp[i] - z0)/dz;
		//			distancez[1] = (z1 - zp[i])/dz;
		//			//get volumn fraction
		//			double volumn[8]={0};
		//			//index
		//			//0:west-south-down; 1:east-south-down; 2:east-north-down; 3:west-north-down
		//			//4:west-south-up;   5:east-south-up;   6:east-north-up;   7:east-north-up
		//			volumn[0] = distancex[0] * distancey[0] * distancez[0];
		//			volumn[1] = distancex[1] * distancey[0] * distancez[0];
		//			volumn[2] = distancex[1] * distancey[1] * distancez[0];
		//			volumn[3] = distancex[0] * distancey[1] * distancez[0];
		//			volumn[4] = distancex[0] * distancey[0] * distancez[1];
		//			volumn[5] = distancex[1] * distancey[0] * distancez[1];
		//			volumn[6] = distancex[1] * distancey[1] * distancez[1];
		//			volumn[7] = distancex[0] * distancey[1] * distancez[1];
		//			double volumn_inverse[8];
		//			for (int i = 0; i<8; i++)
		//			{
		//				volumn_inverse[i] = 1;
		//				for (int j = 0; j<8; j++)
		//				{
		//					if (i == j) continue;
		//					else
		//						volumn_inverse[i] = volumn_inverse[i] * volumn[j];
		//				}
		//			}
		//			double total=0;
		//			for (int i = 0; i<8; i++)
		//				total = total + volumn_inverse[i];
		//			for (int i = 0; i<8; i++)
		//				volumn_inverse[i] = 2 * volumn_inverse[i] / total;
		//
		//			//deposit value
		//			int index;
		//			index = Index3D(icell+1,jcell+1,kcell,(n1-n0+1),(m1-m0+1),(l1-l0+1));
		//			scalar[index] += volumn_inverse[2];
		//			index = Index3D(icell+1,jcell+1,kcell+1,(n1-n0+1),(m1-m0+1),(l1-l0+1));
		//			scalar[index] += volumn_inverse[6];
		//		}
		//		if (icell == -1 &&
		//				jcell == (m1-m0) &&
		//				kcell>=0 && kcell <= (l1-l0-1))
		//		{
		//			double x0,x1,y0,y1,z0,z1;
		//			x0 = origin[0] + icell * dx;
		//			x1 = x0 + dx;
		//			y0 = origin[1] + jcell * dy;
		//			y1 = y0 + dy;
		//			z0 = origin[2] + kcell * dz;
		//			z1 = z0 + dz;
		//			double distancex[2],distancey[2],distancez[2]; //to faces
		//			distancex[0] = (xp[i] - x0)/dx;
		//			distancex[1] = (x1 - xp[i])/dx;
		//			distancey[0] = (yp[i] - y0)/dy;
		//			distancey[1] = (y1 - yp[i])/dy;
		//			distancez[0] = (zp[i] - z0)/dz;
		//			distancez[1] = (z1 - zp[i])/dz;
		//			//get volumn fraction
		//			double volumn[8]={0};
		//			//index
		//			//0:west-south-down; 1:east-south-down; 2:east-north-down; 3:west-north-down
		//			//4:west-south-up;   5:east-south-up;   6:east-north-up;   7:east-north-up
		//			volumn[0] = distancex[0] * distancey[0] * distancez[0];
		//			volumn[1] = distancex[1] * distancey[0] * distancez[0];
		//			volumn[2] = distancex[1] * distancey[1] * distancez[0];
		//			volumn[3] = distancex[0] * distancey[1] * distancez[0];
		//			volumn[4] = distancex[0] * distancey[0] * distancez[1];
		//			volumn[5] = distancex[1] * distancey[0] * distancez[1];
		//			volumn[6] = distancex[1] * distancey[1] * distancez[1];
		//			volumn[7] = distancex[0] * distancey[1] * distancez[1];
		//			double volumn_inverse[8];
		//			for (int i = 0; i<8; i++)
		//			{
		//				volumn_inverse[i] = 1;
		//				for (int j = 0; j<8; j++)
		//				{
		//					if (i == j) continue;
		//					else
		//						volumn_inverse[i] = volumn_inverse[i] * volumn[j];
		//				}
		//			}
		//			double total=0;
		//			for (int i = 0; i<8; i++)
		//				total = total + volumn_inverse[i];
		//			for (int i = 0; i<8; i++)
		//				volumn_inverse[i] = 2 * volumn_inverse[i] / total;
		//
		//			//deposit value
		//			int index;
		//			index = Index3D(icell+1,jcell,kcell,(n1-n0+1),(m1-m0+1),(l1-l0+1));
		//			scalar[index] += volumn_inverse[1];
		//			index = Index3D(icell+1,jcell,kcell+1,(n1-n0+1),(m1-m0+1),(l1-l0+1));
		//			scalar[index] += volumn_inverse[5];
		//		}
		//		if (icell == n1-n0 &&
		//				jcell == -1 &&
		//				kcell>=0 && kcell <= (l1-l0-1))
		//		{
		//			double x0,x1,y0,y1,z0,z1;
		//			x0 = origin[0] + icell * dx;
		//			x1 = x0 + dx;
		//			y0 = origin[1] + jcell * dy;
		//			y1 = y0 + dy;
		//			z0 = origin[2] + kcell * dz;
		//			z1 = z0 + dz;
		//			double distancex[2],distancey[2],distancez[2]; //to faces
		//			distancex[0] = (xp[i] - x0)/dx;
		//			distancex[1] = (x1 - xp[i])/dx;
		//			distancey[0] = (yp[i] - y0)/dy;
		//			distancey[1] = (y1 - yp[i])/dy;
		//			distancez[0] = (zp[i] - z0)/dz;
		//			distancez[1] = (z1 - zp[i])/dz;
		//			//get volumn fraction
		//			double volumn[8]={0};
		//			//index
		//			//0:west-south-down; 1:east-south-down; 2:east-north-down; 3:west-north-down
		//			//4:west-south-up;   5:east-south-up;   6:east-north-up;   7:east-north-up
		//			volumn[0] = distancex[0] * distancey[0] * distancez[0];
		//			volumn[1] = distancex[1] * distancey[0] * distancez[0];
		//			volumn[2] = distancex[1] * distancey[1] * distancez[0];
		//			volumn[3] = distancex[0] * distancey[1] * distancez[0];
		//			volumn[4] = distancex[0] * distancey[0] * distancez[1];
		//			volumn[5] = distancex[1] * distancey[0] * distancez[1];
		//			volumn[6] = distancex[1] * distancey[1] * distancez[1];
		//			volumn[7] = distancex[0] * distancey[1] * distancez[1];
		//			double volumn_inverse[8];
		//			for (int i = 0; i<8; i++)
		//			{
		//				volumn_inverse[i] = 1;
		//				for (int j = 0; j<8; j++)
		//				{
		//					if (i == j) continue;
		//					else
		//						volumn_inverse[i] = volumn_inverse[i] * volumn[j];
		//				}
		//			}
		//			double total=0;
		//			for (int i = 0; i<8; i++)
		//				total = total + volumn_inverse[i];
		//			for (int i = 0; i<8; i++)
		//				volumn_inverse[i] = 2 * volumn_inverse[i] / total;
		//
		//			//deposit value
		//			int index;
		//			index = Index3D(icell,jcell+1,kcell,(n1-n0+1),(m1-m0+1),(l1-l0+1));
		//			scalar[index] += volumn_inverse[3];
		//			index = Index3D(icell,jcell+1,kcell+1,(n1-n0+1),(m1-m0+1),(l1-l0+1));
		//			scalar[index] += volumn_inverse[7];
		//		}
		//		if (icell == (n1-n0) &&
		//				jcell == (m1-m0) &&
		//				kcell>=0 && kcell <= (l1-l0-1))
		//		{
		//			double x0,x1,y0,y1,z0,z1;
		//			x0 = origin[0] + icell * dx;
		//			x1 = x0 + dx;
		//			y0 = origin[1] + jcell * dy;
		//			y1 = y0 + dy;
		//			z0 = origin[2] + kcell * dz;
		//			z1 = z0 + dz;
		//			double distancex[2],distancey[2],distancez[2]; //to faces
		//			distancex[0] = (xp[i] - x0)/dx;
		//			distancex[1] = (x1 - xp[i])/dx;
		//			distancey[0] = (yp[i] - y0)/dy;
		//			distancey[1] = (y1 - yp[i])/dy;
		//			distancez[0] = (zp[i] - z0)/dz;
		//			distancez[1] = (z1 - zp[i])/dz;
		//			//get volumn fraction
		//			double volumn[8]={0};
		//			//index
		//			//0:west-south-down; 1:east-south-down; 2:east-north-down; 3:west-north-down
		//			//4:west-south-up;   5:east-south-up;   6:east-north-up;   7:east-north-up
		//			volumn[0] = distancex[0] * distancey[0] * distancez[0];
		//			volumn[1] = distancex[1] * distancey[0] * distancez[0];
		//			volumn[2] = distancex[1] * distancey[1] * distancez[0];
		//			volumn[3] = distancex[0] * distancey[1] * distancez[0];
		//			volumn[4] = distancex[0] * distancey[0] * distancez[1];
		//			volumn[5] = distancex[1] * distancey[0] * distancez[1];
		//			volumn[6] = distancex[1] * distancey[1] * distancez[1];
		//			volumn[7] = distancex[0] * distancey[1] * distancez[1];
		//			double volumn_inverse[8];
		//			for (int i = 0; i<8; i++)
		//			{
		//				volumn_inverse[i] = 1;
		//				for (int j = 0; j<8; j++)
		//				{
		//					if (i == j) continue;
		//					else
		//						volumn_inverse[i] = volumn_inverse[i] * volumn[j];
		//				}
		//			}
		//			double total=0;
		//			for (int i = 0; i<8; i++)
		//				total = total + volumn_inverse[i];
		//			for (int i = 0; i<8; i++)
		//				volumn_inverse[i] = 2 * volumn_inverse[i] / total;
		//
		//			//deposit value
		//			int index;
		//			index = Index3D(icell,jcell,kcell,(n1-n0+1),(m1-m0+1),(l1-l0+1));
		//			scalar[index] += volumn_inverse[0];
		//			index = Index3D(icell,jcell,kcell+1,(n1-n0+1),(m1-m0+1),(l1-l0+1));
		//			scalar[index] += volumn_inverse[4];
		//		}
	}


	//start print out
	char filename[200];
	FILE *outfile;
	sprintf(filename,"vtkoutput_states_%s",right_flush(step,7));
	sprintf(filename,"%s_%s",filename,outputname);
	if (PN * PM * PL != 1)
		sprintf(filename,"%s-nd%s",filename,right_flush(id,4));
	sprintf(filename,"%s.vtk",filename);
	if (!exists(filename))
	{
		outfile = fopen(filename,"w");
		fprintf(outfile,"# vtk DataFile Version 3.0\n");
		fprintf(outfile,"The actual time is %.8f\n",time);
		fprintf(outfile,"ASCII\n");
		fprintf(outfile,"DATASET STRUCTURED_POINTS\n");
		fprintf(outfile,"DIMENSIONS %d %d %d\n",dimension[0],dimension[1],dimension[2]);
		fprintf(outfile,"SPACING %lf %lf %lf\n",dx,dy,dz);
		fprintf(outfile,"ORIGIN %lf %lf %lf\n", origin[0]+0.5*dx,origin[1]+0.5*dy,origin[2]+0.5*dz);
		fprintf(outfile,"POINT_DATA %d\n",(dimension[0])*(dimension[1])*(dimension[2]));
	}
	else
		outfile = fopen(filename,"a");
	if (statesname[0] == 'D' || statesname[0] == 'd')
		fprintf(outfile,"SCALARS Density double 1\n");
	else if (statesname[0] == 'P' || statesname[0] == 'p')
		fprintf(outfile,"SCALARS Pressure double 1\n");
	else
		fprintf(outfile,"SCALARS Interface double 1\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");

	for (int k = 0;k<(l1-l0+1);k++)
		for (int j = 0;j<(m1-m0+1); j++)
			for (int i = 0; i<(n1-n0+1); i++)
			{
				int index = Index3D(i,j,k,(n1-n0+1),(m1-m0+1),(l1-l0+1));
				fprintf(outfile,"%.16g\n",scalar[index]);
			}
	fclose(outfile);
}

int STORAGE::exists(const char *fname)
{
	FILE *file;
	if ((file = fopen(fname, "r")))
	{
		fclose(file);
		return 1;
	}
	return 0;
}

void STORAGE::particlegathering(int step)
{
	//		if (nprint != step)
	//			return;
	size_t num_array[4];
	int gsize =0;
	int root = 0;
	size_t number = np - nb;
	size_t sum=0;
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gather(&number, 1,MPI_UNSIGNED_LONG, num_array, 1, MPI_UNSIGNED_LONG, root, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(num_array, 4, MPI_UNSIGNED_LONG, root, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	for (int i=0;i<4;i++)
		sum += num_array[i];
	delete[] p;
	delete[] xpdot;
	delete[] ypdot;
	delete[] zpdot;
	delete[] updot;
	delete[] vpdot;
	delete[] wpdot;
	delete[] rhodot;
	delete[] ax;
	delete[] ay;
	delete[] az;
	delete[] ar;
	delete[] ux;
	delete[] vx;
	delete[] wx;
	delete[] aTE;
	delete[] xcor;
	delete[] ycor;
	delete[] zcor;
	delete [] cs;
	delete [] bflag;
	for (size_t i = 0; i < np_max; ++i)
	{
		for (int j = 0; j < 4; ++j)
			delete [] A_MLS[i][j];
		delete [] A_MLS[i];
	}
	delete [] A_MLS;
	delete[] rho_sum_MLS;
	delete[] sum_wab;
	delete[] rho_sum;
	delete[] beta0_MLS;
	delete[] beta1x_MLS;
	delete[] beta1y_MLS;
	delete[] beta1z_MLS;
	delete [] up;
	delete [] vp;
	delete [] wp;
	delete [] TE;
	delete [] rho;
	delete [] pm;
	delete [] upo;
	delete [] vpo;
	delete [] wpo;
	delete [] TEo;
	delete [] rhoo;
	delete [] xpo;
	delete [] ypo;
	delete [] zpo;
	if (id == root)
	{
		xpo = new double[sum];
		ypo = new double[sum];
		zpo = new double[sum];
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (id == root+2)
	{
		sum = num_array[root+2]+num_array[root+3];
		xpo = new double[sum];
		ypo = new double[sum];
		zpo = new double[sum];
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (id == root)
	{
		for (size_t i = 0; i<np-nb; i++)
		{
			xpo[i] = xp[i+nb];
			ypo[i] = yp[i+nb];
			zpo[i] = zp[i+nb];
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (id == root+2)
	{
		for (size_t i = nb; i<np; i++)
		{
			xpo[i-nb] = xp[i];
			ypo[i-nb] = yp[i];
			zpo[i-nb] = zp[i];
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (id == root || id == root+1)
	{
		if (id == root)
		{
			MPI_Status status;
			int receivenumber = num_array[root+1];
			MPI_Recv(&xpo[num_array[root]],receivenumber, MPI_DOUBLE, root+1, 1, MPI_COMM_WORLD, &status);
		}
		else
		{
			int sendnumber = number;
			MPI_Send(&xp[nb],sendnumber, MPI_DOUBLE, root, 1, MPI_COMM_WORLD);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (id == root+2 || id == root+3)
	{
		if (id == root+2)
		{
			MPI_Status status;
			int receivenumber = num_array[root+3];
			MPI_Recv(&xpo[num_array[root+2]],receivenumber, MPI_DOUBLE, root+3, 2, MPI_COMM_WORLD, &status);
		}
		else
		{
			int sendnumber = number;
			MPI_Send(&xp[nb],sendnumber, MPI_DOUBLE, root+2, 2, MPI_COMM_WORLD);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (id == root || id == root+2)
	{
		if (id == root)
		{
			MPI_Status status;
			int receivenumber = num_array[root+2]+num_array[root+3];
			MPI_Recv(&xpo[num_array[root]+num_array[root+1]],receivenumber, MPI_DOUBLE, root+2, 3, MPI_COMM_WORLD, &status);
		}
		else
		{
			int sendnumber = num_array[root+2]+num_array[root+3];
			MPI_Send(xpo,sendnumber, MPI_DOUBLE, root, 3, MPI_COMM_WORLD);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (id == root || id == root+1)
	{
		if (id == root)
		{
			MPI_Status status;
			int receivenumber = num_array[root+1];
			MPI_Recv(&ypo[num_array[root]],receivenumber, MPI_DOUBLE, root+1, 1, MPI_COMM_WORLD, &status);
		}
		else
		{
			int sendnumber = number;
			MPI_Send(&yp[nb],sendnumber, MPI_DOUBLE, root, 1, MPI_COMM_WORLD);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (id == root+2 || id == root+3)
	{
		if (id == root+2)
		{
			MPI_Status status;
			int receivenumber = num_array[root+3];
			MPI_Recv(&ypo[num_array[root+2]],receivenumber, MPI_DOUBLE, root+3, 2, MPI_COMM_WORLD, &status);
		}
		else
		{
			int sendnumber = number;
			MPI_Send(&yp[nb],sendnumber, MPI_DOUBLE, root+2, 2, MPI_COMM_WORLD);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (id == root || id == root+2)
	{
		if (id == root)
		{
			MPI_Status status;
			int receivenumber = num_array[root+2]+num_array[root+3];
			MPI_Recv(&ypo[num_array[root]+num_array[root+1]],receivenumber, MPI_DOUBLE, root+2, 3, MPI_COMM_WORLD, &status);
		}
		else
		{
			int sendnumber = num_array[root+2]+num_array[root+3];
			MPI_Send(ypo,sendnumber, MPI_DOUBLE, root, 3, MPI_COMM_WORLD);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (id == root || id == root+1)
	{
		if (id == root)
		{
			MPI_Status status;
			int receivenumber = num_array[root+1];
			MPI_Recv(&zpo[num_array[root]],receivenumber, MPI_DOUBLE, root+1, 1, MPI_COMM_WORLD, &status);
		}
		else
		{
			int sendnumber = number;
			MPI_Send(&zp[nb],sendnumber, MPI_DOUBLE, root, 1, MPI_COMM_WORLD);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (id == root+2 || id == root+3)
	{
		if (id == root+2)
		{
			MPI_Status status;
			int receivenumber = num_array[root+3];
			MPI_Recv(&zpo[num_array[root+2]],receivenumber, MPI_DOUBLE, root+3, 2, MPI_COMM_WORLD, &status);
		}
		else
		{
			int sendnumber = number;
			MPI_Send(&zp[nb],sendnumber, MPI_DOUBLE, root+2, 2, MPI_COMM_WORLD);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (id == root || id == root+2)
	{
		if (id == root)
		{
			MPI_Status status;
			int receivenumber = num_array[root+2]+num_array[root+3];
			MPI_Recv(&zpo[num_array[root]+num_array[root+1]],receivenumber, MPI_DOUBLE, root+2, 3, MPI_COMM_WORLD, &status);
		}
		else
		{
			int sendnumber = num_array[root+2]+num_array[root+3];
			MPI_Send(zpo,sendnumber, MPI_DOUBLE, root, 3, MPI_COMM_WORLD);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	delete [] scalar;
	delete [] xp;
	delete [] yp;
	delete [] zp;
	MPI_Barrier(MPI_COMM_WORLD);
	if (id == root)
	{
		scalar = new double[(N+1)*(M+1)*(L+1)];
		int dimension[3] = {N+1, M+1,L+1};
		//	scalar = new double**[n1-n0+1];
		//	for (int i = 0; i<(n1-n0+1); i++)
		//	{
		//		scalar[i] = new double*[m1-m0+1];
		//		for (int j = 0;j<(m1-m0+1); j++)
		//			scalar[i][j] = new double[l1-l0+1];
		//	}
		for (int i = 0; i<(N+1); i++)
			for (int j = 0;j<(M+1); j++)
				for (int k = 0;k<(L+1);k++)
				{
					int index = Index3D(i,j,k,(N+1),(M+1),(L+1));
					scalar[index] = -1;
				}
		double origin[3];
		origin[0] = vlx[0];
		origin[1] = vly[0];
		origin[2] = vlz[0];
		//finish mesh building

		//locate position
		for (size_t i = 0;i<sum;i++)
		{
			double distx,disty,distz; //to the origin
			distx = xpo[i] - origin[0];
			disty = ypo[i] - origin[1];
			distz = zpo[i] - origin[2];
			int icell,jcell,kcell;
			icell = (int)(distx / dx);
			jcell = (int)(disty / dy);
			kcell = (int)(distz / dz);
			double x0,x1,y0,y1,z0,z1;
			x0 = origin[0] + icell * dx;
			x1 = x0 + dx;
			y0 = origin[1] + jcell * dy;
			y1 = y0 + dy;
			z0 = origin[2] + kcell * dz;
			z1 = z0 + dz;
			double distancex[2],distancey[2],distancez[2]; //to faces
			distancex[0] = (xpo[i] - x0)/dx;
			distancex[1] = (x1 - xpo[i])/dx;
			distancey[0] = (ypo[i] - y0)/dy;
			distancey[1] = (y1 - ypo[i])/dy;
			distancez[0] = (zpo[i] - z0)/dz;
			distancez[1] = (z1 - zpo[i])/dz;
			//get volumn fraction
			double volumn[8]={0};
			//index
			//0:west-south-down; 1:east-south-down; 2:east-north-down; 3:west-north-down
			//4:west-south-up;   5:east-south-up;   6:east-north-up;   7:east-north-up
			volumn[0] = distancex[0] * distancey[0] * distancez[0];
			volumn[1] = distancex[1] * distancey[0] * distancez[0];
			volumn[2] = distancex[1] * distancey[1] * distancez[0];
			volumn[3] = distancex[0] * distancey[1] * distancez[0];
			volumn[4] = distancex[0] * distancey[0] * distancez[1];
			volumn[5] = distancex[1] * distancey[0] * distancez[1];
			volumn[6] = distancex[1] * distancey[1] * distancez[1];
			volumn[7] = distancex[0] * distancey[1] * distancez[1];
			double volumn_inverse[8];
			for (int count1 = 0; count1<8; count1++)
			{
				volumn_inverse[count1] = 1;
				for (int count2 = 0; count2<8; count2++)
				{
					if (count1 == count2) continue;
					else
						volumn_inverse[count1] = volumn_inverse[count1] * volumn[count2];
				}
			}
			double total=0;
			for (int count1 = 0; count1<8; count1++)
				total = total + volumn_inverse[count1];
			for (int count1 = 0; count1<8; count1++)
				volumn_inverse[count1] = 2 * volumn_inverse[count1] / total;

			//deposit value
			int index;
			index = Index3D(icell,jcell,kcell,(N+1),(M+1),(L+1));
			scalar[index] += volumn_inverse[0];
			index = Index3D(icell+1,jcell,kcell,(N+1),(M+1),(L+1));
			scalar[index] += volumn_inverse[1];
			index = Index3D(icell+1,jcell+1,kcell,(N+1),(M+1),(L+1));
			scalar[index] += volumn_inverse[2];
			index = Index3D(icell,jcell+1,kcell,(N+1),(M+1),(L+1));
			scalar[index] += volumn_inverse[3];
			index = Index3D(icell,jcell,kcell+1,(N+1),(M+1),(L+1));
			scalar[index] += volumn_inverse[4];
			index = Index3D(icell+1,jcell,kcell+1,(N+1),(M+1),(L+1));
			scalar[index] += volumn_inverse[5];
			index = Index3D(icell+1,jcell+1,kcell+1,(N+1),(M+1),(L+1));
			scalar[index] += volumn_inverse[6];
			index = Index3D(icell,jcell+1,kcell+1,(N+1),(M+1),(L+1));
			scalar[index] += volumn_inverse[7];
		}
		char filename[200];
		FILE *outfile;
		sprintf(filename,"vtkoutput_all_%s",right_flush(step,7));
		//		if (PN * PM * PL != 1)
		//			sprintf(filename,"%s-nd%s",filename,right_flush(id,4));
		sprintf(filename,"%s.vtk",filename);
		outfile = fopen(filename,"w");
		fprintf(outfile,"# vtk DataFile Version 3.0\n");
		fprintf(outfile,"The actual time is %.8f\n",-0.1);
		fprintf(outfile,"ASCII\n");
		fprintf(outfile,"DATASET STRUCTURED_POINTS\n");
		fprintf(outfile,"DIMENSIONS %d %d %d\n",dimension[0],dimension[1],dimension[2]);
		fprintf(outfile,"SPACING %lf %lf %lf\n",dx,dy,dz);
		fprintf(outfile,"ORIGIN %lf %lf %lf\n", origin[0],origin[1],origin[2]);
		fprintf(outfile,"POINT_DATA %d\n",(dimension[0])*(dimension[1])*(dimension[2]));
		fprintf(outfile,"SCALARS Interface double 1\n");
		fprintf(outfile,"LOOKUP_TABLE default\n");

		for (int k = 0;k<(L+1);k++)
			for (int j = 0;j<(M+1); j++)
				for (int i = 0; i<(N+1); i++)
				{
					int index = Index3D(i,j,k,(N+1),(M+1),(L+1));
					fprintf(outfile,"%.16g\n",scalar[index]);
				}
		//		cout<<"I am here"<<endl;
		//		size_t i;
		//		char filename[200];
		//		FILE *outfile;
		//		sprintf(filename,"vtkoutput_all_%s",right_flush(step,7));
		//		sprintf(filename,"%s.vtk",filename);
		//		outfile = fopen(filename,"w");
		//		double d = 0;
		//		fprintf(outfile,"# vtk DataFile Version 3.0\n");
		//		fprintf(outfile,"The actual time is %.8f\n",-0.1);
		//		fprintf(outfile,"ASCII\n");
		//		fprintf(outfile,"DATASET POLYDATA\n");
		//		fprintf(outfile,"POINTS %lu double\n",sum);
		//		for (i =0;i<sum;i++)
		//			fprintf(outfile,"%.16g %.16g %.16g\n",xpo[i],ypo[i],zpo[i]);
		//		fprintf(outfile,"POINT_DATA %lu\n",sum);
		//		fprintf(outfile,"SCALARS point double\n");
		//		fprintf(outfile,"LOOKUP_TABLE default\n");
		//		for (i = 0;i<sum;i++)
		//		{
		//			fprintf(outfile,"1\n");
		//		}
		fclose(outfile);
	}
}

void STORAGE::interior_print_out(int step, double dt,double time)
{
	size_t i;
	char filename[200];
	FILE *outfile;
	//	sprintf(filename,"interior_%s",right_flush(step,7));
	sprintf(filename,"interior");
	if (PN * PM * PL != 1)
		sprintf(filename,"%s-nd%s",filename,right_flush(id,4));
	sprintf(filename,"%s.dat",filename);
	outfile = fopen(filename,"w");
	fprintf(outfile,"%lu %lu %lu\n", nbb1, np - nb + nbb1, np_max);
	for (i = 0;i<nbb1;i++)
		fprintf(outfile,"%.24g %.24g %.24g %.24g %.24g %.24g %.24g %.24g %.24g\n", xp[i],yp[i],zp[i],up[i],vp[i],wp[i],TE[i],rho[i],pm[i]);
	for (i = nb;i<np;i++)
		fprintf(outfile,"%.24g %.24g %.24g %.24g %.24g %.24g %.24g %.24g %.24g\n", xp[i],yp[i],zp[i],up[i],vp[i],wp[i],TE[i],rho[i],pm[i]);
	fprintf(outfile,"%.24g\n",dt);
	fprintf(outfile,"%.24g\n",time);
	fprintf(outfile,"%d\n",step);
	fclose(outfile);
}

void STORAGE::interior_read(int &step, double &dt, double &time)
{
	size_t i;
	char filename[200];
	FILE *outfile;
	//	sprintf(filename,"interior_%s",right_flush(step,7));
	sprintf(filename,"interior");
	if (PN * PM * PL != 1)
		sprintf(filename,"%s-nd%s",filename,right_flush(id,4));
	sprintf(filename,"%s.dat",filename);
	outfile = fopen(filename,"r");
	delete [] xp;
	delete [] yp;
	delete [] zp;
	delete [] up;
	delete [] vp;
	delete [] wp;
	delete [] TE;
	delete [] rho;
	delete [] pm;
	delete [] xpo;
	delete [] ypo;
	delete [] zpo;
	delete [] upo;
	delete [] vpo;
	delete [] wpo;
	delete [] TEo;
	delete [] rhoo;
	delete [] bflag;
	delete[] p;
	delete[] xpdot;
	delete[] ypdot;
	delete[] zpdot;
	delete[] updot;
	delete[] vpdot;
	delete[] wpdot;
	delete[] rhodot;
	delete[] ax;
	delete[] ay;
	delete[] az;
	delete[] ar;
	delete[] ux;
	delete[] vx;
	delete[] wx;
	delete[] aTE;
	delete[] xcor;
	delete[] ycor;
	delete[] zcor;
	delete [] cs;
	for (size_t i = 0; i < np_max; ++i)
	{
		for (int j = 0; j < 4; ++j)
			delete [] A_MLS[i][j];
		delete [] A_MLS[i];
	}
	delete [] A_MLS;
	delete[] rho_sum_MLS;
	delete[] sum_wab;
	delete[] rho_sum;
	delete[] beta0_MLS;
	delete[] beta1x_MLS;
	delete[] beta1y_MLS;
	delete[] beta1z_MLS;
	fscanf(outfile,"%lu %lu %lu\n", &nbb1, &np, &np_max);
	nb = nbb1;
	nbb = nbb1;
	xp = new double[np_max];
	yp = new double[np_max];
	zp = new double[np_max];
	up = new double[np_max];
	vp = new double[np_max];
	wp = new double[np_max];
	TE = new double[np_max];
	rho = new double[np_max];
	pm = new double[np_max];
	cs = new double[np_max];
	p = new double[np_max];
	iflag = new int[np_max];
	bflag = new int[np_max];
	xpo = new double[np_max];
	ypo = new double[np_max];
	zpo = new double[np_max];
	upo = new double[np_max];
	vpo = new double[np_max];
	wpo = new double[np_max];
	TEo = new double[np_max];
	rhoo = new double[np_max];
	xpdot = new double[np_max];
	ypdot = new double[np_max];
	zpdot = new double[np_max];
	updot = new double[np_max];
	vpdot = new double[np_max];
	wpdot = new double[np_max];
	TEdot = new double[np_max];
	rhodot = new double[np_max];
	ax = new double[np_max];
	ay = new double[np_max];
	az = new double[np_max];
	ar = new double[np_max];
	ux = new double[np_max];
	vx = new double[np_max];
	wx = new double[np_max];
	aTE = new double[np_max];
	xcor = new double[np_max];
	ycor = new double[np_max];
	zcor = new double[np_max];
	A_MLS = new double**[np_max];
	for (size_t i = 0; i < np_max; ++i)
	{
		A_MLS[i] = new double*[4];
		for (int j = 0; j < 4; ++j)
			A_MLS[i][j] = new double[4];
	}
	rho_sum_MLS = new double[np_max];
	sum_wab = new double[np_max];
	rho_sum = new double[np_max];
	beta0_MLS = new double[np_max];
	beta1x_MLS = new double[np_max];
	beta1y_MLS = new double[np_max];
	beta1z_MLS = new double[np_max];
	for (size_t i=0;i<np_max;i++)
	{
		xp[i]=0;
		yp[i]=0;
		zp[i]=0;
		up[i]=0;
		vp[i]=0;
		wp[i]=0;
		TE[i]=0;
		rho[i]=0;
		pm[i]=0;
		cs[i]=0;
		p[i]=0;
		iflag[i] = 1;
		bflag[i] = 0;
		xpo[i]=0;
		ypo[i]=0;
		zpo[i]=0;
		upo[i]=0;
		vpo[i]=0;
		wpo[i]=0;
		TEo[i]=0;
		rhoo[i]=0;
		xpdot[i] = 0;
		ypdot[i] = 0;
		zpdot[i] = 0;
		updot[i] = 0;
		vpdot[i] = 0;
		wpdot[i] = 0;
		TEdot[i] = 0;
		rhodot[i] = 0;
		xpdot[i]=0;
		ypdot[i]=0;
		zpdot[i]=0;
		updot[i]=0;
		vpdot[i]=0;
		wpdot[i]=0;
		TEdot[i]=0;
		rhodot[i]=0;
		ax[i]=0;
		ay[i]=0;
		az[i]=0;
		ar[i]=0;
		ux[i]=0;
		vx[i]=0;
		wx[i]=0;
		aTE[i]=0;
		xcor[i]=0;
		ycor[i]=0;
		zcor[i]=0;
	}
	for (i = 0;i<np;i++)
		fscanf(outfile,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &xp[i],&yp[i],&zp[i],&up[i],&vp[i],&wp[i],&TE[i],&rho[i],&pm[i]);
	fscanf(outfile,"%lf\n",&dt);
	fscanf(outfile,"%lf\n",&time);
	fprintf(outfile,"%d\n",&step);
	nprint = (int)(time / out);
	next_print_time = (nprint+1) * out;
}


