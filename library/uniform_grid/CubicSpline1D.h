#ifndef CUBICSPLINE_1D_H
#define CUBICSPLINE_1D_H

#include "GRID1D.h"
class CubicSpline1D : public Function_1D
{
protected:
	GRID1D grid      ;
	ARRAY<double> f  ;
	ARRAY<double> fxx;
public:
	//
	void set( const GRID1D& grid_, const ARRAY<double>& f_ )
	{
		grid=grid_;
		f   =   f_; int N=f.size;
		fxx.resize(N);
		// symmetric tridiagonal system
		// d = 2*dx/3, a = dx/6 
		double d = 2./3.;
		double a = 1./6.;
		ARRAY<double> D(N); D[0]=d;
		for(int i=1;i<N;i++) D[i] = d-a*a/D[i-1];
		// RHS
		fxx[0]=0; fxx[N]=0; // the natural condition
		for(int i=1;i<N-1;i++) fxx[i] = (f[i+1]-2*f[i]+f[i-1]);
		// forward and backward substitution
		for(int i=  1;i< N;i++)  fxx[i]-=a*fxx[i-1]/D[i-1];
		for(int i=N-2;i>=0;i--){ fxx[i]-=a*fxx[i+1]; fxx[i]/=D[i];}
	}
	//
	double operator()(double x) const
	{
		int i = int(grid.i_fr_x(x)); i=MAX(MIN(i,grid.isize-2),0);
		double xi   = grid.x_fr_i(i  ), fi  =f[i  ], fxxi  =fxx[i  ];
		double xip1 = grid.x_fr_i(i+1), fip1=f[i+1], fxxip1=fxx[i+1];
		double dx=grid.dx;
		double A=(xip1-x)/dx; double C=A*(A*A-1)/6.;
		double B=(x  -xi)/dx; double D=B*(B*B-1)/6.;
		return A*fi+B*fip1+C*fxxi+D*fxxip1;
	}
};

#endif

