#ifndef GRID1D__H
#define GRID1D__H

//---------------------------------------------------------------------
//
//---------------------------------------------------------------------
struct GRID1D
{
	double xmin,xmax,dx; int isize;
	//---------------------------------------------------------------------
	// basic grid operations
	//---------------------------------------------------------------------
	void set_grid(double xmin_,double xmax_,int isize_)
	{ xmin = xmin_;xmax = xmax_;isize= isize_;dx=(xmax-xmin)/(isize-1);}
	inline double x_fr_i(double i) const{ return xmin+dx*i;}
	inline double i_fr_x(double x) const{ return (x-xmin)/dx;}
	//---------------------------------------------------------------------
	// ENO-2 interpolation
	//---------------------------------------------------------------------
	double ENO2_interpolation( const ARRAY<double>& F, double x ) const
	{
		double i=i_fr_x(x); int i0=int(i)-1; i0=MIN(MAX(i0,0),isize-4);
		return ::ENO2_interpolation(F[i0],F[i0+1],F[i0+2],F[i0+3],i-i0);
	}
};

#endif
