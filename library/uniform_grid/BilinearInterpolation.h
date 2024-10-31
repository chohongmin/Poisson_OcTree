#ifndef TRILINEAR_H
#define TRILINEAR_H

#include "../Macros.h"
#include "GRID2D.h"
#include "../ARRAY2D.h"

class BilinearInterpolation : public Function_2D
{
public:
	       GRID2D           grid;
	const ARRAY2D<double>*  pf  ;
	
	BilinearInterpolation( const GRID2D& grid_, const ARRAY2D<double>& f )
	{
		grid=grid_;
		pf=&f;
	}

	double operator()( double x, double y ) const
	{
		x=MAX(x,grid.xmin); x=MIN(x,grid.xmax);
		y=MAX(y,grid.ymin); y=MIN(y,grid.ymax);

		// move (x,y,z) -> (id,jd,kd)
		// (id,jd,kd) \in [i,i+1]x[j,j+1]x[k,k+1]
		double id = grid.i_fr_x(x); int i = (int)id; if(i<0)i=0; if(i+1>=grid.isize) i=grid.isize-2;
		double jd = grid.j_fr_y(y); int j = (int)jd; if(j<0)j=0; if(j+1>=grid.jsize) j=grid.jsize-2;
		
		// (id,jd) \in [0,1]x[0,1]
		id-=i;
		jd-=j;
		
		const ARRAY2D<double>& f = *pf;
		
		return(1-id)*((1-jd)*f[i  ][j]+(jd)*f[i  ][j+1])
		     +(  id)*((1-jd)*f[i+1][j]+(jd)*f[i+1][j+1]);
	}
};

#endif
