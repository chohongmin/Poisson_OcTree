#ifndef TRILINEAR_CH_H
#define TRILINEAR_CH_H

#include "../Macros.h"
#include "GRID3D.h"
#include "../ARRAY3D.h"

class TrilinearInterpolation : public Function_3D
{
public:
	       GRID3D           grid;
	const ARRAY3D<double>*  pf  ;
	
	TrilinearInterpolation( const GRID3D& grid_, const ARRAY3D<double>& f )
	{
		grid=grid_;
		pf=&f;
	}

	double operator()( double x, double y, double z) const
	{
		x=MAX(x,grid.xmin); x=MIN(x,grid.xmax);
		y=MAX(y,grid.ymin); y=MIN(y,grid.ymax);
		z=MAX(z,grid.zmin); z=MIN(z,grid.zmax);

		// move (x,y,z) -> (id,jd,kd)
		// (id,jd,kd) \in [i,i+1]x[j,j+1]x[k,k+1]
		double id = grid.i_fr_x(x); int i = (int)id; if(i<0)i=0; if(i+1>=grid.isize) i=grid.isize-2;
		double jd = grid.j_fr_y(y); int j = (int)jd; if(j<0)j=0; if(j+1>=grid.jsize) j=grid.jsize-2;
		double kd = grid.k_fr_z(z); int k = (int)kd; if(k<0)k=0; if(k+1>=grid.ksize) k=grid.ksize-2;
		
		// (id,jd,kd) \in [0,1]x[0,1]x[0,1]
		id-=i;
		jd-=j;
		kd-=k;
		
		const ARRAY3D<double>& f = *pf;
		
		return(1-id)*((1-jd)*(f[i  ][j  ][k  ]*(1-kd)+f[i  ][j  ][k+1]*(  kd))
		             +(  jd)*(f[i  ][j+1][k  ]*(1-kd)+f[i  ][j+1][k+1]*(  kd)))
		     +(  id)*((1-jd)*(f[i+1][j  ][k  ]*(1-kd)+f[i+1][j  ][k+1]*(  kd))
				     +(  jd)*(f[i+1][j+1][k  ]*(1-kd)+f[i+1][j+1][k+1]*(  kd)));
		
		
//		return f(i  ,j  ,k  )*(1-id)*(1-jd)*(1-kd)
//		      +f(i  ,j  ,k+1)*(1-id)*(1-jd)*(  kd)
//		      +f(i  ,j+1,k  )*(1-id)*(  jd)*(1-kd)
//		      +f(i  ,j+1,k+1)*(1-id)*(  jd)*(  kd)
//		      +f(i+1,j  ,k  )*(  id)*(1-jd)*(1-kd)
//		      +f(i+1,j  ,k+1)*(  id)*(1-jd)*(  kd)
//		      +f(i+1,j+1,k  )*(  id)*(  jd)*(1-kd)
//		      +f(i+1,j+1,k+1)*(  id)*(  jd)*(  kd);
		
	}
};

#endif
