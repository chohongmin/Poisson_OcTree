#ifndef HEAVISIDE2D_H
#define HEAVISIDE2D_H

#include "../ARRAY2D.h"
#include "../uniform_grid/GRID2D.h"
/*----------------------------------------------------------------------------*\
 * calculate the Heaviside function 
\*----------------------------------------------------------------------------*/
static void calculate_H_Towers(const ARRAY2DwithBUFFER<double>& Phi, 
							         ARRAY2D          <double>& Hiphjph)
{
	int isize = Phi.isize; 
	int jsize = Phi.jsize; 
	Hiphjph.resize(isize-1,jsize-1); 
	// 
	for(int i=0;i<isize-1;i++)
	for(int j=0;j<jsize-1;j++)
	{
		double phi_imhjmh = Phi(i  ,j  );
		double phi_imhjph = Phi(i  ,j+1);
		double phi_iphjmh = Phi(i+1,j  );
		double phi_iphjph = Phi(i+1,j+1);

		     if( phi_imhjmh<=0 && phi_imhjph<=0 && 
			     phi_iphjmh<=0 && phi_iphjph<=0 ) Hiphjph(i,j)=0;
		else if( phi_imhjmh> 0 && phi_imhjph> 0 && 
			     phi_iphjmh> 0 && phi_iphjph> 0 ) Hiphjph(i,j)=1;  
		else
		{
			double r_imhjmh = (phi_imhjmh>0) ? phi_imhjmh : 0;
			double r_imhjph = (phi_imhjph>0) ? phi_imhjph : 0;
			double r_iphjmh = (phi_iphjmh>0) ? phi_iphjmh : 0;
			double r_iphjph = (phi_iphjph>0) ? phi_iphjph : 0;

			Hiphjph(i,j) = ((  r_iphjph-  r_imhjmh)*(phi_iphjph-phi_imhjmh)
				           +(  r_iphjmh-  r_imhjph)*(phi_iphjmh-phi_imhjph))/
						   ((phi_iphjph-phi_imhjmh)*(phi_iphjph-phi_imhjmh)
				           +(phi_iphjmh-phi_imhjph)*(phi_iphjmh-phi_imhjph)+1E-18);
		}
	}
}

static void calculate_H_Towers(const GRID2D& grid,
							   const ARRAY2DwithBUFFER<double>& Phi, 
							         ARRAY2DwithBUFFER<double>& Hiph,
									 ARRAY2DwithBUFFER<double>& Hjph)
{
	int isize = Phi.isize; double dx=grid.dx;
	int jsize = Phi.jsize; double dy=grid.dy;
	Hiph.resize(isize-1,jsize  ,2); Hiph=0;
	Hjph.resize(isize  ,jsize-1,2); Hjph=0;

	// 
	for(int i=0;i<isize-1;i++)
	for(int j=1;j<jsize-1;j++)
	{
		double phi_imh = Phi(i  ,j);
		double phi_iph = Phi(i+1,j);
		double phi_jmh =.25*(Phi(i,j)+Phi(i+1,j)+Phi(i,j-1)+Phi(i+1,j-1));
		double phi_jph =.25*(Phi(i,j)+Phi(i+1,j)+Phi(i,j+1)+Phi(i+1,j+1));

		     if( phi_imh<=0 && phi_iph<=0 && 
			     phi_jmh<=0 && phi_jph<=0 ) Hiph(i,j)=0;
		else if( phi_imh> 0 && phi_iph> 0 && 
			     phi_jmh> 0 && phi_jph> 0 ) Hiph(i,j)=1;  
		else
		{
			double r_imh = (phi_imh>0) ? phi_imh : 0;
			double r_iph = (phi_iph>0) ? phi_iph : 0;
			double r_jmh = (phi_jmh>0) ? phi_jmh : 0;
			double r_jph = (phi_jph>0) ? phi_jph : 0;

			Hiph(i,j) = ((  r_iph-  r_imh)/dx*(phi_iph-phi_imh)/dx 
				        +(  r_jph-  r_jmh)/dy*(phi_jph-phi_jmh)/dy)/
						((phi_iph-phi_imh)/dx*(phi_iph-phi_imh)/dx
						+(phi_jph-phi_jmh)/dy*(phi_jph-phi_jmh)/dy);
		}
	}

	// 
	for(int i=1;i<isize-1;i++)
	for(int j=0;j<jsize-1;j++)
	{
		double phi_jmh = Phi(i,j  );
		double phi_jph = Phi(i,j+1);
		double phi_imh =.25*(Phi(i,j)+Phi(i,j+1)+Phi(i-1,j)+Phi(i-1,j+1));
		double phi_iph =.25*(Phi(i,j)+Phi(i,j+1)+Phi(i+1,j)+Phi(i+1,j+1));

		     if( phi_imh<=0 && phi_iph<=0 && 
			     phi_jmh<=0 && phi_jph<=0 ) Hjph(i,j)=0;
		else if( phi_imh> 0 && phi_iph> 0 && 
			     phi_jmh> 0 && phi_jph> 0 ) Hjph(i,j)=1;  
		else
		{
			double r_imh = (phi_imh>0) ? phi_imh : 0;
			double r_iph = (phi_iph>0) ? phi_iph : 0;
			double r_jmh = (phi_jmh>0) ? phi_jmh : 0;
			double r_jph = (phi_jph>0) ? phi_jph : 0;

			Hjph(i,j) = ((  r_iph-  r_imh)/dx*(phi_iph-phi_imh)/dx 
				        +(  r_jph-  r_jmh)/dy*(phi_jph-phi_jmh)/dy)/
						((phi_iph-phi_imh)/dx*(phi_iph-phi_imh)/dx
						+(phi_jph-phi_jmh)/dy*(phi_jph-phi_jmh)/dy);
		}
	}
}


/*----------------------------------------------------------------------------*\
 * calculate the Heaviside function 
\*----------------------------------------------------------------------------*/
static void calculate_H_as_Lengthfraction(const GRID2D& grid,
	                                      const ARRAY2D<double>& Phi, 
										        ARRAY2D<double>& Hiph,
												ARRAY2D<double>& Hjph)
{
	int isize = Phi.isize; double dx=grid.dx;
	int jsize = Phi.jsize; double dy=grid.dy;
	Hiph.resize(isize-1,jsize  );
	Hjph.resize(isize  ,jsize-1);

	// 
	for(int i=0;i<isize-1;i++)
	for(int j=1;j<jsize-1;j++)
	{
		double phi0 =.25*(Phi(i,j)+Phi(i+1,j)+Phi(i,j-1)+Phi(i+1,j-1));
		double phi1 =.25*(Phi(i,j)+Phi(i+1,j)+Phi(i,j+1)+Phi(i+1,j+1));
		     if(phi0<=0 && phi1<=0) Hiph(i,j)=0;
		else if(phi0> 0 && phi1<=0) Hiph(i,j)=phi0/(phi0-phi1);
		else if(phi0<=0 && phi1> 0) Hiph(i,j)=phi1/(phi1-phi0);
		else                        Hiph(i,j)=1;
	}

	// 
	for(int i=1;i<isize-1;i++)
	for(int j=0;j<jsize-1;j++)
	{
		double phi0 =.25*(Phi(i,j)+Phi(i,j+1)+Phi(i-1,j)+Phi(i-1,j+1));
		double phi1 =.25*(Phi(i,j)+Phi(i,j+1)+Phi(i+1,j)+Phi(i+1,j+1));
		     if(phi0<=0 && phi1<=0) Hjph(i,j)=0;
		else if(phi0> 0 && phi1<=0) Hjph(i,j)=phi0/(phi0-phi1);
		else if(phi0<=0 && phi1> 0) Hjph(i,j)=phi1/(phi1-phi0);
		else                        Hjph(i,j)=1;
	}
}

/*----------------------------------------------------------------------------*\
 * calculate the Heaviside function by the area_positive_region
\*----------------------------------------------------------------------------*/
static double area_fraction_positive_region( double phi0, double phi1, double phi2 )
{
	const double eps = 1E-10;
	if(ABS(phi0)<eps) phi0 = (phi0>0) ? eps : -eps;
	if(ABS(phi1)<eps) phi1 = (phi1>0) ? eps : -eps;
	if(ABS(phi2)<eps) phi2 = (phi2>0) ? eps : -eps;
	     if(phi0<0 && phi1<0 && phi2<0) return 0;
	else if(phi0<0 && phi1<0 && phi2>0) return   phi2*phi2/(phi2-phi0)/(phi2-phi1);
	else if(phi0<0 && phi1>0 && phi2<0) return   phi1*phi1/(phi1-phi0)/(phi1-phi2);
	else if(phi0<0 && phi1>0 && phi2>0) return 1-phi0*phi0/(phi0-phi1)/(phi0-phi2);
	else if(phi0>0 && phi1<0 && phi2<0) return   phi0*phi0/(phi0-phi1)/(phi0-phi2);
	else if(phi0>0 && phi1<0 && phi2>0) return 1-phi1*phi1/(phi1-phi0)/(phi1-phi2);
	else if(phi0>0 && phi1>0 && phi2<0) return 1-phi2*phi2/(phi2-phi0)/(phi2-phi1);
	else                                return 1;
}

static double area_fraction_positive_region( double phi00, 
											 double phi01, 
											 double phi10,
											 double phi11 )
{
	return area_fraction_positive_region(phi00,phi01,phi11)/2
		  +area_fraction_positive_region(phi00,phi10,phi11)/2;
}

static void calculate_H_as_subcellarea(const GRID2D& grid,
	                                   const ARRAY2D<double>& Phi, 
										     ARRAY2D<double>& Hiph,
											 ARRAY2D<double>& Hjph)
{
	int isize = Phi.isize; double dx=grid.dx;
	int jsize = Phi.jsize; double dy=grid.dy;
	Hiph.resize(isize-1,jsize  );
	Hjph.resize(isize  ,jsize-1);

	// 
	for(int i=0;i<isize-1;i++)
	for(int j=1;j<jsize-1;j++)
	{
		double phi00 =.5*(Phi(i  ,j  )+Phi(i  ,j-1));
		double phi01 =.5*(Phi(i  ,j+1)+Phi(i  ,j  ));
		double phi10 =.5*(Phi(i+1,j  )+Phi(i+1,j-1));
		double phi11 =.5*(Phi(i+1,j+1)+Phi(i+1,j  ));

		Hiph(i,j)=area_fraction_positive_region(phi00,phi01,phi10,phi11);
	}

	// 
	for(int i=1;i<isize-1;i++)
	for(int j=0;j<jsize-1;j++)
	{
		double phi00 =.5*(Phi(i  ,j  )+Phi(i-1,j  ));
		double phi01 =.5*(Phi(i+1,j  )+Phi(i  ,j  ));
		double phi10 =.5*(Phi(i  ,j+1)+Phi(i-1,j+1));
		double phi11 =.5*(Phi(i+1,j+1)+Phi(i  ,j+1));

		Hjph(i,j)=area_fraction_positive_region(phi00,phi01,phi10,phi11);
	}
}

#endif
