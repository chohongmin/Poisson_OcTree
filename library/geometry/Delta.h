#ifndef Delta_H
#define Delta_H

#include "../uniform_grid/GRID2D.h"
#include "../uniform_grid/GRID3D.h"
#include "../Macros.h"
#include "../ARRAY2D.h"
#include "../ARRAY3D.h"
#include "POINT2.h"
#include "POINT3.h"

#define FPHI_IN_GRID_3D(i,j) (1./(1.-phi##j/phi##i))
#define    P_IN_GRID_3D(i,j) (P##i*f##j##i+P##j*f##i##j)


/*--------------------------------------------------------------------
 * Discretization of the multi-dimensional delta function
 *------------------------------------------------------------------*/
class Delta 
{
	public:
	//---------------------------------------------------------------------------------
	// discretization of the delta function
	//---------------------------------------------------------------------------------
	static double d0( 	const POINT2& P0,
						const POINT2& P1,
						const POINT2& P2, 	double phi0, 
											double phi1, 
											double phi2) 
	{
		// perturbation
		double eps = 1E-10;
		if(ABS(phi0)<eps) phi0 = (phi0>0) ? eps : -eps;
		if(ABS(phi1)<eps) phi1 = (phi0>0) ? eps : -eps;
		if(ABS(phi2)<eps) phi2 = (phi0>0) ? eps : -eps;

		//
		if( phi0 <0 && phi1 <0 && phi2 <0) return 0; //---
		if( phi0 >0 && phi1 >0 && phi2 >0) return 0; //+++
		if( (phi0 <0 && phi1>0 && phi2>0)||         //-++
			(phi0 >0 && phi1<0 && phi2<0))          //+--
		{
			double f10 = phi1/(phi1-phi0);
			double f20 = phi2/(phi2-phi0);

			POINT2 P01 = P0*f10 + P1*(1-f10);
			POINT2 P02 = P0*f20 + P2*(1-f20); 

			return POINT2::length(P01,P02)/2*(f10+f20);
		}
		if( (phi0>0 && phi1 <0 && phi2>0)||         //+-+
			(phi0<0 && phi1 >0 && phi2<0))          //-+-
		{
			double f01 = phi0/(phi0-phi1);
			double f21 = phi2/(phi2-phi1);

			POINT2 P10 = P1*f01 + P0*(1-f01);
			POINT2 P12 = P1*f21 + P2*(1-f21); 

			return POINT2::length(P10,P12)/2*(1-f01);
		}
		if( (phi0>0 && phi1>0 && phi2<0)||         //++-,
			(phi0<0 && phi1<0 && phi2>0))          //--+
		{
			double f20 = phi2/(phi2-phi0);
			double f21 = phi2/(phi2-phi1);

			POINT2 P02 = P0*f20 + P2*(1-f20);
			POINT2 P12 = P1*f21 + P2*(1-f21); 

			return POINT2::length(P02,P12)/2*f20;
		}
		else
		{
			assert(false);
				return 0;
		}
	}

	static double d1( 	const POINT2& P0,
						const POINT2& P1,
						const POINT2& P2, 	double phi0, 
											double phi1, 
											double phi2){ return d0(P1,P0,P2,phi1,phi0,phi2);} 
	static double d2( 	const POINT2& P0,
						const POINT2& P1,
						const POINT2& P2, 	double phi0, 
											double phi1, 
											double phi2){ return d0(P2,P0,P1,phi2,phi0,phi1);}

	public:
	static void  do_in_2D( const GRID2D& grid,
					       const ARRAY2D<double>& Phi,
		        			     ARRAY2D<double>& D  )
	{
		double dx = grid.dx; int imin=0, imax=grid.isize-1;
		double dy = grid.dy; int jmin=0, jmax=grid.jsize-1;
		
		// 
		POINT2 P   (  0,  0);
		POINT2 Pip1( dx,  0);
		POINT2 Pjp1(  0, dy);
		POINT2 Pim1(-dx,  0);
		POINT2 Pjm1(  0,-dy);
		POINT2 Pip1jp1( dx, dy);
		POINT2 Pim1jm1(-dx,-dy);
		
		// iteration for each cell [i,i+1]x[j,j+1]
		for(int i=imin;i<=imax;i++)
		for(int j=jmin;j<=jmax;j++)
		{
			double& dij = D(i,j); dij=0;

			if(i>imin && j>imin) dij += d0(P,Pim1,Pim1jm1,Phi(i,j),Phi(i-1,j  ),Phi(i-1,j-1));
			if(i>imin && j>imin) dij += d0(P,Pjm1,Pim1jm1,Phi(i,j),Phi(i  ,j-1),Phi(i-1,j-1));
			if(i<imax && j<jmax) dij += d0(P,Pip1,Pip1jp1,Phi(i,j),Phi(i+1,j  ),Phi(i+1,j+1));
			if(i<imax && j<jmax) dij += d0(P,Pjp1,Pip1jp1,Phi(i,j),Phi(i  ,j+1),Phi(i+1,j+1));
			if(i>imin && j<jmax) dij += d0(P,Pim1,Pjp1   ,Phi(i,j),Phi(i-1,j  ),Phi(i  ,j+1));
			if(i<imax && j>jmin) dij += d0(P,Pjm1,Pip1   ,Phi(i,j),Phi(i  ,j-1),Phi(i+1,j  ));

			dij/=(dx*dy);
		} 
	}


	/*-----------------------------------------------------------------------*\
	 * 3D
	\*-----------------------------------------------------------------------*/
	static double D0( const POINT3& P0,
					  const POINT3& P1,
					  const POINT3& P2, 
					  const POINT3& P3, double phi0, 
										double phi1, 
										double phi2,
										double phi3) 
	{
		// perturbation
		const double eps = 1E-10;

		if(ABS(phi0)<eps) phi0 = (phi0>0) ? eps : -eps;
		if(ABS(phi1)<eps) phi1 = (phi1>0) ? eps : -eps;
		if(ABS(phi2)<eps) phi2 = (phi2>0) ? eps : -eps;
		if(ABS(phi3)<eps) phi3 = (phi3>0) ? eps : -eps;

		// cases
		int case_number	= 8*((phi0>0)?1:0)
						+ 4*((phi1>0)?1:0)
						+ 2*((phi2>0)?1:0)
						+   ((phi3>0)?1:0);

		switch(case_number)
		{
			case 0 : /*----*/ 
			case 15: /*++++*/ 
				return 0;
			case 1 : //---+
			case 14: //+++-
			{ 
				double f30 = FPHI_IN_GRID_3D(3,0); double f03 = 1-f30;	POINT3 P03 = P_IN_GRID_3D(0,3);
				double f31 = FPHI_IN_GRID_3D(3,1); double f13 = 1-f31;	POINT3 P13 = P_IN_GRID_3D(1,3);
				double f32 = FPHI_IN_GRID_3D(3,2); double f23 = 1-f32;	POINT3 P23 = P_IN_GRID_3D(2,3);
				
				return  POINT3::area(P03,P13,P23)*f30/3.;
			}
			case 2 : //--+-
			case 13: //++-+
			{ 
				double f20 = FPHI_IN_GRID_3D(2,0); double f02 = 1-f20;	POINT3 P02 = P_IN_GRID_3D(0,2);
				double f21 = FPHI_IN_GRID_3D(2,1); double f12 = 1-f21;	POINT3 P12 = P_IN_GRID_3D(1,2);
				double f23 = FPHI_IN_GRID_3D(2,3); double f32 = 1-f23;	POINT3 P32 = P_IN_GRID_3D(3,2);
				
				return  POINT3::area(P02,P12,P32)*f20/3.;
			}
			case 4 : //-+--
			case 11: //+-++
			{ 
				double f10 = FPHI_IN_GRID_3D(1,0); double f01 = 1-f10;	POINT3 P01 = P_IN_GRID_3D(0,1);
				double f12 = FPHI_IN_GRID_3D(1,2); double f21 = 1-f12;	POINT3 P21 = P_IN_GRID_3D(2,1);
				double f13 = FPHI_IN_GRID_3D(1,3); double f31 = 1-f13;	POINT3 P31 = P_IN_GRID_3D(3,1);
				
				return  POINT3::area(P01,P21,P31)*f10/3.;
			}
			case 8: //+---
			case 7: //-+++
			{ 
				double f01 = FPHI_IN_GRID_3D(0,1); double f10 = 1-f01;	POINT3 P10 = P_IN_GRID_3D(1,0);
				double f02 = FPHI_IN_GRID_3D(0,2); double f20 = 1-f02;	POINT3 P20 = P_IN_GRID_3D(2,0);
				double f03 = FPHI_IN_GRID_3D(0,3); double f30 = 1-f03;	POINT3 P30 = P_IN_GRID_3D(3,0);
				
				return  POINT3::area(P10,P20,P30)*(f10+f20+f30)/3.;
			}
			case 3: //--++
			case 12://++--
			{
				double f02 = FPHI_IN_GRID_3D(0,2); double f20 = 1-f02;	POINT3 P02 = P_IN_GRID_3D(0,2);
				double f03 = FPHI_IN_GRID_3D(0,3); double f30 = 1-f03;	POINT3 P03 = P_IN_GRID_3D(0,3);
				double f12 = FPHI_IN_GRID_3D(1,2); double f21 = 1-f12;	POINT3 P12 = P_IN_GRID_3D(1,2);
				double f13 = FPHI_IN_GRID_3D(1,3); double f31 = 1-f13;	POINT3 P13 = P_IN_GRID_3D(1,3);
				
				return  POINT3::area(P02,P03,P13)*(f20+f30)/3.
					  + POINT3::area(P02,P12,P13)*(f20    )/3.;
			}
			case 5: //-+-+
			case 10://+-+-
			{
				double f01 = FPHI_IN_GRID_3D(0,1); double f10 = 1-f01;	POINT3 P01 = P_IN_GRID_3D(0,1);
				double f03 = FPHI_IN_GRID_3D(0,3); double f30 = 1-f03;	POINT3 P03 = P_IN_GRID_3D(0,3);
				double f21 = FPHI_IN_GRID_3D(2,1); double f12 = 1-f21;	POINT3 P21 = P_IN_GRID_3D(2,1);
				double f23 = FPHI_IN_GRID_3D(2,3); double f32 = 1-f23;	POINT3 P23 = P_IN_GRID_3D(2,3);
				
				return  POINT3::area(P01,P03,P23)*(f10+f30)/3.
					  + POINT3::area(P01,P21,P23)*(f10    )/3.;
			}
			case 6: //-++-
			case 9: //+--+
			{
				double f01 = FPHI_IN_GRID_3D(0,1); double f10 = 1-f01;	POINT3 P01 = P_IN_GRID_3D(0,1);
				double f02 = FPHI_IN_GRID_3D(0,2); double f20 = 1-f02;	POINT3 P02 = P_IN_GRID_3D(0,2);
				double f31 = FPHI_IN_GRID_3D(3,1); double f13 = 1-f31;	POINT3 P31 = P_IN_GRID_3D(3,1);
				double f32 = FPHI_IN_GRID_3D(3,2); double f23 = 1-f32;	POINT3 P32 = P_IN_GRID_3D(3,2);
				
				return  POINT3::area(P01,P02,P32)*(f10+f20)/3.
					  + POINT3::area(P01,P31,P32)*(f10    )/3.;
			}
			default:
				assert(false);
				return 0;
		}	
	}
	
	static	void do_in_3D(  const GRID3D        & grid,
		                    const ARRAY3DwithBUFFER<double>& Phi, 
								  ARRAY3D          <double>& D ) 
	{
		double dx = grid.dx; int imin=0, imax=grid.isize-1;
		double dy = grid.dy; int jmin=0, jmax=grid.jsize-1;
		double dz = grid.dz; int kmin=0, kmax=grid.ksize-1;

		D.resize(Phi.isize,Phi.jsize,Phi.ksize);
		
		//---------------------------------------------------------------------------------
		// Neighborhood
		//---------------------------------------------------------------------------------
		POINT3 Pmmm(-dx,-dy,-dz);	POINT3 P0mm(  0,-dy,-dz);	POINT3 Ppmm( dx,-dy,-dz);
		POINT3 Pmm0(-dx,-dy,  0);	POINT3 P0m0(  0,-dy,  0);	POINT3 Ppm0( dx,-dy,  0);
		POINT3 Pmmp(-dx,-dy, dz);	POINT3 P0mp(  0,-dy, dz);	POINT3 Ppmp( dx,-dy, dz);
		POINT3 Pm0m(-dx,  0,-dz);	POINT3 P00m(  0,  0,-dz);	POINT3 Pp0m( dx,  0,-dz);
		POINT3 Pm00(-dx,  0,  0);	POINT3 P000(  0,  0,  0);	POINT3 Pp00( dx,  0,  0);
		POINT3 Pm0p(-dx,  0, dz);	POINT3 P00p(  0,  0, dz);	POINT3 Pp0p( dx,  0, dz);
		POINT3 Pmpm(-dx, dy,-dz);	POINT3 P0pm(  0, dy,-dz);	POINT3 Pppm( dx, dy,-dz);
		POINT3 Pmp0(-dx, dy,  0);	POINT3 P0p0(  0, dy,  0);	POINT3 Ppp0( dx, dy,  0);
		POINT3 Pmpp(-dx, dy, dz);	POINT3 P0pp(  0, dy, dz);	POINT3 Pppp( dx, dy, dz);

		for(int i=imin;i<=imax;i++)
		for(int j=jmin;j<=jmax;j++)
		for(int k=kmin;k<=kmax;k++)
		{
			double sum = 0;

			//---------------------------------------------------------------------------------
			// Neighborhood
			//---------------------------------------------------------------------------------
			double phimmm=Phi(i-1,j-1,k-1), phi0mm=Phi(i,j-1,k-1), phipmm=Phi(i+1,j-1,k-1);
			double phimm0=Phi(i-1,j-1,k  ), phi0m0=Phi(i,j-1,k  ), phipm0=Phi(i+1,j-1,k  );
			double phimmp=Phi(i-1,j-1,k+1), phi0mp=Phi(i,j-1,k+1), phipmp=Phi(i+1,j-1,k+1);
			double phim0m=Phi(i-1,j  ,k-1), phi00m=Phi(i,j  ,k-1), phip0m=Phi(i+1,j  ,k-1);
			double phim00=Phi(i-1,j  ,k  ), phi000=Phi(i,j  ,k  ), phip00=Phi(i+1,j  ,k  );
			double phim0p=Phi(i-1,j  ,k+1), phi00p=Phi(i,j  ,k+1), phip0p=Phi(i+1,j  ,k+1);
			double phimpm=Phi(i-1,j+1,k-1), phi0pm=Phi(i,j+1,k-1), phippm=Phi(i+1,j+1,k-1);
			double phimp0=Phi(i-1,j+1,k  ), phi0p0=Phi(i,j+1,k  ), phipp0=Phi(i+1,j+1,k  );
			double phimpp=Phi(i-1,j+1,k+1), phi0pp=Phi(i,j+1,k+1), phippp=Phi(i+1,j+1,k+1);

			if(i+1<=imax && j+1<=jmax && k+1<=kmax) sum+=D0(P000,Pp00,P0p0,P00p,phi000,phip00,phi0p0,phi00p);
			if(i+1<=imax && j-1>=jmin && k-1>=kmin) sum+=D0(P000,Pp00,P00m,P0m0,phi000,phip00,phi00m,phi0m0);
			if(i-1>=imin && j+1<=jmax && k-1>=kmin) sum+=D0(P000,P00m,P0p0,Pm00,phi000,phi00m,phi0p0,phim00);
			if(i-1>=imin && j-1>=jmin && k+1<=kmax) sum+=D0(P000,P0m0,Pm00,P00p,phi000,phi0m0,phim00,phi00p);
			if(i-1>=imin && j-1>=jmin && k-1>=kmin){sum+=D0(P000,P0mm,Pm0m,Pmm0,phi000,phi0mm,phim0m,phimm0);
													sum+=D0(P000,P0mm,Pm0m,P00m,phi000,phi0mm,phim0m,phi00m);
													sum+=D0(P000,P0mm,P0m0,Pmm0,phi000,phi0mm,phi0m0,phimm0);
													sum+=D0(P000,Pm00,Pm0m,Pmm0,phi000,phim00,phim0m,phimm0);}
			if(i-1>=imin && j+1<=jmax && k+1<=kmax){sum+=D0(P000,Pm00,Pmp0,Pm0p,phi000,phim00,phimp0,phim0p);
													sum+=D0(P000,P0pp,Pmp0,Pm0p,phi000,phi0pp,phimp0,phim0p);
													sum+=D0(P000,P0pp,Pmp0,P0p0,phi000,phi0pp,phimp0,phi0p0);
													sum+=D0(P000,P0pp,P00p,Pm0p,phi000,phi0pp,phi00p,phim0p);}
			if(i+1<=imax && j-1>=jmin && k+1<=kmax){sum+=D0(P000,Ppm0,P0m0,P0mp,phi000,phipm0,phi0m0,phi0mp);
													sum+=D0(P000,Ppm0,Pp0p,P0mp,phi000,phipm0,phip0p,phi0mp);
													sum+=D0(P000,Ppm0,Pp0p,Pp00,phi000,phipm0,phip0p,phip00);
													sum+=D0(P000,P00p,Pp0p,P0mp,phi000,phi00p,phip0p,phi0mp);}
			if(i+1<=imax && j+1<=jmax && k-1>=kmin){sum+=D0(P000,Pp0m,P0pm,P00m,phi000,phip0m,phi0pm,phi00m);
													sum+=D0(P000,Pp0m,P0pm,Ppp0,phi000,phip0m,phi0pm,phipp0);
													sum+=D0(P000,Pp0m,Pp00,Ppp0,phi000,phip0m,phip00,phipp0);
													sum+=D0(P000,P0p0,P0pm,Ppp0,phi000,phi0p0,phi0pm,phipp0);}
			D(i,j,k)=sum/(dx*dy*dz);
		}	
	}
}; 
 
#endif
              
			  
