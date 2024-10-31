#ifndef GEOMETRIC_INTEGRAL_H
#define GEOMETRIC_INTEGRAL_H

#include "POINT3.h"
#include "POINT2.h"
#include "../Macros.h"
#include "../ARRAY3D.h"
#include "../uniform_grid/GRID3D.h"

/*-------------------------------------------------------------------------------------*\
 * calculates integral over interface(phi=0) or over irregular domain(phi<=0)
 * on a cell of uniform grid
\*-------------------------------------------------------------------------------------*/
class GeometricIntegral
{
		
protected:
	/*-------------------------------------------------------------------------------------*\
	 * linear-interpolation based integrals over simplices
	\*-------------------------------------------------------------------------------------*/
	static double integral_over_interface_mpp(	const POINT2& P0,
												const POINT2& P1,
												const POINT2& P2, 	double phi0, double f0, 
																	double phi1, double f1, 
																	double phi2, double f2 )
	{
		double r01 = phi0/(phi0-phi1); double r10 = 1.-r01;
		double r02 = phi0/(phi0-phi2); double r20 = 1.-r02;
		
		POINT2 P01(P0.x*r10 + P1.x*r01, P0.y*r10 + P1.y*r01); double f01 = f0*r10 + f1*r01;
		POINT2 P02(P0.x*r20 + P2.x*r02, P0.y*r20 + P2.y*r02); double f02 = f0*r20 + f2*r02;
		
		POINT2& P01_02 = P02;
		P01_02.x -= P01.x;
		P01_02.y -= P01.y;
		
		return .5*(f01+f02)*sqrt(P01_02.x*P01_02.x
								+P01_02.y*P01_02.y);
	}
	
	//
	static double integral_over_interface_mppp(  const POINT3& P0,
												 const POINT3& P1,
												 const POINT3& P2,
												 const POINT3& P3, 	double phi0, double f0, 
																	double phi1, double f1, 
																	double phi2, double f2,
																	double phi3, double f3)
	{
		double r01 = phi0/(phi0-phi1); double r10 = 1.-r01;
		double r02 = phi0/(phi0-phi2); double r20 = 1.-r02;
		double r03 = phi0/(phi0-phi3); double r30 = 1.-r03;
		
		POINT3 P01(	P0.x*r10 + P1.x*r01, 
					P0.y*r10 + P1.y*r01, 
					P0.z*r10 + P1.z*r01); double f01 = f0*r10 + f1*r01;
		POINT3 P02(	P0.x*r20 + P2.x*r02, 
					P0.y*r20 + P2.y*r02, 
					P0.z*r20 + P2.z*r02); double f02 = f0*r20 + f2*r02;
		POINT3 P03(	P0.x*r30 + P3.x*r03, 
					P0.y*r30 + P3.y*r03, 
					P0.z*r30 + P3.z*r03); double f03 = f0*r30 + f3*r03;
		
		return (f01+f02+f03)/3*POINT3::area(P01,P02,P03);
	}
	
	//
	static double integral_over_interface_mmpp(  const POINT3& P0,
												 const POINT3& P1,
												 const POINT3& P2,
												 const POINT3& P3, 	double phi0, double f0, 
																	double phi1, double f1, 
																	double phi2, double f2,
																	double phi3, double f3)
	{
		double r02 = phi0/(phi0-phi2); double r20 = 1.-r02;
		double r03 = phi0/(phi0-phi3); double r30 = 1.-r03;
		double r12 = phi1/(phi1-phi2); double r21 = 1.-r12;
		double r13 = phi1/(phi1-phi3); double r31 = 1.-r13;
		
		POINT3 P02(	P0.x*r20 + P2.x*r02, 
					P0.y*r20 + P2.y*r02, 
					P0.z*r20 + P2.z*r02); double f02 = f0*r20 + f2*r02;
		POINT3 P03(	P0.x*r30 + P3.x*r03, 
					P0.y*r30 + P3.y*r03, 
					P0.z*r30 + P3.z*r03); double f03 = f0*r30 + f3*r03;
		POINT3 P12(	P1.x*r21 + P2.x*r12, 
					P1.y*r21 + P2.y*r12, 
					P1.z*r21 + P2.z*r12); double f12 = f1*r21 + f2*r12;
		POINT3 P13(	P1.x*r31 + P3.x*r13, 
					P1.y*r31 + P3.y*r13, 
					P1.z*r31 + P3.z*r13); double f13 = f1*r31 + f3*r13;
		
		return (f02+f03+f13)/3.*POINT3::area(P02,P03,P13)+
		       (f02+f12+f13)/3.*POINT3::area(P02,P12,P13);
	}
	
	//
	static double integral_over_minusdomain_mppp(const POINT3& P0,
												 const POINT3& P1,
												 const POINT3& P2,
												 const POINT3& P3, 	double phi0, double f0, 
																	double phi1, double f1, 
																	double phi2, double f2,
																	double phi3, double f3)
	{
		double r01 = phi0/(phi0-phi1); double r10 = 1.-r01;
		double r02 = phi0/(phi0-phi2); double r20 = 1.-r02;
		double r03 = phi0/(phi0-phi3); double r30 = 1.-r03;
		
		POINT3 P01(P0.x*r10 + P1.x*r01, P0.y*r10 + P1.y*r01, P0.z*r10 + P1.z*r01); double f01 = f0*r10 + f1*r01;
		POINT3 P02(P0.x*r20 + P2.x*r02, P0.y*r20 + P2.y*r02, P0.z*r20 + P2.z*r02); double f02 = f0*r20 + f2*r02;
		POINT3 P03(P0.x*r30 + P3.x*r03, P0.y*r30 + P3.y*r03, P0.z*r30 + P3.z*r03); double f03 = f0*r30 + f3*r03;
		
		return .25*(f0+f01+f02+f03)*POINT3::volume(P0,P01,P02,P03);
	}
	
	//
	static double integral_over_minusdomain_mmpp(const POINT3& P0,
												 const POINT3& P1,
												 const POINT3& P2,
												 const POINT3& P3, 	double phi0, double f0, 
																	double phi1, double f1, 
																	double phi2, double f2,
																	double phi3, double f3)
	{
		double r02 = phi0/(phi0-phi2); double r20 = 1.-r02;
		double r03 = phi0/(phi0-phi3); double r30 = 1.-r03;
		double r12 = phi1/(phi1-phi2); double r21 = 1.-r12;
		double r13 = phi1/(phi1-phi3); double r31 = 1.-r13;
		
		POINT3 P02(P0.x*r20 + P2.x*r02, P0.y*r20 + P2.y*r02, P0.z*r20 + P2.z*r02); double f02 = f0*r20 + f2*r02;
		POINT3 P03(P0.x*r30 + P3.x*r03, P0.y*r30 + P3.y*r03, P0.z*r30 + P3.z*r03); double f03 = f0*r30 + f3*r03;
		POINT3 P12(P1.x*r21 + P2.x*r12, P1.y*r21 + P2.y*r12, P1.z*r21 + P2.z*r12); double f12 = f1*r21 + f2*r12;
		POINT3 P13(P1.x*r31 + P3.x*r13, P1.y*r31 + P3.y*r13, P1.z*r31 + P3.z*r13); double f13 = f1*r31 + f3*r13;
		
		return .25*(f0+f02+f03+f13)*POINT3::volume(P0,P02,P03,P13)+
		       .25*(f1+f02+f12+f13)*POINT3::volume(P1,P02,P12,P13)+
			   .25*(f0+f1 +f02+f13)*POINT3::volume(P0,P1 ,P02,P13);
	}
	
	//
	static double integral_over_minusdomain_pmmm(const POINT3& P0,
												 const POINT3& P1,
												 const POINT3& P2,
												 const POINT3& P3, 	double phi0, double f0, 
																	double phi1, double f1, 
																	double phi2, double f2,
																	double phi3, double f3)
	{
		double r01 = phi0/(phi0-phi1); double r10 = 1.-r01;
		double r02 = phi0/(phi0-phi2); double r20 = 1.-r02;
		double r03 = phi0/(phi0-phi3); double r30 = 1.-r03;
		
		POINT3 P01(P0.x*r10 + P1.x*r01, P0.y*r10 + P1.y*r01, P0.z*r10 + P1.z*r01); double f01 = f0*r10 + f1*r01;
		POINT3 P02(P0.x*r20 + P2.x*r02, P0.y*r20 + P2.y*r02, P0.z*r20 + P2.z*r02); double f02 = f0*r20 + f2*r02;
		POINT3 P03(P0.x*r30 + P3.x*r03, P0.y*r30 + P3.y*r03, P0.z*r30 + P3.z*r03); double f03 = f0*r30 + f3*r03;
		
		return .25*(f0+f1 +f2 +f3 )*POINT3::volume(P0,P1 ,P2 ,P3 )-
		       .25*(f0+f01+f02+f03)*POINT3::volume(P0,P01,P02,P03);
	}
public:
	/*----------------------------------------------------------------------------------
	 * integral over interface in 2D
	 *--------------------------------------------------------------------------------*/
	static double integral_over_interface( double dx, 
										   double dy, 
										   double phi00, double f00,
										   double phi01, double f01,
										   double phi10, double f10,
										   double phi11, double f11 )
	{
		// perturbation 
		if(ABS(phi00)<eps_machine) phi00 = (phi00>0) ? eps_machine : -eps_machine;
		if(ABS(phi01)<eps_machine) phi01 = (phi01>0) ? eps_machine : -eps_machine;
		if(ABS(phi10)<eps_machine) phi10 = (phi10>0) ? eps_machine : -eps_machine;
		if(ABS(phi11)<eps_machine) phi11 = (phi11>0) ? eps_machine : -eps_machine;
		
		if(phi00<0 && phi01<0 && phi10<0 && phi11<0) return 0;
		if(phi00>0 && phi01>0 && phi10>0 && phi11>0) return 0;
		
		POINT2 P00( 0, 0);
		POINT2 P01( 0,dy);
		POINT2 P10(dx, 0);
		POINT2 P11(dx,dy);
		
		const POINT2& P0 = P00; const double& phi0 = phi00; const double& f0 = f00;
		const POINT2& P2 = P11; const double& phi2 = phi11; const double& f2 = f11;
		
		// Kuhn Triangulation
		double sum=0;
		for(int n=0;n<2;n++)
		{
			const POINT2&   P1 = (n==0) ?   P01 :   P10;
			const double&   f1 = (n==0) ?   f01 :   f10;
			const double& phi1 = (n==0) ? phi01 : phi10;
			
			// there are 8 types
					if(phi0<0 && phi1<0 && phi2<0) ;
			else if(phi0<0 && phi1<0 && phi2>0) sum += integral_over_interface_mpp(P2,P0,P1,phi2,f2,phi0,f0,phi1,f1);
			else if(phi0<0 && phi1>0 && phi2<0) sum += integral_over_interface_mpp(P1,P0,P2,phi1,f1,phi0,f0,phi2,f2);
			else if(phi0<0 && phi1>0 && phi2>0) sum += integral_over_interface_mpp(P0,P1,P2,phi0,f0,phi1,f1,phi2,f2);
			else if(phi0>0 && phi1<0 && phi2<0) sum += integral_over_interface_mpp(P0,P1,P2,phi0,f0,phi1,f1,phi2,f2);
			else if(phi0>0 && phi1<0 && phi2>0) sum += integral_over_interface_mpp(P1,P0,P2,phi1,f1,phi0,f0,phi2,f2);
			else if(phi0>0 && phi1>0 && phi2<0) sum += integral_over_interface_mpp(P2,P0,P1,phi2,f2,phi0,f0,phi1,f1);
			else 								;
			
		}
		
		return sum;
	}
protected:
	/*----------------------------------------------------------------------------------
	 * trilinear interpolation
	 *--------------------------------------------------------------------------------*/
	static double trilinear( double xi, double dx, double x,
							 double yj, double dy, double y,
							 double zk, double dz, double z, double f[2][2][2] )
	{
		x-=xi; x/=dx;
		y-=yj; y/=dy;
		z-=zk; z/=dz;
		return f[0][0][0]*(1-x)*(1-y)*(1-z)
			  +f[0][0][1]*(1-x)*(1-y)*(  z)
			  +f[0][1][0]*(1-x)*(  y)*(1-z)
			  +f[0][1][1]*(1-x)*(  y)*(  z)
			  +f[1][0][0]*(  x)*(1-y)*(1-z)
			  +f[1][0][1]*(  x)*(1-y)*(  z)
			  +f[1][1][0]*(  x)*(  y)*(1-z)
			  +f[1][1][1]*(  x)*(  y)*(  z);
	}

public:
	/*----------------------------------------------------------------------------------
	 * integral over interface in 3D
	 * 
	 * verified to be 2nd order, March 28th 2010
	 *--------------------------------------------------------------------------------*/
	static double integral_over_interface( const GRID3D& grid,
										   const ARRAY3D<double>& Phi, 
										   const ARRAY3D<double>& F  )
	{
		double sum=0;
		
		POINT3 P000( 0, 0, 0);
		POINT3 P001( 0, 0,dz);
		POINT3 P010( 0,dy, 0);
		POINT3 P011( 0,dy,dz);
		POINT3 P100(dx, 0, 0);
		POINT3 P101(dx, 0,dz);
		POINT3 P110(dx,dy, 0);
		POINT3 P111(dx,dy,dz);
		
		for(int i=0;i<grid.isize-1;i++)
		for(int j=0;j<grid.jsize-1;j++)
		for(int k=0;k<grid.ksize-1;k++)
		{
			if(ABS(phi000)<eps_machine) phi000 = (phi000>0) ? eps_machine : -eps_machine;
			if(ABS(phi001)<eps_machine) phi001 = (phi001>0) ? eps_machine : -eps_machine;
			if(ABS(phi010)<eps_machine) phi010 = (phi010>0) ? eps_machine : -eps_machine;
			if(ABS(phi011)<eps_machine) phi011 = (phi011>0) ? eps_machine : -eps_machine;
			if(ABS(phi100)<eps_machine) phi100 = (phi100>0) ? eps_machine : -eps_machine;
			if(ABS(phi101)<eps_machine) phi101 = (phi101>0) ? eps_machine : -eps_machine;
			if(ABS(phi110)<eps_machine) phi110 = (phi110>0) ? eps_machine : -eps_machine;
			if(ABS(phi111)<eps_machine) phi111 = (phi111>0) ? eps_machine : -eps_machine;
			                                         
			if( phi000<0 && phi001<0 && phi010<0 && phi011<0 &&
				phi100<0 && phi101<0 && phi110<0 && phi111<0 ) return 0;
			if( phi000>0 && phi001>0 && phi010>0 && phi011>0 &&
				phi100>0 && phi101>0 && phi110>0 && phi111>0 ) return 0;
					
			POINT3 P0 = P000; double phi0 = Phi(i,j,k); double f0 = f000;
			POINT3 P3 = P111; double phi3 = phi111; double f3 = f111;
			
			// Kuhn Triangulation
			double sum=0;
			for(int n=0;n<6;n++)
			{
				POINT3 P1; double phi1, f1;
				POINT3 P2; double phi2, f2;
				
					 if(n==0){ P1=P001; phi1=phi001; f1=f001; P2=P011; phi2=phi011; f2=f011; }
				else if(n==1){ P1=P001; phi1=phi001; f1=f001; P2=P101; phi2=phi101; f2=f101; }
				else if(n==2){ P1=P010; phi1=phi010; f1=f010; P2=P011; phi2=phi011; f2=f011; }
				else if(n==3){ P1=P010; phi1=phi010; f1=f010; P2=P110; phi2=phi110; f2=f110; }
				else if(n==4){ P1=P100; phi1=phi100; f1=f100; P2=P101; phi2=phi101; f2=f101; }
				else         { P1=P100; phi1=phi100; f1=f100; P2=P110; phi2=phi110; f2=f110; }
				
				// there are 16 types
					 if(phi0<0 && phi1<0 && phi2<0 && phi3<0) ;
				else if(phi0<0 && phi1<0 && phi2<0 && phi3>0) sum += integral_over_interface_mppp(P3,P0,P1,P2,phi3,f3,phi0,f0,phi1,f1,phi2,f2);
				else if(phi0>0 && phi1>0 && phi2>0 && phi3<0) sum += integral_over_interface_mppp(P3,P0,P1,P2,phi3,f3,phi0,f0,phi1,f1,phi2,f2);
				else if(phi0<0 && phi1<0 && phi2>0 && phi3<0) sum += integral_over_interface_mppp(P2,P0,P1,P3,phi2,f2,phi0,f0,phi1,f1,phi3,f3);
				else if(phi0>0 && phi1>0 && phi2<0 && phi3>0) sum += integral_over_interface_mppp(P2,P0,P1,P3,phi2,f2,phi0,f0,phi1,f1,phi3,f3);
				else if(phi0<0 && phi1>0 && phi2<0 && phi3<0) sum += integral_over_interface_mppp(P1,P0,P2,P3,phi1,f1,phi0,f0,phi2,f2,phi3,f3);
				else if(phi0>0 && phi1<0 && phi2>0 && phi3>0) sum += integral_over_interface_mppp(P1,P0,P2,P3,phi1,f1,phi0,f0,phi2,f2,phi3,f3);
				else if(phi0<0 && phi1>0 && phi2>0 && phi3>0) sum += integral_over_interface_mppp(P0,P1,P2,P3,phi0,f0,phi1,f1,phi2,f2,phi3,f3);
				else if(phi0>0 && phi1<0 && phi2<0 && phi3<0) sum += integral_over_interface_mppp(P0,P1,P2,P3,phi0,f0,phi1,f1,phi2,f2,phi3,f3);
				else if(phi0<0 && phi1<0 && phi2>0 && phi3>0) sum += integral_over_interface_mmpp(P0,P1,P2,P3,phi0,f0,phi1,f1,phi2,f2,phi3,f3);
				else if(phi0>0 && phi1>0 && phi2<0 && phi3<0) sum += integral_over_interface_mmpp(P0,P1,P2,P3,phi0,f0,phi1,f1,phi2,f2,phi3,f3);
				else if(phi0<0 && phi1>0 && phi2<0 && phi3>0) sum += integral_over_interface_mmpp(P0,P2,P1,P3,phi0,f0,phi2,f2,phi1,f1,phi3,f3);
				else if(phi0>0 && phi1<0 && phi2>0 && phi3<0) sum += integral_over_interface_mmpp(P0,P2,P1,P3,phi0,f0,phi2,f2,phi1,f1,phi3,f3);
				else if(phi0<0 && phi1>0 && phi2>0 && phi3<0) sum += integral_over_interface_mmpp(P0,P3,P1,P2,phi0,f0,phi3,f3,phi1,f1,phi2,f2);
				else if(phi0>0 && phi1<0 && phi2<0 && phi3>0) sum += integral_over_interface_mmpp(P0,P3,P1,P2,phi0,f0,phi3,f3,phi1,f1,phi2,f2);
				else                                          ;                                                                                 
			}                                                                                                                  
		}	                                                                                                                   
		return sum;   
	}                   
	
	/*----------------------------------------------------------------------------------
	 * integral over interface in 3D
	 *--------------------------------------------------------------------------------*/
	static double integral_over_minusdomain(double dx, 
										   double dy, 
										   double dz,
										   double phi000, double f000,
										   double phi001, double f001,
										   double phi010, double f010,
										   double phi011, double f011,
										   double phi100, double f100,
										   double phi101, double f101,
										   double phi110, double f110,
										   double phi111, double f111)
	{
		// perturbation 
		if(ABS(phi000)<eps_machine) phi000 = (phi000>0) ? eps_machine : -eps_machine;
		if(ABS(phi001)<eps_machine) phi001 = (phi001>0) ? eps_machine : -eps_machine;
		if(ABS(phi010)<eps_machine) phi010 = (phi010>0) ? eps_machine : -eps_machine;
		if(ABS(phi011)<eps_machine) phi011 = (phi011>0) ? eps_machine : -eps_machine;
		if(ABS(phi100)<eps_machine) phi100 = (phi100>0) ? eps_machine : -eps_machine;
		if(ABS(phi101)<eps_machine) phi101 = (phi101>0) ? eps_machine : -eps_machine;
		if(ABS(phi110)<eps_machine) phi110 = (phi110>0) ? eps_machine : -eps_machine;
		if(ABS(phi111)<eps_machine) phi111 = (phi111>0) ? eps_machine : -eps_machine;
		                                         
		
		if( phi000<0 && phi001<0 && phi010<0 && phi011<0 &&
			phi100<0 && phi101<0 && phi110<0 && phi111<0 ) return dx*dy*dz*.125*(f000+f001+f010+f011+f100+f101+f110+f111);
		if( phi000>0 && phi001>0 && phi010>0 && phi011>0 &&
			phi100>0 && phi101>0 && phi110>0 && phi111>0 ) return 0;
				
		POINT3 P000( 0, 0, 0);
		POINT3 P001( 0, 0,dz);
		POINT3 P010( 0,dy, 0);
		POINT3 P011( 0,dy,dz);
		POINT3 P100(dx, 0, 0);
		POINT3 P101(dx, 0,dz);
		POINT3 P110(dx,dy, 0);
		POINT3 P111(dx,dy,dz);
		
		POINT3 P0 = P000; double phi0 = phi000; double f0 = f000;
		POINT3 P3 = P111; double phi3 = phi111; double f3 = f111;
		
		// Kuhn Triangulation
		double sum=0;
		for(int n=0;n<6;n++)
		{
			POINT3 P1; double phi1, f1;
			POINT3 P2; double phi2, f2;
			
			if(n==0){ P1=P001; phi1=phi001; f1=f001; P2=P011; phi2=phi011; f2=f011; }
			if(n==1){ P1=P001; phi1=phi001; f1=f001; P2=P101; phi2=phi101; f2=f101; }
			if(n==2){ P1=P010; phi1=phi010; f1=f010; P2=P011; phi2=phi011; f2=f011; }
			if(n==3){ P1=P010; phi1=phi010; f1=f010; P2=P110; phi2=phi110; f2=f110; }
			if(n==4){ P1=P100; phi1=phi100; f1=f100; P2=P101; phi2=phi101; f2=f101; }
			if(n==5){ P1=P100; phi1=phi100; f1=f100; P2=P110; phi2=phi110; f2=f110; }
			
			// there are 16 types
			/*----*/      if(phi0<0 && phi1<0 && phi2<0 && phi3<0) sum += .25*(f0+f1+f2+f3)*POINT3::volume(P0,P1,P2,P3);
			/*---+*/ else if(phi0<0 && phi1<0 && phi2<0 && phi3>0) sum += integral_over_minusdomain_pmmm(P3,P0,P1,P2,phi3,f3,phi0,f0,phi1,f1,phi2,f2);
			/*+++-*/ else if(phi0>0 && phi1>0 && phi2>0 && phi3<0) sum += integral_over_minusdomain_mppp(P3,P0,P1,P2,phi3,f3,phi0,f0,phi1,f1,phi2,f2);
			/*--+-*/ else if(phi0<0 && phi1<0 && phi2>0 && phi3<0) sum += integral_over_minusdomain_pmmm(P2,P0,P1,P3,phi2,f2,phi0,f0,phi1,f1,phi3,f3);
			/*++-+*/ else if(phi0>0 && phi1>0 && phi2<0 && phi3>0) sum += integral_over_minusdomain_mppp(P2,P0,P1,P3,phi2,f2,phi0,f0,phi1,f1,phi3,f3);
			/*-+--*/ else if(phi0<0 && phi1>0 && phi2<0 && phi3<0) sum += integral_over_minusdomain_pmmm(P1,P0,P2,P3,phi1,f1,phi0,f0,phi2,f2,phi3,f3);
			/*+-++*/ else if(phi0>0 && phi1<0 && phi2>0 && phi3>0) sum += integral_over_minusdomain_mppp(P1,P0,P2,P3,phi1,f1,phi0,f0,phi2,f2,phi3,f3);
			/*-+++*/ else if(phi0<0 && phi1>0 && phi2>0 && phi3>0) sum += integral_over_minusdomain_mppp(P0,P1,P2,P3,phi0,f0,phi1,f1,phi2,f2,phi3,f3);
			/*+---*/ else if(phi0>0 && phi1<0 && phi2<0 && phi3<0) sum += integral_over_minusdomain_pmmm(P0,P1,P2,P3,phi0,f0,phi1,f1,phi2,f2,phi3,f3);
			/*--++*/ else if(phi0<0 && phi1<0 && phi2>0 && phi3>0) sum += integral_over_minusdomain_mmpp(P0,P1,P2,P3,phi0,f0,phi1,f1,phi2,f2,phi3,f3);
			/*++--*/ else if(phi0>0 && phi1>0 && phi2<0 && phi3<0) sum += integral_over_minusdomain_mmpp(P2,P3,P0,P1,phi2,f2,phi3,f3,phi0,f0,phi1,f1);
			/*-+-+*/ else if(phi0<0 && phi1>0 && phi2<0 && phi3>0) sum += integral_over_minusdomain_mmpp(P0,P2,P1,P3,phi0,f0,phi2,f2,phi1,f1,phi3,f3);
			/*+-+-*/ else if(phi0>0 && phi1<0 && phi2>0 && phi3<0) sum += integral_over_minusdomain_mmpp(P1,P3,P0,P2,phi1,f1,phi3,f3,phi0,f0,phi2,f2);
			/*-++-*/ else if(phi0<0 && phi1>0 && phi2>0 && phi3<0) sum += integral_over_minusdomain_mmpp(P0,P3,P1,P2,phi0,f0,phi3,f3,phi1,f1,phi2,f2);
			/*+--+*/ else if(phi0>0 && phi1<0 && phi2<0 && phi3>0) sum += integral_over_minusdomain_mmpp(P1,P2,P0,P3,phi1,f1,phi2,f2,phi0,f0,phi3,f3);
			/*++++*/ else                                          ;                                                                                 
		}                                                                                                                  
		                                                                                                                   
		return sum;                                                                                                        
	}             
};

#endif
