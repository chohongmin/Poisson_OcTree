#ifndef SIMPLICIAL_ISOSURFACING_H
#define SIMPLICIAL_ISOSURFACING_H

#include "../geometry/POINT2.h"
#include "../geometry/POINT3.h"
#include "../Macros.h"
#include "GRID3D.h"

#define INTERFACE_2D_L POINT2(xi  ,(yj*phi01-yjp1*phi00)/(phi01-phi00))
#define INTERFACE_2D_R POINT2(xip1,(yj*phi11-yjp1*phi10)/(phi11-phi10))
#define INTERFACE_2D_B POINT2((xi*phi10-xip1*phi00)/(phi10-phi00),yj  )
#define INTERFACE_2D_T POINT2((xi*phi11-xip1*phi01)/(phi11-phi01),yjp1)

class SimplicialIsosurfacing
{
	public:
	/*--------------------------------------------------------------------------------
	* find the isosurface, phi=0 on [xi,xip1]x[yj,yjp1]
	* the isosurface is represented as a union of line-segments
	* return value is the number of triangles,n 
	* (Q[2*i], Q[2*i+1]) is the ith triangle with normal n[i]
	* in maximum Q[4], n[2]
	*------------------------------------------------------------------------------*/
	static int isosurface( double xi, double xip1,
						   double yj, double yjp1, double phi00, double phi01, 
												   double phi10, double phi11, POINT2* Qs,
																			   POINT2* ns)
	{
		// perturbation 
		const double eps=1E-15;
		if(ABS(phi00)<eps) phi00 = (phi00>0)?eps:-eps;
		if(ABS(phi01)<eps) phi01 = (phi01>0)?eps:-eps;
		if(ABS(phi10)<eps) phi10 = (phi10>0)?eps:-eps;
		if(ABS(phi11)<eps) phi11 = (phi11>0)?eps:-eps;
		
		// enforce that phi11>0
		if(phi11<0)
		{
			phi11=-phi11;
			phi00=-phi00;
			phi01=-phi01;
			phi10=-phi10;
		}

		// now there are 8 cases
		int type=4*(phi10>0?1:0)
			    +2*(phi01>0?1:0)
			    +  (phi00>0?1:0);

		if(type==7) return 0;
		
		// normal
		POINT2 N=POINT2((phi11+phi10-phi01-phi00)/(2*(xip1-xi)),
			            (phi11-phi10+phi01-phi00)/(2*(yjp1-yj))); 
		
		switch(type)
		{
		case 0://---+, 00,01,10,11
			*Qs=INTERFACE_2D_T; Qs++;
			*Qs=INTERFACE_2D_R; Qs++;
			*ns=N;              ns++;
			return 1;
		case 1://+--+
			*Qs=INTERFACE_2D_B; Qs++;
			*Qs=INTERFACE_2D_L; Qs++;
			*Qs=INTERFACE_2D_L; Qs++;
			*Qs=INTERFACE_2D_T; Qs++;
			*Qs=INTERFACE_2D_T; Qs++;
			*Qs=INTERFACE_2D_R; Qs++;
			*Qs=INTERFACE_2D_R; Qs++;
			*Qs=INTERFACE_2D_B; Qs++;
			*ns=N;              ns++;
			*ns=N;              ns++;
			*ns=N;              ns++;
			*ns=N;              ns++;
			return 4;
		case 2://-+-+
			*Qs=INTERFACE_2D_L; Qs++;
			*Qs=INTERFACE_2D_R; Qs++;
			*ns=N;              ns++;
			return 1;
		case 3://++-+
			*Qs=INTERFACE_2D_B; Qs++;
			*Qs=INTERFACE_2D_R; Qs++;
			*ns=N;              ns++;
			return 1;
		case 4://--++
			*Qs=INTERFACE_2D_B; Qs++;
			*Qs=INTERFACE_2D_T; Qs++;
			*ns=N;              ns++;
			return 1;
		case 5://+-++
			*Qs=INTERFACE_2D_T; Qs++;
			*Qs=INTERFACE_2D_L; Qs++;
			*ns=N;              ns++;
			return 1;
		case 6://-+++
			*Qs=INTERFACE_2D_B; Qs++;
			*Qs=INTERFACE_2D_L; Qs++;
			*ns=N;              ns++;
			return 1;
		default://++++
			return 0;
		};
	}
	
	/*--------------------------------------------------------------------------------
	* isosurfacing on uniform grid
	* Vs : vertices of the triangulation
	* Ns : normals  of the triangulation
	* Thus size(Vs)=2*size(Ns)
	*------------------------------------------------------------------------------*/
	static void isosurface( const GRID2D& grid, 
	                        const ARRAY2D<double>& Phi, 
								ARRAY<POINT2>& Vs_arr, 
							    ARRAY<POINT2>& Ns_arr,
								   double level=0) {
		int size=0, buff_size=1000;
		POINT2* Vs = new POINT2[buff_size*2];
		POINT2* Ns = new POINT2[buff_size  ];
		POINT2 Vs_ij[8];
		POINT2 Ns_ij[4];
		
		int isize = grid.isize; double dx = grid.dx;
		int jsize = grid.jsize; double dy = grid.dy;
		for(int i=0;i<isize-1;i++){ double xi=grid.x_fr_i(i);
		for(int j=0;j<jsize-1;j++){ double yj=grid.y_fr_j(j);
			int size_ij =isosurface(xi,xi+dx,
									yj,yj+dy,Phi[i  ][j]-level,Phi[i  ][j+1]-level,
									         Phi[i+1][j]-level,Phi[i+1][j+1]-level,Vs_ij,Ns_ij);
			if(size+size_ij>=buff_size){
				POINT2* Vs_temp = new POINT2[buff_size*2*2]; 
				POINT2* Ns_temp = new POINT2[buff_size*2  ];
				memcpy(Vs_temp,Vs,sizeof(POINT2)*buff_size*2); delete[] Vs; Vs=Vs_temp;
				memcpy(Ns_temp,Ns,sizeof(POINT2)*buff_size  ); delete[] Ns; Ns=Ns_temp; buff_size*=2;
			}
			
			memcpy(Vs+2*size,Vs_ij,sizeof(POINT2)*size_ij*2);
			memcpy(Ns+  size,Ns_ij,sizeof(POINT2)*size_ij  ); size+=size_ij;
		}}
		
		Vs_arr.resize(size*2); memcpy((POINT2*)Vs_arr,Vs,sizeof(POINT2)*size*2); delete[]Vs;
		Ns_arr.resize(size  ); memcpy((POINT2*)Ns_arr,Ns,sizeof(POINT2)*size  ); delete[]Ns;
	}
	
	/*------------------------------------------------------------------------------------*\
	 * 3D isosurfacing
	\*------------------------------------------------------------------------------------*/
	protected:
	//
	static void isosurfacing_mppp(const POINT3& P0, double phi0,
								  const POINT3& P1,	double phi1,
								  const POINT3& P2,	double phi2,
								  const POINT3& P3,	double phi3, POINT3*& Qs)
	{
		double r01 = phi0/(phi0-phi1); double r10 = 1.-r01;
		double r02 = phi0/(phi0-phi2); double r20 = 1.-r02;
		double r03 = phi0/(phi0-phi3); double r30 = 1.-r03;
		
		(*Qs)=P0*r10 + P1*r01; Qs++; 
		(*Qs)=P0*r20 + P2*r02; Qs++;
		(*Qs)=P0*r30 + P3*r03; Qs++;
	}

	//
	static void isosurfacing_mmpp(const POINT3& P0, double phi0,
								  const POINT3& P1,	double phi1,
								  const POINT3& P2,	double phi2,
								  const POINT3& P3,	double phi3, POINT3*& Qs)
	{
		double r02 = phi0/(phi0-phi2); double r20 = 1.-r02;
		double r03 = phi0/(phi0-phi3); double r30 = 1.-r03;
		double r12 = phi1/(phi1-phi2); double r21 = 1.-r12;
		double r13 = phi1/(phi1-phi3); double r31 = 1.-r13;
		
		(*Qs)=P2*r02 + P0*r20; Qs++; 
		(*Qs)=P3*r03 + P0*r30; Qs++;
		(*Qs)=P3*r13 + P1*r31; Qs++;
		
		(*Qs)=P2*r02 + P0*r20; Qs++; 
		(*Qs)=P2*r12 + P1*r21; Qs++;
		(*Qs)=P3*r13 + P1*r31; Qs++;
	}
	
	public:
	/*--------------------------------------------------------------------------------
	* find the isosurface, phi=0 on [xi,xip1]x[yj,yjp1]x[zk,zkp1]
	* the isosurface is represented as a union of triangles
	* return value is the number of triangles,n 
	* (Q[3*i], Q[3*i+1],Q[3*i+2]) is the ith triangle with normal n[i]
	* in maximum Q[36], n[12]
	*------------------------------------------------------------------------------*/
	static int isosurface( double xi, double xip1,
						   double yj, double yjp1,
						   double zk, double zkp1, double phi000, double phi001, 
												   double phi010, double phi011,
												   double phi100, double phi101,
												   double phi110, double phi111, POINT3* Qs,
																			     POINT3* ns)
	{
		// perturbation 
		const double eps=1E-15;
		if(ABS(phi000)<eps) phi000 = (phi000>0)?eps:-eps;
		if(ABS(phi001)<eps) phi001 = (phi001>0)?eps:-eps;
		if(ABS(phi010)<eps) phi010 = (phi010>0)?eps:-eps;
		if(ABS(phi011)<eps) phi011 = (phi011>0)?eps:-eps;
		if(ABS(phi100)<eps) phi100 = (phi100>0)?eps:-eps;
		if(ABS(phi101)<eps) phi101 = (phi101>0)?eps:-eps;
		if(ABS(phi110)<eps) phi110 = (phi110>0)?eps:-eps;
		if(ABS(phi111)<eps) phi111 = (phi111>0)?eps:-eps;
		
		// 
		int count_pos=0;
		if(phi000>0) count_pos++;
		if(phi001>0) count_pos++;
		if(phi010>0) count_pos++;
		if(phi011>0) count_pos++;
		if(phi100>0) count_pos++;
		if(phi101>0) count_pos++;
		if(phi110>0) count_pos++;
		if(phi111>0) count_pos++;

		//
		if(count_pos==0 || count_pos==8) return 0;
		else
		{
			double dx=xip1-xi;
			double dy=yjp1-yj;
			double dz=zkp1-zk;

			POINT3 P000(xi  ,yj  ,zk  );
			POINT3 P001(xi  ,yj  ,zkp1);
			POINT3 P010(xi  ,yjp1,zk  );
			POINT3 P011(xi  ,yjp1,zkp1);
			POINT3 P100(xip1,yj  ,zk  );
			POINT3 P101(xip1,yj  ,zkp1);
			POINT3 P110(xip1,yjp1,zk  );
			POINT3 P111(xip1,yjp1,zkp1);
		
			const POINT3& P0 = P000; double phi0 = phi000; 
			const POINT3& P3 = P111; double phi3 = phi111; 
				  POINT3  P1       ; double phi1         ;
				  POINT3  P2       ; double phi2         ;
				  POINT3  N;
			
			// Kuhn Triangulation
			int number_of_trs=0;
			for(int n=0;n<6;n++)
			{
				//
					 if(n==0){ P1=P001; phi1=phi001; P2=P011; phi2=phi011; }
				else if(n==1){ P1=P001; phi1=phi001; P2=P101; phi2=phi101; }
				else if(n==2){ P1=P010; phi1=phi010; P2=P011; phi2=phi011; }
				else if(n==3){ P1=P010; phi1=phi010; P2=P110; phi2=phi110; }
				else if(n==4){ P1=P100; phi1=phi100; P2=P101; phi2=phi101; }
				else         { P1=P100; phi1=phi100; P2=P110; phi2=phi110; }

				//
				     if(n==0){ N=POINT3((phi3-phi2)/dx,(phi2-phi1)/dy,(phi1-phi0)/dz); }
				else if(n==1){ N=POINT3((phi2-phi1)/dx,(phi3-phi2)/dy,(phi1-phi0)/dz); }
				else if(n==2){ N=POINT3((phi3-phi2)/dx,(phi1-phi0)/dy,(phi2-phi1)/dz); }
				else if(n==3){ N=POINT3((phi2-phi1)/dx,(phi1-phi0)/dy,(phi3-phi2)/dz); }
				else if(n==4){ N=POINT3((phi1-phi0)/dx,(phi3-phi2)/dy,(phi2-phi1)/dz); }
				else         { N=POINT3((phi1-phi0)/dx,(phi2-phi1)/dy,(phi3-phi2)/dz); }

				//
				int type = ((phi0>0)?8:0)
					     + ((phi1>0)?4:0)
					     + ((phi2>0)?2:0)
					     + ((phi3>0)?1:0);

				if(type>=8) type=15-type;
			
				// there are 8 types
				if(type==0) ;
				else
				{
					*ns=N;ns++;
					     if(type==1)/*---+*/{ isosurfacing_mppp(P3,phi3,P0,phi0,P1,phi1,P2,phi2,Qs);number_of_trs++ ;}
					else if(type==2)/*--+-*/{ isosurfacing_mppp(P2,phi2,P0,phi0,P1,phi1,P3,phi3,Qs);number_of_trs++ ;}
					else if(type==3)/*--++*/{ isosurfacing_mmpp(P0,phi0,P1,phi1,P2,phi2,P3,phi3,Qs);number_of_trs+=2;*ns=N;ns++;}
					else if(type==4)/*-+--*/{ isosurfacing_mppp(P1,phi1,P0,phi0,P2,phi2,P3,phi3,Qs);number_of_trs++ ;}
					else if(type==5)/*-+-+*/{ isosurfacing_mmpp(P0,phi0,P2,phi2,P1,phi1,P3,phi3,Qs);number_of_trs+=2;*ns=N;ns++;}
					else if(type==6)/*-++-*/{ isosurfacing_mmpp(P1,phi1,P2,phi2,P0,phi0,P3,phi3,Qs);number_of_trs+=2;*ns=N;ns++;}
					else            /*-+++*/{ isosurfacing_mppp(P0,phi0,P1,phi1,P2,phi2,P3,phi3,Qs);number_of_trs++ ;}
				}
			}
			return number_of_trs;
		}
	}
	
	/*--------------------------------------------------------------------------------
	* isosurfacing on uniform grid
	* Vs : vertices of the triangulation
	* Ns : normals  of the triangulation
	* Thus size(Vs)=3*size(Ns)
	*------------------------------------------------------------------------------*/
	static void isosurface( const GRID3D& grid, 
	                        const ARRAY3D<double>& Phi, ARRAY<POINT3>& Vs_arr, 
							                            ARRAY<POINT3>& Ns_arr ) {
		int size=0, buff_size=1000;
		POINT3* Vs = (POINT3*)malloc(sizeof(POINT3)*buff_size*3);
		POINT3* Ns = (POINT3*)malloc(sizeof(POINT3)*buff_size  );
		POINT3 Vs_ijk[36];
		POINT3 Ns_ijk[12];
		
		int isize = grid.isize; double dx = grid.dx;
		int jsize = grid.jsize; double dy = grid.dy;
		int ksize = grid.ksize; double dz = grid.dz;
		for(int i=0;i<isize-1;i++){ double xi=grid.x_fr_i(i);
		for(int j=0;j<jsize-1;j++){ double yj=grid.y_fr_j(j);
		for(int k=0;k<ksize-1;k++){ double zk=grid.z_fr_k(k);
			int size_ijk=isosurface(xi,xi+dx,
									yj,yj+dy,
									zk,zk+dz,Phi[i  ][j  ][k],Phi[i  ][j  ][k+1],
									         Phi[i  ][j+1][k],Phi[i  ][j+1][k+1],
											 Phi[i+1][j  ][k],Phi[i+1][j  ][k+1],
									         Phi[i+1][j+1][k],Phi[i+1][j+1][k+1],Vs_ijk,Ns_ijk);
			if(size+size_ijk>=buff_size){
				POINT3* Vs_temp = (POINT3*)malloc(sizeof(POINT3)*buff_size*2*3); 
				POINT3* Ns_temp = (POINT3*)malloc(sizeof(POINT3)*buff_size*2  );
				memcpy(Vs_temp,Vs,sizeof(POINT3)*buff_size*3); free(Vs); Vs=Vs_temp;
				memcpy(Ns_temp,Ns,sizeof(POINT3)*buff_size  ); free(Ns); Ns=Ns_temp; buff_size*=2;
			}
			
			memcpy(Vs+3*size,Vs_ijk,sizeof(POINT3)*size_ijk*3);
			memcpy(Ns+  size,Ns_ijk,sizeof(POINT3)*size_ijk  ); size+=size_ijk;
		}}}
		
		Vs_arr.resize(size*3); memcpy((POINT3*)Vs_arr,Vs,sizeof(POINT3)*size*3); free(Vs);
		Ns_arr.resize(size  ); memcpy((POINT3*)Ns_arr,Ns,sizeof(POINT3)*size  ); free(Ns);
	}
};

#endif

