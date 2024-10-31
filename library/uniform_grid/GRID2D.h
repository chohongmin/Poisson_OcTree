#ifndef GRID_2D__H
#define GRID_2D__H

#include "../geometry/RBF.h"

//
struct GRID2D {
	double xmin,xmax,dx; int isize;
	double ymin,ymax,dy; int jsize;

	// basic grid operations
	void set_grid(double xmin_,double xmax_,int isize_,
				  double ymin_,double ymax_,int jsize_) {
		xmin = xmin_;xmax = xmax_;isize= isize_;dx=(xmax-xmin)/(isize-1);
		ymin = ymin_;ymax = ymax_;jsize= jsize_;dy=(ymax-ymin)/(jsize-1); }

	inline double x_fr_i(double i) const{ return xmin+dx*i;}
	inline double y_fr_j(double j) const{ return ymin+dy*j;}
	inline double i_fr_x(double x) const{ return (x-xmin)/dx;}
	inline double j_fr_y(double y) const{ return (y-ymin)/dy;}

	// Finite Differences
	inline double get_value_with_linear_ext( const ARRAY2D<double>& F, int i, int j) const
	{
		int i0=MIN(MAX(i,0),F.isize-2);
		int j0=MIN(MAX(j,0),F.jsize-2);
		if(i==i0 && j==j0) return F[i][j];
		else			   return F[i0  ][j0  ]*(1-(i-i0))*(1-(j-j0))+
			                      F[i0  ][j0+1]*(1-(i-i0))*(  (j-j0))+
								  F[i0+1][j0  ]*(  (i-i0))*(1-(j-j0))+
								  F[i0+1][j0+1]*(  (i-i0))*(  (j-j0));
	}
	void Dx_2nd( const ARRAY2D<double>& F, int i, int j, double& fx_m, double& fx_p ) const {
		double fim2 = get_value_with_linear_ext(F,i-2,j);
		double fim1 = get_value_with_linear_ext(F,i-1,j); double Dfim2 = fim1-fim2;
		double fi   = get_value_with_linear_ext(F,i  ,j); double Dfim1 = fi  -fim1; double DDfim2 = .5*(Dfim1-Dfim2);
		double fip1 = get_value_with_linear_ext(F,i+1,j); double Dfi   = fip1-fi  ; double DDfim1 = .5*(Dfi  -Dfim1);
		double fip2 = get_value_with_linear_ext(F,i+2,j); double Dfip1 = fip2-fip1; double DDfi   = .5*(Dfip1-Dfi  );
		fx_m = (Dfim1+MINMOD(DDfim2,DDfim1))/dx;
		fx_p = (Dfi  -MINMOD(DDfim1,DDfi  ))/dx; }

	void Dy_2nd( const ARRAY2D<double>& F, int i, int j, double& fy_m, double& fy_p ) const {
		double fjm2 = get_value_with_linear_ext(F,i,j-2);
		double fjm1 = get_value_with_linear_ext(F,i,j-1); double Dfjm2 = fjm1-fjm2;
		double fj   = get_value_with_linear_ext(F,i,j  ); double Dfjm1 = fj  -fjm1; double DDfjm2 = .5*(Dfjm1-Dfjm2);
		double fjp1 = get_value_with_linear_ext(F,i,j+1); double Dfj   = fjp1-fj  ; double DDfjm1 = .5*(Dfj  -Dfjm1);
		double fjp2 = get_value_with_linear_ext(F,i,j+2); double Dfjp1 = fjp2-fjp1; double DDfj   = .5*(Dfjp1-Dfj  );
		fy_m = (Dfjm1+MINMOD(DDfjm2,DDfjm1))/dy;
		fy_p = (Dfj  -MINMOD(DDfjm1,DDfj  ))/dy; }

	inline void Dx_2nd_subcell( const ARRAY2D<double>& Phi,
								const ARRAY2D<double>& F  , 
								double f_interface, int i, int j,
								double& dx_m, double& dx_p, 
								double& fx_m, double& fx_p) const {
		double fim2 = get_value_with_linear_ext(F,i-2,j); 
		double fim1 = get_value_with_linear_ext(F,i-1,j), Dfim2 = fim1-fim2;
		double fi   = get_value_with_linear_ext(F,i  ,j), Dfim1 = fi  -fim1, DDfim2 = .5*(Dfim1-Dfim2);
		double fip1 = get_value_with_linear_ext(F,i+1,j), Dfi   = fip1-fi  , DDfim1 = .5*(Dfi  -Dfim1);
		double fip2 = get_value_with_linear_ext(F,i+2,j), Dfip1 = fip2-fip1, DDfi   = .5*(Dfip1-Dfi  );
		double pim2 = get_value_with_linear_ext(Phi,i-2,j);
		double pim1 = get_value_with_linear_ext(Phi,i-1,j), Dpim2 = pim1-pim2;
		double pi   = get_value_with_linear_ext(Phi,i  ,j), Dpim1 = pi  -pim1, DDpim2 = .5*(Dpim1-Dpim2);
		double pip1 = get_value_with_linear_ext(Phi,i+1,j), Dpi   = pip1-pi  , DDpim1 = .5*(Dpi  -Dpim1);
		double pip2 = get_value_with_linear_ext(Phi,i+2,j), Dpip1 = pip2-pip1, DDpi   = .5*(Dpip1-Dpi  );
		dx_m = (pi*pim1<0) ? dx*subcell_resolution(pi,pim1,2*MINMOD(DDpim2,DDpim1)) : dx;
		dx_p = (pi*pip1<0) ? dx*subcell_resolution(pi,pip1,2*MINMOD(DDpi  ,DDpim1)) : dx;
		fx_m = ((pi*pim1<0)?(fi-f_interface):(fi-fim1))/dx_m+dx_m*MINMOD(DDfim2,DDfim1)/dx/dx; 
		fx_p = ((pi*pip1<0)?(f_interface-fi):(fip1-fi))/dx_p-dx_p*MINMOD(DDfi  ,DDfim1)/dx/dx; }
	
	//
	inline void Dy_2nd_subcell( const ARRAY2D<double>& Phi,
								const ARRAY2D<double>& F  , 
								double f_interface, int i, int j,
								double& dy_m, double& dy_p, 
								double& fy_m, double& fy_p) const {
		double fjm2 = get_value_with_linear_ext(F  ,i,j-2); 
		double fjm1 = get_value_with_linear_ext(F  ,i,j-1), Dfjm2 = fjm1-fjm2;
		double fj   = get_value_with_linear_ext(F  ,i,j  ), Dfjm1 = fj  -fjm1, DDfjm2 = .5*(Dfjm1-Dfjm2);
		double fjp1 = get_value_with_linear_ext(F  ,i,j+1), Dfj   = fjp1-fj  , DDfjm1 = .5*(Dfj  -Dfjm1);
		double fjp2 = get_value_with_linear_ext(F  ,i,j+2), Dfjp1 = fjp2-fjp1, DDfj   = .5*(Dfjp1-Dfj  );
		double pjm2 = get_value_with_linear_ext(Phi,i,j-2);
		double pjm1 = get_value_with_linear_ext(Phi,i,j-1), Dpjm2 = pjm1-pjm2;
		double pj   = get_value_with_linear_ext(Phi,i,j  ), Dpjm1 = pj  -pjm1, DDpjm2 = .5*(Dpjm1-Dpjm2);
		double pjp1 = get_value_with_linear_ext(Phi,i,j+1), Dpj   = pjp1-pj  , DDpjm1 = .5*(Dpj  -Dpjm1);
		double pjp2 = get_value_with_linear_ext(Phi,i,j+2), Dpjp1 = pjp2-pjp1, DDpj   = .5*(Dpjp1-Dpj  );
		dy_m = (pj*pjm1<0) ? dy*subcell_resolution(pj,pjm1,2*MINMOD(DDpjm2,DDpjm1)) : dy;
		dy_p = (pj*pjp1<0) ? dy*subcell_resolution(pj,pjp1,2*MINMOD(DDpj  ,DDpjm1)) : dy;
		fy_m = ((pj*pjm1<0)?(fj-f_interface):(fj-fjm1))/dy_m+dy_m*MINMOD(DDfjm2,DDfjm1)/dy/dy; 
		fy_p = ((pj*pjp1<0)?(f_interface-fj):(fjp1-fj))/dy_p-dy_p*MINMOD(DDfj  ,DDfjm1)/dy/dy; }

	//---------------------------------------------------------------------
	// reinitialization
	//---------------------------------------------------------------------
	void reinitialize_GaussSeidal(const ARRAY2D<double>& P0,
								        ARRAY2D<double>& P , int it){
		bool iflip = (it%4)>=2 ? true : false;
		bool jflip = (it%2)>=1 ? true : false;
	#pragma omp parallel for 
		for(int i_=0;i_<P.isize;i_++)
		for(int j_=0;j_<P.jsize;j_++) {
			int i= iflip ? P.isize-1-i_ : i_;
			int j= jflip ? P.jsize-1-j_ : j_;

			double dx_m,dx_p, fx_m,fx_p; Dx_2nd_subcell(P0,P,0,i,j,dx_m,dx_p,fx_m,fx_p);
			double dy_m,dy_p, fy_m,fy_p; Dy_2nd_subcell(P0,P,0,i,j,dy_m,dy_p,fy_m,fy_p);
			// Godunov Hamiltonian
			double phi0 = P0[i][j];
			if(phi0>0) {
				if(fx_m<0) fx_m=0; if(fx_p>0) fx_p=0;
				if(fy_m<0) fy_m=0; if(fy_p>0) fy_p=0; }
			else {
				if(fx_m>0) fx_m=0; if(fx_p<0) fx_p=0;
				if(fy_m>0) fy_m=0; if(fy_p<0) fy_p=0; }
			//
			const double EPS=1E-10;
			if( dx_m<EPS || dx_p<EPS ||
				dy_m<EPS || dy_p<EPS) P[i][j]= 0;
			else {
				double H = SQRT(MAX(fx_m*fx_m, fx_p*fx_p)
							   +MAX(fy_m*fy_m, fy_p*fy_p));
	    		double dt = .45 *MIN(MIN(dx_m,dx_p),MIN(dy_m,dy_p));
				double sgn = (phi0>0)?1:-1;
				double pnp1 = P[i][j] - dt*sgn*(H-1);
				if(phi0*pnp1>0) P[i][j] = pnp1;
				else			P[i][j] =-pnp1; }}}
	//
	void reinitialize_Euler(const ARRAY2D<double>& P0,
							const ARRAY2D<double>& Pn,
								  ARRAY2D<double>& Pnp1){
	#pragma omp parallel for 
		for(int i=0;i<isize;i++)
		for(int j=0;j<jsize;j++) {
			double phi0 = P0[i][j];
			double phi  = Pn[i][j];
			double dx_m,dx_p, fx_m,fx_p; Dx_2nd_subcell(P0,Pn,0,i,j,dx_m,dx_p,fx_m,fx_p);
			double dy_m,dy_p, fy_m,fy_p; Dy_2nd_subcell(P0,Pn,0,i,j,dy_m,dy_p,fy_m,fy_p);
			// Godunov Hamiltonian
			if(phi0>0) {
				if(fx_m<0) fx_m=0; if(fx_p>0) fx_p=0;
				if(fy_m<0) fy_m=0; if(fy_p>0) fy_p=0; }
			else {
				if(fx_m>0) fx_m=0; if(fx_p<0) fx_p=0;
				if(fy_m>0) fy_m=0; if(fy_p<0) fy_p=0; }
			//
			const double EPS=1E-10;
			if( dx_m<EPS || dx_p<EPS ||
				dy_m<EPS || dy_p<EPS) Pnp1[i][j]= 0;
			else {
				double H = SQRT(MAX(fx_m*fx_m, fx_p*fx_p)
							   +MAX(fy_m*fy_m, fy_p*fy_p));
	    		double dt = .45 *MIN(MIN(dx_m,dx_p),MIN(dy_m,dy_p));
				double sgn = (phi0>0)?1:-1;
				double pnp1 = Pn[i][j] - dt*sgn*(H-1);
				if(phi0*pnp1>0) Pnp1[i][j] = pnp1;
				else			Pnp1[i][j] =-pnp1; }}}

	void reinitialize_RK2( const ARRAY2D<double>& P0,
		                   const ARRAY2D<double>& Pn,
						         ARRAY2D<double>& Pnp1){
		ARRAY2D<double> Pnp2(P0.isize,P0.jsize);
		reinitialize_Euler(P0,Pn  ,Pnp1);
		reinitialize_Euler(P0,Pnp1,Pnp2); 
		for(int i=0;i<P0.isize;i++)
		for(int j=0;j<P0.jsize;j++)
			Pnp1[i][j] = .5*(Pn[i][j]+Pnp2[i][j]);}

	void reinitialize_RK2( ARRAY2D<double>& P, int iteration ){
		ARRAY2D<double> P0=P;
		ARRAY2D<double> P1(P.isize,P.jsize);
		ARRAY2D<double> P2(P.isize,P.jsize);
		for(int n=0;n<iteration;n++) {
			reinitialize_Euler(P0,P ,P1); 
			reinitialize_Euler(P0,P1,P2); 
			for(int i=0;i<P.isize;i++)
			for(int j=0;j<P.jsize;j++)
				P[i][j] = .5*(P[i][j]+P2[i][j]);}}

	void reinitialize_GS( ARRAY2D<double>& P, int iteration ){
		ARRAY2D<double> P0=P;
		for(int n=0;n<iteration;n++) 
			reinitialize_GaussSeidal(P0,P,n);}
	//---------------------------------------------------------------------
	// Extrapolation of F from {F:defined} to {F:undefined} by the
	// thin-plate spline RBF interpolation
	//---------------------------------------------------------------------
	void extrapolate_RBF( const ARRAY2D<bool  >& F_defined,
		                        ARRAY2D<double>& F  , int band_width )
	{
		RBF_2D rbf;
		double* xs = (double*)malloc(sizeof(double)*SQR(2*band_width+1));
		double* ys = (double*)malloc(sizeof(double)*SQR(2*band_width+1));
		double* fs = (double*)malloc(sizeof(double)*SQR(2*band_width+1)); int N;
		int isize_=F.isize;
		int jsize_=F.jsize;
		for(int i=0;i<isize_;i++)
		for(int j=0;j<jsize_;j++)
			if(F_defined[i][j]==false)
			{
				int a_m = MAX(0,i-band_width); int a_p=MIN(isize_-1,i+band_width);
				int b_m = MAX(0,j-band_width); int b_p=MIN(jsize_-1,j+band_width); N=0;
				for(int a=a_m;a<=a_p;a++)
				for(int b=b_m;b<=b_p;b++)
				{
					if( SQR(i-a)+SQR(j-b)<=SQR(band_width) && F_defined[a][b]==true)
					{
						xs[N]=a;
						ys[N]=b;
						fs[N]=F[a][b]; N++;
					}
				}

				if(N>=3)
				{
					rbf.set_interpolation(xs,ys,fs,N);
					F[i][j]=rbf.evaluate(i,j);
				}
		}
		free(xs);
		free(ys);
		free(fs);
	}
	//---------------------------------------------------------------------
	// Extrapolation of F from {F:defined} to {F:undefined} by the
	// thin-plate spline RBF interpolation
	//---------------------------------------------------------------------
	void extrapolate_Shepard( const ARRAY2D<bool  >& F_defined,
		                            ARRAY2D<double>& F  , int band_width )
	{
		for(int i=0;i<isize;i++)
		for(int j=0;j<jsize;j++)
			if(F_defined[i][j]==false)
			{
				double fw_sum=0;
				double  w_sum=0;
				int a_m = MAX(0,i-band_width); int a_p=MIN(isize-1,i+band_width);
				int b_m = MAX(0,j-band_width); int b_p=MIN(jsize-1,j+band_width); 
				for(int a=a_m;a<=a_p;a++)
				for(int b=b_m;b<=b_p;b++)
				{
					if( SQR(i-a)+SQR(j-b)<=SQR(band_width) && F_defined[a][b]==true)
					{
						double w = 1/(SQR(dx*(i-a))+SQR(dy*(j-b)));
						 w_sum +=        w;
						fw_sum += F[a][b]*w;
					}
				}

				if(w_sum!=0)
					F[i][j] = fw_sum/w_sum;
			}
	}
	//---------------------------------------------------------------------
	// ENO-2 interpolation
	//---------------------------------------------------------------------
	double ENO2_interpolation( const ARRAY2D<double>& F, double x, double y ) const
	{
		double i=i_fr_x(x); int i0=int(i)-1; i0=MIN(MAX(i0,0),isize-4);
		double j=j_fr_y(y); int j0=int(j)-1; j0=MIN(MAX(j0,0),jsize-4);
		double f[4][4];
		for(int a=0;a<4;a++)
		for(int b=0;b<4;b++) f[a][b]=F[i0+a][j0+b];
		return ::ENO2_interpolation(f,i-i0,j-j0);
	}
};

#endif
