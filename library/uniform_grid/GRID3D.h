#ifndef GRID_3D__H
#define GRID_3D__H

#include "../ARRAY3D.h"
#include <omp.h>

struct GRID3D
{
	double xmin,xmax,dx; int isize;
	double ymin,ymax,dy; int jsize;
	double zmin,zmax,dz; int ksize;

	// basic grid operations
	void set_grid(double xmin_,double xmax_,int isize_,
				  double ymin_,double ymax_,int jsize_,
				  double zmin_,double zmax_,int ksize_) {
		xmin = xmin_;xmax = xmax_;isize= isize_;dx=(xmax-xmin)/(isize-1);
		ymin = ymin_;ymax = ymax_;jsize= jsize_;dy=(ymax-ymin)/(jsize-1);
		zmin = zmin_;zmax = zmax_;ksize= ksize_;dz=(zmax-zmin)/(ksize-1);}

	inline double x_fr_i(double i) const{ return xmin+dx*i;}
	inline double y_fr_j(double j) const{ return ymin+dy*j;}
	inline double z_fr_k(double k) const{ return zmin+dz*k;}
	inline double i_fr_x(double x) const{ return (x-xmin)/dx;}
	inline double j_fr_y(double y) const{ return (y-ymin)/dy;}
	inline double k_fr_z(double z) const{ return (z-zmin)/dz;}

	// Finite Differences
	inline double get_value_with_linear_ext( const ARRAY3D<double>& F, int i, int j, int k) const
	{
		int i0=MIN(MAX(i,0),F.isize-2); i-=i0;
		int j0=MIN(MAX(j,0),F.jsize-2); j-=j0;
		int k0=MIN(MAX(k,0),F.ksize-2); k-=k0;
		if(i==0 && j==0 && k==0 ) return F[i0][j0][k0];
		else	return F[i0  ][j0  ][k0  ]*(1-i)*(1-j)*(1-k)+
                       F[i0  ][j0+1][k0  ]*(1-i)*(  j)*(1-k)+
					   F[i0+1][j0  ][k0  ]*(  i)*(1-j)*(1-k)+
					   F[i0+1][j0+1][k0  ]*(  i)*(  j)*(1-k)+
					   F[i0  ][j0  ][k0+1]*(1-i)*(1-j)*(  k)+
                       F[i0  ][j0+1][k0+1]*(1-i)*(  j)*(  k)+
					   F[i0+1][j0  ][k0+1]*(  i)*(1-j)*(  k)+
					   F[i0+1][j0+1][k0+1]*(  i)*(  j)*(  k);
	}

	void Dx_2nd( const ARRAY3D<double>& F, int i, int j, int k, double& fx_m, double& fx_p ) const
	{
		double fim2 = get_value_with_linear_ext(F,i-2,j,k);
		double fim1 = get_value_with_linear_ext(F,i-1,j,k); double Dfim2 = fim1-fim2;
		double fi   = get_value_with_linear_ext(F,i  ,j,k); double Dfim1 = fi  -fim1; double DDfim2 = .5*(Dfim1-Dfim2);
		double fip1 = get_value_with_linear_ext(F,i+1,j,k); double Dfi   = fip1-fi  ; double DDfim1 = .5*(Dfi  -Dfim1);
		double fip2 = get_value_with_linear_ext(F,i+2,j,k); double Dfip1 = fip2-fip1; double DDfi   = .5*(Dfip1-Dfi  );

		fx_m = (Dfim1+MINMOD(DDfim2,DDfim1))/dx;
		fx_p = (Dfi  -MINMOD(DDfim1,DDfi  ))/dx;
	}

	void Dy_2nd( const ARRAY3D<double>& F, int i, int j, int k, double& fy_m, double& fy_p ) const
	{
		double fjm2 = get_value_with_linear_ext(F,i,j-2,k);
		double fjm1 = get_value_with_linear_ext(F,i,j-1,k); double Dfjm2 = fjm1-fjm2;
		double fj   = get_value_with_linear_ext(F,i,j  ,k); double Dfjm1 = fj  -fjm1; double DDfjm2 = .5*(Dfjm1-Dfjm2);
		double fjp1 = get_value_with_linear_ext(F,i,j+1,k); double Dfj   = fjp1-fj  ; double DDfjm1 = .5*(Dfj  -Dfjm1);
		double fjp2 = get_value_with_linear_ext(F,i,j+2,k); double Dfjp1 = fjp2-fjp1; double DDfj   = .5*(Dfjp1-Dfj  );

		fy_m = (Dfjm1+MINMOD(DDfjm2,DDfjm1))/dy;
		fy_p = (Dfj  -MINMOD(DDfjm1,DDfj  ))/dy;
	}

	void Dz_2nd( const ARRAY3D<double>& F, int i, int j, int k, double& fz_m, double& fz_p ) const
	{
		double fkm2 = get_value_with_linear_ext(F,i,j,k-2);
		double fkm1 = get_value_with_linear_ext(F,i,j,k-1); double Dfkm2 = fkm1-fkm2;
		double fk   = get_value_with_linear_ext(F,i,j,k  ); double Dfkm1 = fk  -fkm1; double DDfkm2 = .5*(Dfkm1-Dfkm2);
		double fkp1 = get_value_with_linear_ext(F,i,j,k+1); double Dfk   = fkp1-fk  ; double DDfkm1 = .5*(Dfk  -Dfkm1);
		double fkp2 = get_value_with_linear_ext(F,i,j,k+2); double Dfkp1 = fkp2-fkp1; double DDfk   = .5*(Dfkp1-Dfk  );

		fz_m = (Dfkm1+MINMOD(DDfkm2,DDfkm1))/dz;
		fz_p = (Dfk  -MINMOD(DDfkm1,DDfk  ))/dz;
	}

	inline void Dx_2nd_subcell( const ARRAY3D<double>& Phi,
								const ARRAY3D<double>& F  , 
								double f_interface, int i, int j, int k,
								double& dx_m, double& dx_p, 
								double& fx_m, double& fx_p) const {
		double fm2 = get_value_with_linear_ext(F  ,i-2,j,k); 
		double fm1 = get_value_with_linear_ext(F  ,i-1,j,k), Dfm2 = fm1-fm2;
		double f   = get_value_with_linear_ext(F  ,i  ,j,k), Dfm1 = f  -fm1, DDfm2 = .5*(Dfm1-Dfm2);
		double fp1 = get_value_with_linear_ext(F  ,i+1,j,k), Df   = fp1-f  , DDfm1 = .5*(Df  -Dfm1);
		double fp2 = get_value_with_linear_ext(F  ,i+2,j,k), Dfp1 = fp2-fp1, DDf   = .5*(Dfp1-Df  );
		double pm2 = get_value_with_linear_ext(Phi,i-2,j,k);
		double pm1 = get_value_with_linear_ext(Phi,i-1,j,k), Dpm2 = pm1-pm2;
		double p   = get_value_with_linear_ext(Phi,i  ,j,k), Dpm1 = p  -pm1, DDpm2 = .5*(Dpm1-Dpm2);
		double pp1 = get_value_with_linear_ext(Phi,i+1,j,k), Dp   = pp1-p  , DDpm1 = .5*(Dp  -Dpm1);
		double pp2 = get_value_with_linear_ext(Phi,i+2,j,k), Dpp1 = pp2-pp1, DDp   = .5*(Dpp1-Dp  );
		dx_m = (p*pm1<0) ? dx*subcell_resolution(p,pm1,2*MINMOD(DDpm2,DDpm1)) : dx;
		dx_p = (p*pp1<0) ? dx*subcell_resolution(p,pp1,2*MINMOD(DDp  ,DDpm1)) : dx;
		fx_m = ((p*pm1<0)?(f-f_interface):(f-fm1))/dx_m+dx_m*MINMOD(DDfm2,DDfm1)/dx/dx; 
		fx_p = ((p*pp1<0)?(f_interface-f):(fp1-f))/dx_p-dx_p*MINMOD(DDf  ,DDfm1)/dx/dx; }

	inline void Dy_2nd_subcell( const ARRAY3D<double>& Phi,
								const ARRAY3D<double>& F  , 
								double f_interface, int i, int j, int k,
								double& dy_m, double& dy_p, 
								double& fy_m, double& fy_p) const {
		double fm2 = get_value_with_linear_ext(F  ,i,j-2,k); 
		double fm1 = get_value_with_linear_ext(F  ,i,j-1,k), Dfm2 = fm1-fm2;
		double f   = get_value_with_linear_ext(F  ,i,j  ,k), Dfm1 = f  -fm1, DDfm2 = .5*(Dfm1-Dfm2);
		double fp1 = get_value_with_linear_ext(F  ,i,j+1,k), Df   = fp1-f  , DDfm1 = .5*(Df  -Dfm1);
		double fp2 = get_value_with_linear_ext(F  ,i,j+2,k), Dfp1 = fp2-fp1, DDf   = .5*(Dfp1-Df  );
		double pm2 = get_value_with_linear_ext(Phi,i,j-2,k);
		double pm1 = get_value_with_linear_ext(Phi,i,j-1,k), Dpm2 = pm1-pm2;
		double p   = get_value_with_linear_ext(Phi,i,j  ,k), Dpm1 = p  -pm1, DDpm2 = .5*(Dpm1-Dpm2);
		double pp1 = get_value_with_linear_ext(Phi,i,j+1,k), Dp   = pp1-p  , DDpm1 = .5*(Dp  -Dpm1);
		double pp2 = get_value_with_linear_ext(Phi,i,j+2,k), Dpp1 = pp2-pp1, DDp   = .5*(Dpp1-Dp  );
		dy_m = (p*pm1<0) ? dy*subcell_resolution(p,pm1,2*MINMOD(DDpm2,DDpm1)) : dy;
		dy_p = (p*pp1<0) ? dy*subcell_resolution(p,pp1,2*MINMOD(DDp  ,DDpm1)) : dy;
		fy_m = ((p*pm1<0)?(f-f_interface):(f-fm1))/dy_m+dy_m*MINMOD(DDfm2,DDfm1)/dy/dy; 
		fy_p = ((p*pp1<0)?(f_interface-f):(fp1-f))/dy_p-dy_p*MINMOD(DDf  ,DDfm1)/dy/dy; }

	inline void Dz_2nd_subcell( const ARRAY3D<double>& Phi,
								const ARRAY3D<double>& F  , 
								double f_interface, int i, int j, int k,
								double& dz_m, double& dz_p, 
								double& fz_m, double& fz_p) const {
		double fm2 = get_value_with_linear_ext(F  ,i,j,k-2); 
		double fm1 = get_value_with_linear_ext(F  ,i,j,k-1), Dfm2 = fm1-fm2;
		double f   = get_value_with_linear_ext(F  ,i,j,k  ), Dfm1 = f  -fm1, DDfm2 = .5*(Dfm1-Dfm2);
		double fp1 = get_value_with_linear_ext(F  ,i,j,k+1), Df   = fp1-f  , DDfm1 = .5*(Df  -Dfm1);
		double fp2 = get_value_with_linear_ext(F  ,i,j,k+2), Dfp1 = fp2-fp1, DDf   = .5*(Dfp1-Df  );
		double pm2 = get_value_with_linear_ext(Phi,i,j,k-2);
		double pm1 = get_value_with_linear_ext(Phi,i,j,k-1), Dpm2 = pm1-pm2;
		double p   = get_value_with_linear_ext(Phi,i,j,k  ), Dpm1 = p  -pm1, DDpm2 = .5*(Dpm1-Dpm2);
		double pp1 = get_value_with_linear_ext(Phi,i,j,k+1), Dp   = pp1-p  , DDpm1 = .5*(Dp  -Dpm1);
		double pp2 = get_value_with_linear_ext(Phi,i,j,k+2), Dpp1 = pp2-pp1, DDp   = .5*(Dpp1-Dp  );
		dz_m = (p*pm1<0) ? dz*subcell_resolution(p,pm1,2*MINMOD(DDpm2,DDpm1)) : dz;
		dz_p = (p*pp1<0) ? dz*subcell_resolution(p,pp1,2*MINMOD(DDp  ,DDpm1)) : dz;
		fz_m = ((p*pm1<0)?(f-f_interface):(f-fm1))/dz_m+dz_m*MINMOD(DDfm2,DDfm1)/dz/dz; 
		fz_p = ((p*pp1<0)?(f_interface-f):(fp1-f))/dz_p-dz_p*MINMOD(DDf  ,DDfm1)/dz/dz; }
	//---------------------------------------------------------------------
	// Extrapolation of F from {F:defined} to {F:undefined} by the
	// thin-plate spline RBF interpolation
	//---------------------------------------------------------------------
	static void extrapolate_RBF( const ARRAY3D<bool  >& F_defined,
		                               ARRAY3D<double>& F  , int band_width )
	{
		int isize=F.isize;
		int jsize=F.jsize;
		int ksize=F.ksize;
		// calculate the nodes to extrapolate
		ARRAY<int> i_nodes(isize); i_nodes=0;
		for(int i=0;i<isize;i++)
		for(int j=0;j<jsize;j++)
		for(int k=0;k<ksize;k++)
			if(F_defined[i][j][k]==false)
				i_nodes[i]++;
		// parallel
		int NT;
		int ibeg_T[16];
		int iend_T[16];
		#pragma omp parallel shared(i_nodes,ibeg_T,iend_T,NT) firstprivate(isize,jsize,ksize,band_width)
		{
			#pragma omp single
			{ 
				NT=omp_get_num_threads();
				load_distribution(NT,isize,i_nodes,ibeg_T,iend_T);
			}
			int T=omp_get_thread_num();

			RBF_3D rbf;
			ARRAY<double> xs(sizeof(double)*CUBE(2*band_width+1));
			ARRAY<double> ys(sizeof(double)*CUBE(2*band_width+1));
			ARRAY<double> zs(sizeof(double)*CUBE(2*band_width+1));
			ARRAY<double> fs(sizeof(double)*CUBE(2*band_width+1)); int N;
			
			for(int i=ibeg_T[T];i<=iend_T[T];i++)
			for(int j=0;j<jsize;j++)
			for(int k=0;k<ksize;k++)
				if(F_defined[i][j][k]==false)
				{
					int a_m = MAX(0,i-band_width); int a_p=MIN(isize-1,i+band_width);
					int b_m = MAX(0,j-band_width); int b_p=MIN(jsize-1,j+band_width);
					int c_m = MAX(0,k-band_width); int c_p=MIN(ksize-1,k+band_width); N=0;

					for(int a=a_m;a<=a_p;a++)
					for(int b=b_m;b<=b_p;b++)
					for(int c=c_m;c<=c_p;c++)
					{
						if( SQR(i-a)+SQR(j-b)+SQR(k-c)<=SQR(band_width) && F_defined[a][b][c])
						{
							xs[N]=a;
							ys[N]=b;
							zs[N]=c;
							fs[N]=F[a][b][c]; N++;
						}
					}

					if(N>=3)
						F[i][j][k]=rbf.interpolate((double*)xs,
												   (double*)ys,
												   (double*)zs,
												   (double*)fs,i,j,k,N);
			}
		}
	}

	//---------------------------------------------------------------------
	// Extrapolation of F from {F:defined} to {F:undefined} by the
	// thin-plate spline RBF interpolation
	//---------------------------------------------------------------------
	void extrapolate_Shepard( const ARRAY3D<bool  >& F_defined,
		                            ARRAY3D<double>& F  , int band_width )
	{
		for(int i=0;i<isize;i++)
		for(int j=0;j<jsize;j++)
		for(int k=0;k<ksize;k++)
			if(F_defined[i][j][k]==false)
			{
				double fw_sum=0;
				double  w_sum=0;
				int a_m = MAX(0,i-band_width); int a_p=MIN(isize-1,i+band_width);
				int b_m = MAX(0,j-band_width); int b_p=MIN(jsize-1,j+band_width); 
				int c_m = MAX(0,k-band_width); int c_p=MAX(ksize-1,k+band_width);
				for(int a=a_m;a<=a_p;a++)
				for(int b=b_m;b<=b_p;b++)
				for(int c=c_m;c<=c_p;c++)
				{
					if( SQR(i-a)+SQR(j-b)+SQR(k-c)<=SQR(band_width) && F_defined[a][b][c]==true)
					{
						double w = 1/(SQR(dx*(i-a))+SQR(dy*(j-b))+SQR(dz*(k-c)));
						 w_sum +=            w;
						fw_sum += F[a][b][c]*w;
					}
				}

				if(w_sum!=0)
					F[i][j][k] = fw_sum/w_sum;
			}
	}
	//---------------------------------------------------------------------
	// ENO-2 interpolation
	//---------------------------------------------------------------------
	double ENO2_interpolation( const ARRAY3D<double>& F, double x, double y, double z ) const
	{
		double i=i_fr_x(x); int i0=int(i)-1; i0=MIN(MAX(i0,0),isize-4);
		double j=j_fr_y(y); int j0=int(j)-1; j0=MIN(MAX(j0,0),jsize-4);
		double k=k_fr_z(z); int k0=int(k)-1; k0=MIN(MAX(k0,0),ksize-4);
		double f[4][4][4];
		for(int a=0;a<4;a++)
		for(int b=0;b<4;b++)
		for(int c=0;c<4;c++) f[a][b][c]=F[i0+a][j0+b][k0+c];
		return ::ENO2_interpolation(f,i-i0,j-j0,k-k0);
	}
};

#endif
