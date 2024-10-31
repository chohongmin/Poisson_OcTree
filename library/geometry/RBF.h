#ifndef EXTRAPOLATION_RBF__H
#define EXTRAPOLATION_RBF__H

#include "../algebra/MATRIX.h"

/*-------------------------------------------------------------------------------------------*\
 * The thin-plate spline interpolation in 2D
 * f(x) = sum_i ai*phi(|x-xi|) + sum_j bj*p_j(x)
\*-------------------------------------------------------------------------------------------*/
class RBF_2D 
{
protected:
	int N;
	double *xs, *ys;         // data points
	double *as;              // coeffs ai
	double  b0,b1,b2, x0,y0; // coeff bj expanded at x0,y0

	double get_phi(double r) const{ return (r>1E-7)?r*r*log(r):0;}
	double inner_product( const double* u, 
						  const double* v ) const 
	{ 
		double sum=0; for(int i=N;i>0;i--,u++,v++) sum+=(*u)*(*v); return sum; 
	}
	double inner_product_with_one( const double* u ) const 
	{
		double sum=0; for(int i=N;i>0;i--,u++    ) sum+=(*u)     ; return sum; 
	}
	// K=(B^TB)^(-1)
	double k00,k01,k02; 
	double k10,k11,k12;
	double k20,k21,k22;
	// P = I - BKB^T, B=(X,Y,1)	
	void project(double* v) {
		double bv0 = inner_product(xs,v); 
		double bv1 = inner_product(ys,v); 
		double bv2 = inner_product_with_one(v); 
		double kbv0 = k00*bv0 + k01*bv1 + k02*bv2;
		double kbv1 = k10*bv0 + k11*bv1 + k12*bv2;
		double kbv2 = k20*bv0 + k21*bv1 + k22*bv2;
		for(int i=0;i<N;i++) v[i] -= xs[i]*kbv0 + ys[i]*kbv1 + kbv2; }
	
public:
	RBF_2D(){N=0;xs=ys=as=0;}
	~RBF_2D(){if(N!=0){free(xs);free(ys);free(as);}}

	/*-------------------------------------------------------------------------------------------*\
	 * the expansion point (x0,y0) is the average of the points (xi,yi).
	 * calculate the coefficients ai and bj from the linear system
	 * | A   B | |a| = |f|
	 * | B^T 0 | |b|   |0|,   Aik=phi(|xi-xk|), Bij=p_j(xi-x0)
	 * The linear system is solved by CG on PAPc=Pf, where P=I-B(B^TB)^(-1)B^T
	 * then set a=Pc and b=(B^TB)^(-1)B^T(f-Aa)
	\*-------------------------------------------------------------------------------------------*/
	void set_interpolation( const double* xs_,
		                    const double* ys_,
							const double* fs , int N_)
	{
		//
		if(N<N_){free(xs);free(ys);free(as);
		         xs=(double*)malloc(sizeof(double)*N_);
				 ys=(double*)malloc(sizeof(double)*N_);
				 as=(double*)malloc(sizeof(double)*N_);}
		N=N_;
		memcpy(xs,xs_,sizeof(double)*N);
		memcpy(ys,ys_,sizeof(double)*N);
		// find the average (x0,y0)
		// p0(x,y)=x-x0, p1(x,y)=y-y0, p2(x,y)=1
		x0=y0=0;
		for(int i=0;i<N;i++){x0+=xs[i];y0+=ys[i];} x0/=N; y0/=N;
		for(int i=0;i<N;i++){xs[i]-=x0;ys[i]-=y0;}
		// calculate K=(B^TB)^(-1)
		double bb00 = inner_product(xs,xs);
		double bb11 = inner_product(ys,ys); 
		double bb22 = N;
		double bb01 = inner_product      (xs,ys); double bb10=bb01;
		double bb02 = inner_product_with_one(xs); double bb20=bb02;
		double bb12 = inner_product_with_one(ys); double bb21=bb12;

		double det = bb00*(bb11*bb22-bb12*bb21)
			        -bb01*(bb10*bb22-bb12*bb20)
					+bb02*(bb10*bb21-bb11*bb20);
		// in the case when the points are collinear, constant function
		if(ABS(det)<1E-12) 
		{ 
			b0=b1=b2=0;
			for(int i=0;i<N;i++){ as[i]=0; b2+=fs[i]; } b2/=N;
		}
		else
		{
			// calculate the matrix A
			// B=[xs,ys,1]
			double** A=(double**)malloc(sizeof(double*)*N);
			for(int i=0;i<N;i++) A[i]=(double*)malloc(sizeof(double)*N);
			for(int i=0;i<N;i++)
			for(int k=i;k<N;k++){
				double rik = SQRT(SQR(xs[i]-xs[k])+SQR(ys[i]-ys[k]));
				A[i][k]=A[k][i]=get_phi(rik);
			}
			
			k00=(bb11*bb22-bb12*bb21)/det;
			k11=(bb00*bb22-bb02*bb20)/det;
			k22=(bb00*bb11-bb01*bb10)/det; 
			k01=(bb02*bb21-bb01*bb22)/det; k10=k01;
			k02=(bb01*bb12-bb02*bb11)/det; k20=k02;
			k12=(bb02*bb10-bb00*bb12)/det; k21=k12;

			// CG on PAP(X)=Pf
			double* X =as;
			double* R =(double*)malloc(sizeof(double)*N);
			double* P =(double*)malloc(sizeof(double)*N);
			double* AP=(double*)malloc(sizeof(double)*N);
			memset(X, 0,sizeof(double)*N);              //X=0 
			memcpy(R,fs,sizeof(double)*N); project(R);  //R=Pf
			memcpy(P, R,sizeof(double)*N);              //P=R
			double rr=inner_product(R,R); bool CG_success=false; int it;
			for(it=0;it<N;it++) {
				if(rr<1E-12) { CG_success=true; break; }
				project(P);
				for(int i=0;i<N;i++) AP[i]=inner_product(A[i],P); 
				project(AP);
				double pap=inner_product(P,AP);
				double a = rr/pap;
				for(int i=0;i<N;i++) {
					X[i] += a* P[i];
					R[i] -= a*AP[i]; }
				double rr_new=inner_product(R,R);
				double b = rr_new/rr; rr=rr_new;
				for(int i=0;i<N;i++) P[i] = R[i] + b*P[i];
			}
			if(CG_success)
			{
				//printf("CG with nullspace converged in %d iterations for %d size\n",it,N);
				// a=P(X) and b=(B^TB)^(-1)B^T(f-Aa)
				for(int i=0;i<N;i++)
					R[i] = fs[i]- inner_product(A[i],X);

				double br0 = inner_product(xs,R); 
				double br1 = inner_product(ys,R); 
				double br2 = inner_product_with_one(R); 
				b0 = k00*br0+k01*br1+k02*br2;
				b1 = k10*br0+k11*br1+k12*br2;
				b2 = k20*br0+k21*br1+k22*br2;
			}
			else
			{
				//printf("Warning! CG failed in RBF with N=%f\n",N);
				b0=b1=b2=0;
				for(int i=0;i<N;i++){ as[i]=0; b2+=fs[i]; } b2/=N;
			}

			// memory cleanup
			free(R); free(P); free(AP);
			for(int i=0;i<N;i++) free(A[i]); free(A);
		}
	}

	double evaluate( double x, double y ) const
	{
		x-=x0; y-=y0; double sum=b0*x+b1*y+b2;
		for(int i=0;i<N;i++)
			sum += as[i]*get_phi(SQRT(SQR(xs[i]-x)+SQR(ys[i]-y)));
		return sum;
	}
};








/*-------------------------------------------------------------------------------------------*\
 * The poly-harmonic interpolation in 3D
 * f(x) = sum_i ai*phi(|x-xi|) + sum_j bj*p_j(x)
\*-------------------------------------------------------------------------------------------*/
class RBF_3D
{
protected:
	ARRAY2D<double> A; 
	VECTOR          X,P,R,AP;
	inline double get_phi(double r) const{ return r;}
	double inner_product( const double* u,const double* v,int N)const{double sum=0;for(int i=N;i>0;i--,u++,v++)sum+=(*u)*(*v);return sum;}
	double inner_product_with_one(        const double* u,int N)const{double sum=0;for(int i=N;i>0;i--,u++    )sum+=(*u)     ;return sum;}
	// K=(B^TB)^(-1)
	double k00,k01,k02,k03; 
	double k10,k11,k12,k13;
	double k20,k21,k22,k23;
	double k30,k31,k32,k33;
	// P = I - BKB^T, B=(X,Y,Z,1)	
	void project(const double* xs, 
		         const double* ys, 
				 const double* zs, double* v, int N) {
		double bv0=0,bv1=0,bv2=0,bv3=0;
		for(int i=N;i>0;i--,xs++,ys++,zs++,v++) { double v_=*v;
			bv0 += (*xs)*v_;
			bv1 += (*ys)*v_;
			bv2 += (*zs)*v_;
			bv3 +=       v_;} xs-=N;ys-=N;zs-=N;v-=N;
		double kbv0 = k00*bv0 + k01*bv1 + k02*bv2 + k03*bv3;
		double kbv1 = k10*bv0 + k11*bv1 + k12*bv2 + k13*bv3;
		double kbv2 = k20*bv0 + k21*bv1 + k22*bv2 + k23*bv3;
		double kbv3 = k30*bv0 + k31*bv1 + k32*bv2 + k33*bv3;
		for(int i=N;i>0;i--,xs++,ys++,zs++,v++) (*v) -= (*xs)*kbv0 + (*ys)*kbv1 + (*zs)*kbv2 + kbv3; }
	
public:
	/*-------------------------------------------------------------------------------------------*\
	 * the expansion point (x0,y0) is the average of the points (xi,yi).
	 * calculate the coefficients ai and bj from the linear system
	 * | A   B | |a| = |f|
	 * | B^T 0 | |b|   |0|,   Aik=phi(|xi-xk|), Bij=p_j(xi-x0)
	 * The linear system is solved by CG on PAPc=Pf, where P=I-B(B^TB)^(-1)B^T
	 * then set a=Pc and b=(B^TB)^(-1)B^T(f-Aa)
	\*-------------------------------------------------------------------------------------------*/
	double interpolate( double* xs,
		                double* ys,
					    double* zs,
					    double* fs, double x, double y, double z, int N)
	{
		//
		if(A.isize<N){A.resize(N,N);X.resize(N);R.resize(N);P.resize(N);AP.resize(N);}
		// find the average (x0,y0,z0) and translate
		double x0=0,y0=0,z0=0;
		for(int i=0;i<N;i++){x0+=xs[i];y0+=ys[i];z0+=zs[i];} x0/=N; y0/=N; z0/=N;
		for(int i=0;i<N;i++){xs[i]-=x0;ys[i]-=y0;zs[i]-=z0;} x-=x0; y-=y0; z-=z0;
		// calculate K=(B^TB)^(-1)
		double bb00 = inner_product(xs,xs,N);
		double bb11 = inner_product(ys,ys,N); 
		double bb22 = inner_product(zs,zs,N);
		double bb33 = N;
		double bb01 = inner_product      (xs,ys,N); double bb10=bb01;
		double bb02 = inner_product      (xs,zs,N); double bb20=bb02;
		double bb12 = inner_product      (ys,zs,N); double bb21=bb12;
		double bb03 = inner_product_with_one(xs,N); double bb30=bb03;
		double bb13 = inner_product_with_one(ys,N); double bb31=bb13;
		double bb23 = inner_product_with_one(zs,N); double bb32=bb23;
		double det = bb00*(bb11*(bb22*bb33 - bb23*bb32) - bb12*(bb21*bb33 - bb23*bb31) + bb13*(bb21*bb32 - bb22*bb31)) 
				   - bb01*(bb10*(bb22*bb33 - bb23*bb32) - bb12*(bb20*bb33 - bb23*bb30) + bb13*(bb20*bb32 - bb22*bb30))
				   + bb02*(bb10*(bb21*bb33 - bb23*bb31) - bb11*(bb20*bb33 - bb23*bb30) + bb13*(bb20*bb31 - bb21*bb30)) 
				   - bb03*(bb10*(bb21*bb32 - bb22*bb31) - bb11*(bb20*bb32 - bb22*bb30) + bb12*(bb20*bb31 - bb21*bb30));
		// in the case when the points are collinear, constant function
		if(ABS(det)<1E-12) { double sum=0; for(int i=0;i<N;i++) sum+=fs[i]; sum/=N; return sum;}
		// calculate the matrix A
		for(int i=0;i<N;i++)
		for(int k=i;k<N;k++){
			double rik = SQRT(SQR(xs[i]-xs[k])
				             +SQR(ys[i]-ys[k])
							 +SQR(zs[i]-zs[k]));
			A[i][k]=A[k][i]=get_phi(rik);
		}
		// calculate K
		k00=( bb11*(bb22*bb33-bb23*bb32)-bb12*(bb21*bb33-bb23*bb31)+bb13*(bb21*bb32-bb22*bb31))/det;
		k11=( bb00*(bb22*bb33-bb23*bb32)-bb02*(bb20*bb33-bb23*bb30)+bb03*(bb20*bb32-bb22*bb30))/det;
		k22=( bb00*(bb11*bb33-bb13*bb31)-bb01*(bb10*bb33-bb13*bb30)+bb03*(bb10*bb31-bb11*bb30))/det; 
		k33=( bb00*(bb11*bb22-bb12*bb21)-bb01*(bb10*bb22-bb12*bb20)+bb02*(bb10*bb21-bb11*bb20))/det;
		k01=(-bb01*(bb22*bb33-bb23*bb32)+bb02*(bb21*bb33-bb23*bb31)-bb03*(bb21*bb32-bb22*bb31))/det; k10=k01;
		k02=( bb01*(bb12*bb33-bb13*bb32)-bb02*(bb11*bb33-bb13*bb31)+bb03*(bb11*bb32-bb12*bb31))/det; k20=k02;
		k03=(-bb01*(bb12*bb23-bb13*bb22)+bb02*(bb11*bb23-bb13*bb21)-bb03*(bb11*bb22-bb12*bb21))/det; k30=k03;
		k12=(-bb00*(bb12*bb33-bb13*bb32)+bb02*(bb10*bb33-bb13*bb30)-bb03*(bb10*bb32-bb12*bb30))/det; k21=k12;
		k13=( bb00*(bb12*bb23-bb13*bb22)-bb02*(bb10*bb23-bb13*bb20)+bb03*(bb10*bb22-bb12*bb20))/det; k31=k13;
		k23=(-bb00*(bb11*bb23-bb13*bb21)+bb01*(bb10*bb23-bb13*bb20)-bb03*(bb10*bb21-bb11*bb20))/det; k32=k23;
		// CG on PAP(X)=Pf
		memset(X, 0,sizeof(double)*N);              //X=0 
		memcpy(R,fs,sizeof(double)*N); project(xs,ys,zs,R,N);  //R=Pf
		memcpy(P, R,sizeof(double)*N);              //P=R
		double rr=inner_product(R,R,N); int it; bool CG_success=false;
		for(it=0;it<N;it++) {
			if(rr<1E-9) { CG_success=true; break; }
			for(int i=0;i<N;i++) AP[i]=inner_product(A[i],P,N); project(xs,ys,zs,AP,N);
			double pap=inner_product(P,AP,N);
			double a = rr/pap;
			for(int i=0;i<N;i++) {
				X[i] += a* P[i];
				R[i] -= a*AP[i]; }
			double rr_new=inner_product(R,R,N);
			double b = rr_new/rr; rr=rr_new;
			for(int i=0;i<N;i++) P[i] = R[i] + b*P[i];
		}
		if(CG_success)
		{
			//printf("CG with nullspace converged in %d iterations for %d size\n",it,N);
			// a=P(X) and b=(B^TB)^(-1)B^T(f-Aa)
			for(int i=0;i<N;i++) R[i] = fs[i]- inner_product(A[i],X,N);
			double br0 = inner_product(xs,R,N); 
			double br1 = inner_product(ys,R,N); 
			double br2 = inner_product(zs,R,N);
			double br3 = inner_product_with_one(R,N); 
			double b0 = k00*br0+k01*br1+k02*br2+k03*br3;
			double b1 = k10*br0+k11*br1+k12*br2+k13*br3;
			double b2 = k20*br0+k21*br1+k22*br2+k23*br3;
			double b3 = k30*br0+k31*br1+k32*br2+k33*br3;
			double sum=b0*x+b1*y+b2*z+b3;
			for(int i=0;i<N;i++)
				sum += X[i]*get_phi( SQRT(SQR(xs[i]-x)
										 +SQR(ys[i]-y)
										 +SQR(zs[i]-z)));
			return sum;
		}
		else{double sum=0; for(int i=0;i<N;i++) sum+=fs[i]; sum/=N; return sum;}
	}
};

/*
//-------------------------------------------------------------------------------------------
// The poly-harmonic interpolation in 3D
// f(x) = sum_i ai*phi(|x-xi|) + sum_j bj*p_j(x)
//-------------------------------------------------------------------------------------------
class RBF_3D
{
protected:
	ARRAY2D<double> A ; 
	VECTOR          as;
	inline double get_phi(double r) const{ return r;}
public:
	//-------------------------------------------------------------------------------------------
	// the expansion point (x0,y0) is the average of the points (xi,yi).
	// calculate the coefficients ai and bj from the linear system
	// | A   B | |a| = |f|
	// | B^T 0 | |b|   |0|,   Aik=phi(|xi-xk|), Bij=p_j(xi-x0)
	// The linear system is solved by CG on PAPc=Pf, where P=I-B(B^TB)^(-1)B^T
	// then set a=Pc and b=(B^TB)^(-1)B^T(f-Aa)
	//-------------------------------------------------------------------------------------------
	double interpolate( const double* xs,
		                const double* ys,
					    const double* zs,
					    const double* fs, double x, double y, double z, int N)
	{
		//
		if(A.isize<N+4)
		{ 
			A.resize(N+4,N+4); 
			as.resize(N+4); 
		}
		// translation
		double x0=0,y0=0,z0=0;
		for(int i=0;i<N;i++){x0+=xs[i];
		                     y0+=ys[i];
							 z0+=zs[i];} x0/=N; 
							             y0/=N; 
										 z0/=N;
		x-=x0;
		y-=y0;
		z-=z0;
		// calculate the matrix A
		for(int i=0;i<N+4;i++)
		for(int k=i;k<N+4;k++){
			if(i<N && k<N) {
				double rik = SQRT(SQR(xs[i]-xs[k])
								 +SQR(ys[i]-ys[k])
								 +SQR(zs[i]-zs[k]));
				A[i][k]=get_phi(rik); }
			else if(i<N && k>=N) {
				     if(k==N  ) A[i][k] = xs[i]-x0;
				else if(k==N+1) A[i][k] = ys[i]-y0;
				else if(k==N+2) A[i][k] = zs[i]-z0;
				else            A[i][k] = 1; }
			else A[i][k]=0; 
			A[k][i]=A[i][k];
		}
		// gaussian elimination
		memcpy((double*)as,fs,sizeof(double)*N); as[N]=as[N+1]=as[N+2]=as[N+3]=0;
		if(MATRIX::gaussian_elimination(A.data_2d,as,N+4)==false) {
			// in the case when the points are collinear, constant function
			for(int i=0;i<N;i++){ as[i]=0; as[N+3]+=fs[i]; } as[N+3]/=N;
			as[N]=as[N+1]=as[N+2]=0;
		}
		// evaluation
		double sum=as[N]*x+as[N+1]*y+as[N+2]*z+as[N+3];
		for(int i=0;i<N;i++)
			sum += as[i]*get_phi(SQRT(SQR(xs[i]-x)
			                         +SQR(ys[i]-y)
									 +SQR(zs[i]-z)));
		return sum ; 
	}
};
*/
#endif
