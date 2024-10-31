#ifndef MATRIX_H
#define MATRIX_H

#include <malloc.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include "../Macros.h"
#include "../ARRAY.h"
#include "../ARRAY2D.h"
#include "../geometry/POINT3.h"
#include "VECTOR.h"

class MATRIX : public ARRAY2D<double>
{
public:
	MATRIX(int row=1, int col=1) : ARRAY2D<double>(row,col){}

	//-------------------------------------------------------------------------
	// accessors
	//-------------------------------------------------------------------------
	int row() const{ return isize; }
	int col() const{ return jsize; }

	inline double  operator()(int i, int j) const{ return data[i*jsize+j];	}
	inline double& operator()(int i, int j)      { return data[i*jsize+j];	}

	void print(const char* name) const {
		printf("--------------------%s--------------------\n", name);
		for(int i=0;i<isize;i++) {
			for(int j=0;j<jsize;j++)
				printf("% .3f", data[i*jsize+j]);
			printf("\n"); }
		printf("\n"); }

	//-------------------------------------------------------------------------
	// unary operations : +=,-=,*=,/=
	//-------------------------------------------------------------------------
	void operator+=(const MATRIX& B){ ARRAY2D<double>::operator+=(B);}
	void operator-=(const MATRIX& B){ ARRAY2D<double>::operator-=(B);}
	void operator*=(const double& scalar){ ARRAY2D<double>::operator*=(scalar);}
	void operator/=(const double& scalar){ ARRAY2D<double>::operator/=(scalar);}

	//-------------------------------------------------------------------------
	// binary operations : +,-,*,/
	//-------------------------------------------------------------------------
	MATRIX operator+(const MATRIX& B) const{MATRIX C(*this); C+=B; return C;}
	MATRIX operator-(const MATRIX& B) const{MATRIX C(*this); C-=B; return C;}
	MATRIX operator*(double        s) const{MATRIX C(*this); C*=s; return C;}
	MATRIX operator/(double        s) const{MATRIX C(*this); C/=s; return C;}
	
	MATRIX operator*(const MATRIX& B) const	// matrix = matrix * matrix
	{
		assert( jsize == B.isize ); int ksize=B.jsize;
		MATRIX C(isize, ksize);
		for(int i=0;i<isize;i++){ 
			      double* Ci = C.data + i*ksize;
			const double* Ai =   data + i*jsize;
			
			for(int k=0;k<ksize;k++){
				double sum=0;
				const double* Bk = B.data + k      ;
				for(int j=jsize;j>0;j--,Ai++,Bk+=ksize)
					sum += (*Ai)*(*Bk);
				*Ci=sum; Ci++; Ai-=jsize;
			}
		}
		return C;
	}

	VECTOR operator*(const VECTOR& X) const	{
		assert( jsize == X.size );
		VECTOR B(isize);
		const double* X_ = (const double*)X;
		      double* B_ = (      double*)B;
		for(int i=0;i<isize;i++){
			double sum=0.;
			const double* Ai = data+i*jsize;
			for(int j=jsize;j>0;j--,Ai++,X_++)
				sum += (*Ai)*(*X_);
			*B_=sum; B_++; X_-=jsize;
		}
		return B;
	}

	//-------------------------------------------------------------------------
	// LU decomposition, PA=LU , Algorithm 3.4.2 (Golub) modified to partial pivot
	// L has only one in the diagonal.
	// To save memory, L and U are combined in the one matrix output
	// P represents the permutation : for each k, swap k and P(k)
	//-------------------------------------------------------------------------
	static void LU_decomposition( MATRIX& A, ARRAY<int>& P ) {
		assert(A.isize==A.jsize); int n=A.isize; double* elmt = A.data;
		P.resize(n); 
		double* buff = new double[n];

		for(int k=0;k<n;k++) {
			// pivoting
			{
				//select i>=k to maximize |uik|
				double  max=0; int i; double* U_k=elmt+k*n+k;
				for(int j=k;j<n;j++){
					double ujk = *U_k; U_k+=n; 
					if(ABS(ujk)>max){ max=ujk; i=j;}}
				 P[k]=i;
				// swap
				if(i>k){
					memcpy(buff    ,elmt+k*n,sizeof(double)*n);
					memcpy(elmt+k*n,elmt+i*n,sizeof(double)*n);
					memcpy(elmt+i*n,buff    ,sizeof(double)*n);}
			}

			double ukk = elmt[k*n+k];

			if(ABS(ukk)<1E-20) printf("Warning : matrix in LU is nearly singular!\n");
			else {
				for(int j=k+1;j<n;j++){
					double ljk = elmt[j*n+k]/ukk; elmt[j*n+k] = ljk;
					double* Uj = elmt+j*n+k+1;
					double* Uk = elmt+k*n+k+1; int count=n-(k+1);
					while(count>0){ (*Uj)-=ljk*(*Uk); count--; Uj++; Uk++;}
				}
			}
		}
		delete[] buff;
	}

	//------------------------------------------------------------------------------
	// Gaussian elimination
	//------------------------------------------------------------------------------
	static bool gaussian_elimination( double** A, double* B, int n ) 
	{
		// prepare buffer
		VECTOR buff(n);
		// upper triangulation
		for(int i=0;i<n;i++) {
			// pivoting : A[j][i] = max A[i:n-1][i]
			{   int j=i; double max=ABS(A[i][i]);
				for(int k=i+1;k<n;k++) if(ABS(A[k][i])>max){j=k;max=ABS(A[k][i]);}
				// swap A[i][i:n-1] and A[j][i:n-1]
				memcpy(buff  ,A[i]+i,sizeof(double)*(n-i));
				memcpy(A[i]+i,A[j]+i,sizeof(double)*(n-i));
				memcpy(A[j]+i,buff  ,sizeof(double)*(n-i));
			}
			// upper triangulartion
			double aii = A[i][i];
			if(ABS(aii)<1E-10) return false; // singular matrix
			for(int j=i+1;j<n;j++){
				double* ai = A[i]+i; 
				double* aj = A[j]+i; double aji = *aj;
				for(int k=n-i-1;k>0;k--,ai++,aj++) (*aj)-=(*ai)*aji/aii;
				                                    B[j]-= B[i]*aji/aii;
			}
		}
		// back substitution
		for(int i=n-1;i>=0;i--){
			double sum=B[i]; double* ai=A[i];
			for(int j=i+1;j<n;j++) sum-=ai[j]*B[j];
			B[i]=sum/A[i][i];
		}
		return true;
	}

	static void solve_by_LU( const MATRIX& LU, const ARRAY<int>& P, VECTOR& B ) 
	{
		// B := PB
		int n = LU.isize; 
		   double* B_ = (   double*)B;
		const int* P_ = (const int*)P;
		const double* elmt  = LU.data;
		for(int k=0;k<n;k++) {
			int i=P[k]; 
			if(i!=k){ double t=B_[i];B_[i]=B_[k];B_[k]=t;}}
		
		// B := L^(-1)B
		for(int i=0;i<n;i++){
			double sum=B_[i];
			const double* Li = elmt + i*n;
			for(int j=0;j<i;j++) sum-=Li[j]*B_[j];
			B_[i]=sum;}

		// B := U^(-1)B
		for(int i=n-1;i>=0;i--){
			double sum=B_[i];
			const double* Ui = elmt + i*n;
			for(int j=i+1;j<n;j++) sum-=Ui[j]*B_[j];
			B_[i]=sum/Ui[i];}
	}

	void inverse_by_LU(MATRIX& B) const
	{
		assert(isize==jsize); int n=isize;
		VECTOR X(n); B.resize(n,n); double* elmt = B.data;
		MATRIX LU(*this); ARRAY<int> P; LU_decomposition(LU,P);
		for(int i=0;i<n;i++){
			X=0; X[i]=1; solve_by_LU(LU,P,X);
			      double* B_i = elmt + i;
			const double* X_  = (const double*)X;
			for(int j=n;j>0;j--,B_i+=n,X_++)
				*(B_i) = *(X_);
		}
	}
	
	//-------------------------------------------------------------------------
	// calculate the determinant from the LU decomposition
	//-------------------------------------------------------------------------
	double determinant() const
	{
		// determinant of LU
		MATRIX LU(*this); ARRAY<int> P; LU_decomposition(LU,P);
		double product = 1;
		for(int i=0;i<isize;i++)
			product *= LU(i,i);
		// determinant of P
		for(int i=  0;i<isize;i++)
			if(P[i]!=i) product = -product;
		return product;
	}

	MATRIX transpose() const {
		MATRIX B(jsize, isize);
		for(int j=0;j<jsize;j++){
			      double* Bj = B.data + j*isize;
			const double* A_j=   data + j      ;

			for(int i=isize;i>0;i--,Bj++,A_j+=jsize)
				*Bj = *A_j;
		}
		return B;
	}

	//-------------------------------------------------------------------------
	//  QR factorization
	//-------------------------------------------------------------------------
	/*
	void qr_factorization( MATRIX& Q, MATRIX& R ) const
	{
		assert( isize==jsize ); int N=isize;
		R=*this; // R=A
		Q.resize(N,N); 
		for(int i=0;i<N;i++) 
		for(int j=0;j<N;j++)
			Q[i][j]= (i==j)?1:0; // Q=I

		VECTOR a(N);
		VECTOR v(N);
		VECTOR Rv(N);
		VECTOR Qv(N);
		for(int j=0;j<N-1;j++)
		{
			// a=R(j:N-1,j)
			int M = N-j;
			for(int i=0;i<M;i++) a[i]=R[i+j][j];

			// v = a + sgn(a1)*length(a)e1
			double al=0; for(int i=0;i<M;i++) al+=a[i]*a[i]; al=SQRT(al); v=a; v[0] += SGN(a[0])*al;
			double vl=0; for(int i=0;i<M;i++) vl+=v[i]*v[i]; vl=SQRT(vl); v/=vl;

			// Rv 
			for(int i=0;i<M;i++)
			{
				double sum=0; for(int k=0;k<M;k++) sum+=R[k+j][i+j]*v[k]; Rv[i]=sum;
			}
					
			// R=R-2v(Rv)^T
			for(int i=0;i<M;i++)
			for(int k=0;k<M;k++)
				R[i+j][k+j] -= 2* v[i]*Rv[k];
			
			// Qv
			for(int i=0;i<N;i++)
			{
				double sum=0; for(int k=0;k<M;k++) sum+=Q[i][k+j]*v[k]; Qv[i]=sum;
			}
			
			// Q=Q-2(Qv)v^T
			for(int i=0;i<N;i++)
			for(int k=0;k<M;k++)
				Q[i][k+j] -= 2*Qv[i]*v[k];
		}
	}
	*/

	
};

static MATRIX operator*(const double& scalar, const MATRIX& A){
	MATRIX C(A); C*=scalar; return C;
}

/*-----------------------------------------------------------------------------------------------*\
 * inverse of a 2x2 real matrix
\*-----------------------------------------------------------------------------------------------*/
static MATRIX inverse_2_by_2( const MATRIX& A )
{
	MATRIX B(2,2);
	
	double a11=A(0,0), a12=A(0,1);
	double a21=A(1,0), a22=A(1,1);
	
	double det=a11*a22-a12*a21;
	
	B(0,0)=  (a22)/det; B(0,1)= -(a12)/det;
	B(1,0)= -(a21)/det; B(1,1)=  (a11)/det;
	
	return B;
}

/*-----------------------------------------------------------------------------------------------*\
 * inverse of a 3x3 real matrix
\*-----------------------------------------------------------------------------------------------*/
static MATRIX inverse_3_by_3( const MATRIX& A )
{
	MATRIX B(3,3);
	
	double a11=A(0,0), a12=A(0,1), a13=A(0,2);
	double a21=A(1,0), a22=A(1,1), a23=A(1,2);
	double a31=A(2,0), a32=A(2,1), a33=A(2,2);
	
	double det=(a11*a22-a12*a21)*a33+(a13*a21-a11*a23)*a32+(a12*a23-a13*a22)*a31;
	
	B(0,0)=  (a22*a33-a23*a32)/det;
	B(0,1)= -(a12*a33-a13*a32)/det;
	B(0,2)=  (a12*a23-a13*a22)/det;
	B(1,0)= -(a21*a33-a23*a31)/det;
	B(1,1)=  (a11*a33-a13*a31)/det;
	B(1,2)= -(a11*a23-a13*a21)/det;
	B(2,0)=  (a21*a32-a22*a31)/det;
	B(2,1)= -(a11*a32-a12*a31)/det;
	B(2,2)=  (a11*a22-a12*a21)/det;

	return B;
}

static POINT3 operator*(const MATRIX& A, const POINT3& X )
{
	POINT3 AX;
	AX.x = A(0,0)*X.x + A(0,1)*X.y + A(0,2)*X.z;
	AX.y = A(1,0)*X.x + A(1,1)*X.y + A(1,2)*X.z;
	AX.z = A(2,0)*X.x + A(2,1)*X.y + A(2,2)*X.z;
	return AX;
}

static POINT2 operator*(const MATRIX& A, const POINT2& X )
{
	POINT2 AX;
	AX.x = A(0,0)*X.x + A(0,1)*X.y;
	AX.y = A(1,0)*X.x + A(1,1)*X.y;
	return AX;
}

static void solve_by_LU( const MATRIX& LU, const ARRAY<int>& P, POINT3& B ) 
{
	// B := PB
	double B_[3]; B_[0]=B.x; 
	              B_[1]=B.y; 
				  B_[2]=B.z;
	const int* P_ = (const int*)P;
	const double* elmt  = LU.data;
	for(int k=0;k<3;k++) {
		int i=P[k]; 
		if(i!=k){ double t=B_[i];B_[i]=B_[k];B_[k]=t;}}
	
	// B := L^(-1)B
	for(int i=0;i<3;i++){
		double sum=B_[i];
		const double* Li = elmt + i*3;
		for(int j=0;j<i;j++) sum-=Li[j]*B_[j];
		B_[i]=sum;}

	// B := U^(-1)B
	for(int i=2;i>=0;i--){
		double sum=B_[i];
		const double* Ui = elmt + i*3;
		for(int j=i+1;j<3;j++) sum-=Ui[j]*B_[j];
		B_[i]=sum/Ui[i];}
	B.x=B_[0];
	B.y=B_[1];
	B.z=B_[2];
}
#endif
