#include <omp.h>
#include "../Macros.h"
#include "SparseSymmMatrix.h"

//----------------------------------------------------------------------------
// find index
//----------------------------------------------------------------------------
int SparseSymmMatrix::find_index( int i, int j ) const {
	if(i==j) return i;
	if(j< i){ int temp=i;i=j;j=temp; }
	int ai  =JA[i  ];
	int aip1=JA[i+1];
	for(;ai<aip1 && JA[ai]<j;ai++) ;
	return ai;}
//----------------------------------------------------------------------------
// constructor
//----------------------------------------------------------------------------
SparseSymmMatrix::SparseSymmMatrix(int N_) {
	N=N_; M=0; Mmax=500; imax=0;
	JA.resize(N+Mmax+1); 
	AA.resize(N+Mmax+1); for(int i=0;i<=N+1;i++){ JA[i]=N+1; AA[i]=0; }}

void SparseSymmMatrix::set_dimension (int N_){
	N=N_; M=0; Mmax=500; imax=0;
	JA.resize(N+Mmax+1); 
	AA.resize(N+Mmax+1); for(int i=0;i<=N+1;i++){ JA[i]=N+1; AA[i]=0; }}
//----------------------------------------------------------------------------
// get_val
//----------------------------------------------------------------------------
double SparseSymmMatrix::get_val(int i,int j) const {
	if(i==j) return AA[i];
	else { int a=find_index(i,j); if(a<=N+M && JA[a]==j) return AA[a]; else return 0; }}
//----------------------------------------------------------------------------
// add_element
//----------------------------------------------------------------------------
void SparseSymmMatrix::add_element(int i,int j, double v)
{
	if(v==0) return;
	if(i==j){ AA[i]+=v; return;}
	if(i> j) return;
	// prepare new row
	if(i>imax){for(int k=imax+2;k<=i+1;k++) JA[k]=JA[imax+1]; imax=i; }
	int a = find_index(i,j);
	if(a<=N+M && JA[a]==j){ AA[a]+=v; return; }
	// new index, increase the array when full
	if(M==Mmax) {
		Mmax=Mmax+1+Mmax/2;
		JA.resize_with_existing_data_kept(N+Mmax+1); 
		AA.resize_with_existing_data_kept(N+Mmax+1); }
	for(int n=N+M;n>=a;n--){AA[n+1]=AA[n];JA[n+1]=JA[n];}
	AA[a]=v;JA[a]=j;M++;
	// IA[i+1:N]++
	for(int n=i+1;n<=imax+1;n++) JA[n]++;
}
//----------------------------------------------------------------------------
// apply_A
//----------------------------------------------------------------------------
void SparseSymmMatrix::apply_A( const double* X, double* Y ) const
{
	// diagonal
	int i;
	const double* AA_=AA.data;
	#pragma omp parallel for private(i) firstprivate(X,Y,AA_)
	for(i=0;i<N;i++) 
		Y[i]=AA_[i]*X[i];
	// upper and lower triangular 
	int NT;
	#pragma omp parallel 
	{
		// load distribution
		#pragma omp single
		{NT=omp_get_num_threads();}
		int T = omp_get_thread_num();
		int ibeg_T =     (int(N/NT)+1)* T        ;
		int iend_T = MIN((int(N/NT)+1)*(T+1),N)-1;
		// i : ibeg~iend
		const int*    IA_ = JA+ibeg_T ;
		const int*    JA_ = JA+IA_[0];
		const double* AA_ = AA+IA_[0];
		for(int i=ibeg_T;i<=iend_T;i++)
		{
			double sum=0;
			double xi =X[i];
			int ai   = *IA_; IA_++;
			int aip1 = *IA_; int count=aip1-ai; int j; double aij;
			while(count>0)
			{
				j  =*JA_;JA_++;
				aij=*AA_;AA_++;
				sum  += aij*X[j]; // upper triangular
				Y[j] += aij*xi  ; // lower triangular
				count--; 
			}
			Y[i] += sum;
		}
	}
}
//----------------------------------------------------------------------------
// inner product
//----------------------------------------------------------------------------
double SparseSymmMatrix::inner_product( const double* X, const double* Y, int N ) 
{
	int NT; double sum=0;
	#pragma omp parallel shared(NT,X,Y,N) reduction(+:sum)
	{
		// load distribution
		#pragma omp single
		{NT=omp_get_num_threads();}
		int T = omp_get_thread_num();
		int ibeg_T =     (int(N/NT)+1)* T        ;
		int iend_T = MIN((int(N/NT)+1)*(T+1),N)-1;
		// i : ibeg~iend
		const double* X_ = X+ibeg_T;
		const double* Y_ = Y+ibeg_T;
		for(int i=iend_T-ibeg_T;i>=0;i--,X_++,Y_++) sum+=(*X_)*(*Y_);
	}
	return sum;
}
//----------------------------------------------------------------------------
// calculate the MILU preconditioner
//----------------------------------------------------------------------------
void SparseSymmMatrix::calculate_MILU() 
{
	// D=D0
	D.resize(N);
	memcpy(D,AA,sizeof(double)*N);
	//
	for(int k=0;k<N;k++)
	{
		int ak  =JA[k  ];
		int akp1=JA[k+1];
		double uk_=0; double dkk=D[k];
		for(int i_=ak;i_<akp1;i_++) uk_+=AA[i_];
		for(int i_=ak;i_<akp1;i_++)
		{
			int      i = JA[i_];
			double uki = AA[i_];
			D[i] -= uki*(uki+.97*(uk_-uki))/dkk;
		}
		if(ABS(dkk)<1E-5) dkk=1E-5;
	}
}
//----------------------------------------------------------------------------
// M=(L+D)/D*(U+D), and 1/M=(U+D)^(-1)*(L/D+I)^(-1)
//----------------------------------------------------------------------------
void SparseSymmMatrix::apply_M( double* X ) const
{
	// (L/D+I)^{-1} * X
	for(int j=0;j<N;j++)
	{
		double dj = D[j];
		double xj = X[j]; xj/=dj;
		int aj = JA[j  ]; const int   * JA_=JA+aj;
		int ajp1=JA[j+1]; const double* AA_=AA+aj;
		for(int i_=ajp1-aj;i_>0;i_--,JA_++,AA_++)
			X[*JA_] -= (*AA_)*xj;
	}
    // (U+D)^{-1} * x
    for(int i=N-1;i>=0;i--){
        double sum = X[i];
        int aj   = JA[i  ]; const int   * JA_=JA+aj;
        int ajp1 = JA[i+1]; const double* AA_=AA+aj;
        for(int j_=ajp1-aj; j_>0;j_--,JA_++,AA_++)
            sum -= (*AA_)*X[*JA_];
        X[i] = sum/D[i];}
}
//----------------------------------------------------------------------------
// projection into the sum-zero space
//----------------------------------------------------------------------------
void SparseSymmMatrix::projection_sum_zero( double* X ) const
{
	int NT; double sum=0;
	#pragma omp parallel shared(NT,X) reduction(+:sum)
	{
		// load distribution
		#pragma omp single
		{NT=omp_get_num_threads();}
		int T = omp_get_thread_num();
		int ibeg_T =     (int(N/NT)+1)* T        ;
		int iend_T = MIN((int(N/NT)+1)*(T+1),N)-1;
		// i : ibeg~iend
		const double* X_ = X+ibeg_T;
		for(int i=iend_T-ibeg_T;i>=0;i--,X_++) sum+=(*X_);
	}
	double avg=sum/N;
	#pragma omp parallel shared(NT,X) reduction(+:sum)
	{
		// load distribution
		#pragma omp single
		{NT=omp_get_num_threads();}
		int T = omp_get_thread_num();
		int ibeg_T =     (int(N/NT)+1)* T        ;
		int iend_T = MIN((int(N/NT)+1)*(T+1),N)-1;
		// i : ibeg~iend
		double* X_ = X+ibeg_T;
		for(int i=iend_T-ibeg_T;i>=0;i--,X_++) (*X_)-=avg;
	}
}

//----------------------------------------------------------------------------
// The PCG routine
//----------------------------------------------------------------------------
void SparseSymmMatrix::PCG( const VECTOR& B, VECTOR& X, bool with_initial_guess,
						               bool nonzero_kernel,
                                       int  max_it        , 
						               double tolerance )
{
    calculate_MILU(); // preconditioner
			const double* pB  = (const double*)B;
	              double* pX  = (      double*)X;
	VECTOR  R(N); double* pR  = (      double*)R;
	VECTOR  Z(N); double* pZ  = (      double*)Z;
	VECTOR  P(N); double* pP  = (      double*)P;
	VECTOR AP(N); double* pAP = (      double*)AP;
	//---------------------------------------------------------------------
	// initilialize R,P,R0,  solve MPAx=MPb
	//---------------------------------------------------------------------
	if(with_initial_guess==false) memset((void*)pX,0,sizeof(double)*N);
	apply_A(X,R); for(int i=0;i<N;i++) pR[i] = pB[i] - pR[i];  // R=B-AX
	if(nonzero_kernel) projection_sum_zero(R);
	memcpy(pZ,pR,sizeof(double)*N); apply_M(Z); // Z=invert(M)*R
	memcpy(pP,pZ,sizeof(double)*N);               // P=Z;  
	double  rz = inner_product(R,Z,N); 
	double rz0 = ABS(rz);
	//---------------------------------------------------------------------
	// iteration
	//---------------------------------------------------------------------
	int iteration = 0; 
	while( iteration++ < max_it ) {
		// a
		apply_A(P,AP); if(nonzero_kernel) projection_sum_zero(AP);
		double app = inner_product(AP,P,N);
		double a   = rz/app;
		// check convergence
		if(ABS(rz)<tolerance) {
			printf( "MICCG converged in %d iterations\n", iteration ); return; }
		// x & r
		for(int i=N;i>0;i--,pX++,pP++,pR++,pAP++){
			*pX+= a*(*pP );
			*pR-= a*(*pAP);
		}
		pX-=N; pP -=N; pR-=N; pAP-=N;
		// b
		memcpy(pZ,pR,sizeof(double)*N);
		apply_M(Z);
		double rz_new = inner_product(R,Z,N);
		double b      = rz_new/rz; rz=rz_new;
		for(int i=N;i>0;i--,pZ++,pP++)
			*pP = (*pZ)+b*(*pP);
		pP-=N; pZ-=N;
	}
	printf( "MICCG failed to converge in maximum iteration %d\n", iteration );
}

//---------------------------------------------------------------------
// print
//---------------------------------------------------------------------
void SparseSymmMatrix::print(const char* name) const {
	printf("-------------------SPARSE Symm MATRIX : %s--------------\n",name);
	for(int i=0;i<N;i++) {
		printf("i:%3d | ",i);
		int index_i  =JA[i  ];
		int index_ip1=JA[i+1];
		printf("(%3d,% 3.5f), ", i,AA[i]);
		for(int a=index_i;a<index_ip1;a++){
			int    j = JA[a];
			double v = AA[a];
			printf("(%3d,% 3.5f), ", j,v);}
		printf("\n");}}

void SparseSymmMatrix::print_matlab(const char* filename) const {
	FILE* fp=fopen(filename,"wt");
	fprintf(fp,"A=[\n");
	for(int i=0;i<N;i++) {
		for(int j=0;j<N;j++) {
			if(j<N-1) fprintf(fp,"%f,",get_val(i,j));
			else      fprintf(fp,"%f;",get_val(i,j));
		}
		if(i<N-1) fprintf(fp,"\n");
		else      fprintf(fp,"] ");
	}
	fclose(fp);
}