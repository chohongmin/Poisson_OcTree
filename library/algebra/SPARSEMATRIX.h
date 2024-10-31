#ifndef _SPARSE_MATRIX_CRS_
#define _SPARSE_MATRIX_CRS_

#include "VECTOR.h"

//---------------------------------------------------------------------
// SPARSEMATRIX in the "Compressed Sparse Row" Storage format
//
// Usage :  clear_all();
//          add_element(); to construct the matrix
//          
// Important : For faster creation of large size matrix, start with N=1 
//             and put elements sequencitially from low row to high row.
//
// 0-based index
//
// by Chohong Min, 2006
//---------------------------------------------------------------------
class SPARSEMATRIX 
{
protected:
	int*    IA; int N; int Nmax;// IA:0~N  -> 0~M  , NxN matrix
	int*    JA; int M; int Mmax;// JA:0~M-1-> 0~N-1
	double* AA;                 // AA:0~M-1-> R

public:
	void clear_all(){
		if(IA!=0)delete[] IA; IA=0; N=0; Nmax=0;
		if(JA!=0)delete[] JA; JA=0; M=0; Mmax=0;
		if(AA!=0)delete[] AA; AA=0; }
	
	 SPARSEMATRIX( ) { 
		Nmax=100; N=0; IA=new int   [Nmax+1]; IA[0]=IA[1]=0;
		Mmax=100; M=0; JA=new int   [Mmax  ];
		               AA=new double[Mmax  ]; }
	~SPARSEMATRIX(){ clear_all();}
	
protected:
	void ensure_IA_has_enough_size(int imax) {
		if(N>imax) // already enough
			;
		else { 	// increase the buffer
			if(Nmax==(imax+1)) {
				int Nmaxnew = (imax+1) + ((imax+1)/2); // 50% buffering
				int* IAnew = new int[Nmaxnew+1];
				memcpy(IAnew,IA,sizeof(int)*(N+1));	delete[] IA; 
				IA=IAnew; Nmax=Nmaxnew; }
			// update IA : IA[N:imax+1]=M;
			for(int n=N;n<=imax+1;n++) IA[n]=M;
			N=imax+1; }}

	void ensure_JA_AA_not_full() {
		if(M==Mmax) {
			int Mmaxnew = Mmax+1 + (Mmax/2); // 50% buffering
			int*    JAnew = new int   [Mmaxnew];
			double* AAnew = new double[Mmaxnew];
			memcpy(JAnew,JA,sizeof(int   )*(M));
			memcpy(AAnew,AA,sizeof(double)*(M));
			delete[] JA; JA=JAnew; Mmax=Mmaxnew;
			delete[] AA; AA=AAnew; }}
	
	// i,j    : 0~N-1
	// return : 0~M-1
	inline int find_index( int i, int j ) const {
		int ai   = IA[i  ];
		int aip1 = IA[i+1];
		for(int a=ai;a<aip1;a++) 
			if(JA[a]==j)
				return a;		
		return -1; // this means that the index does not exist
	}
	
public:
	int dim() const{ return N; }
	int number_of_rows    () const {return N;}
    int number_of_elements() const {return M;}
    const    int* get_IA() const { return IA; }
    const    int* get_JA() const { return JA; }
    const double* get_AA() const { return AA; }
    double get_val( int i, int j ) const {
		int a = find_index(i,j);
		if(a==-1) return 0;  // no entry means zero value
		else      return AA[a]; }

	void add_element(int i, int j, double v ){
		if(ABS(v)<1E-10) return;
		// array size of IA
		if(i>=N) ensure_IA_has_enough_size(i);
		// search the insertion location
		int ai   = IA[i  ];
		int aip1 = IA[i+1];
		int a    = ai; // insertion point
		for(;a<aip1 && JA[a]<j; a++) ;
		// existing index 
		if(a<aip1 && JA[a]==j) AA[a]+=v;
		// new index
		else {
			ensure_JA_AA_not_full();
			// JA[a+1:M]=JA[a:M-1]; JA[a]=j;
			// AA[a+1:M]=AA[a:M-1]; AA[a]=v;
			for(int n=M-1;n>=a;n--) {
				JA[n+1] = JA[n];
				AA[n+1] = AA[n];}
			JA[a]=j;M++;
			AA[a]=v;
			
			// IA[i+1:N]++
			for(int n=i+1;n<=N;n++) IA[n]++;}}

	//---------------------------------------------------------------------
	// matrix multiplication
	//---------------------------------------------------------------------
	void apply_A(const VECTOR& X, VECTOR& AX) const {
		assert( X.size==N);
		if(    AX.size!=N) AX.resize(N);
		const double* pX  = (const double*) X;  
		      double* pAX = (      double*)AX; 
		// optimized 
		const int*    IA_ = IA;
		const int*    JA_ = JA;
		const double* AA_ = AA;
		for(int i=N-1;i>=0;i--) {
			double sum = 0;
			int ai   = *IA_; IA_++;
			int aip1 = *IA_; 
			int count = aip1-ai; int j; double aij;
			while(count>0) {
				  j = *JA_; JA_++;
				aij = *AA_; AA_++;
				sum += aij*pX[j];
				count--; }
			*pAX = sum; pAX++;}}

	//---------------------------------------------------------------------
	// matrix multiplication
	//---------------------------------------------------------------------
	bool is_symmetric() const {
		// optimized 
		const int*    IA_ = IA;
		const int*    JA_ = JA;
		const double* AA_ = AA;
		for(int i=0;i<=N-1;i++) {
			int ai   = *IA_; IA_++;
			int aip1 = *IA_; 
			int count = aip1-ai; int j; double aij;
			while(count>0) {
				  j = *JA_; JA_++;
				aij = *AA_; AA_++;
				if(j>i)
				{
					double aji = get_val(j,i);
					if(ABS(aij-aji)>1E-15)
						return false;
				}
				count--; }}
		return true;
	}

	enum Preconditioner
	{
		ILU,
		MILU
	};

protected:
	// A=E+D0+F, M=(E+D)(1/D)(F+D), match the diagonals 
	// Dii = D0ii - sum_{j=0}^{i-1}Eij*Fji/Djj
	void calculate_ILU( VECTOR& D )
	{
		for(int i=0;i<N;i++) {
			int ai   = IA[i  ];
			int aip1 = IA[i+1];
			double dii = 0;
			for(int j_=ai,j=JA[j_];j_<aip1 && j<=i; j_++,j=JA[j_])
			{
				if(j<i)
				{
					double eij = AA[j_];
					double fji = get_val(j,i);
					double djj = D[j];
					dii -= eij*fji/djj;
				}
				else
					dii += AA[j_];
			}
			D[i]=dii; }
	}
	// A=E+D0+F, M=(E+D)(1/D)(F+D), match the row sums
	// Dii = D0ii - sum_{j=0}^{i-1}Eij/Djj*sum_{k=j+1}^{n-1}Fjk
	void calculate_MILU( VECTOR& D )
	{
		for(int i=0;i<N;i++) {
			int ai   = IA[i  ];
			int aip1 = IA[i+1];
			double dii = 0;
			for(int j_=ai,j=JA[j_];j_<aip1 && j<=i; j_++,j=JA[j_])
			{
				if(j<i)
				{
					double eij = AA[j_];
					double djj = D[j];
					double sum_fjk=0;
					double fji      ;
					int aj   = IA[j  ];
					int ajp1 = IA[j+1];
					for(int k_=aj,k=JA[k_];k_<ajp1; k_++,k=JA[k_])
						if     (k==i)     fji  = AA[k_];
						else if(k >j) sum_fjk += AA[k_];

					dii -= eij/djj*(fji+.95*sum_fjk); // interpolated MILU
				}
				else
					dii += AA[j_];
			}
			D[i]=dii; }
	}
	// M=(E+D)(1/D)(F+D) = (E/D+I)(F+D)
	void apply_M( const VECTOR& D, VECTOR& X )
	{
		      double* pX = (      double*)X; // 0-based
		const double* pD = (const double*)D;
        		
        // (E/D+I)^{-1} * B
		for(int i=0;i<N;i++){
            double sum = pX[i];
            int ai   = IA[i  ];
            int aip1 = IA[i+1];
            for(int j_=ai,j=JA[j_]; j_<aip1 && j<i;j_++,j=JA[j_])
                sum -= AA[j_]*pX[j]/pD[j];
            pX[i] = sum;}
        
        // (F+D)^{-1} * x
        for(int i=N-1;i>=0;i--){
            double sum = pX[i];
            int ai   = IA[i  ];
            int aip1 = IA[i+1];
            
            for(int j_=aip1-1,j=JA[j_]; j_>=ai && j>i;j_--,j=JA[j_])
                sum -= AA[j_]*pX[j];
            pX[i] = sum/pD[i];}
	}

	void projection_average_zero(VECTOR& X )
	{
		double* pX = (double*)X; 
		double avg =0; for(int i=N;i>0;i--,pX++) avg+=(*pX); pX--; 
		       avg/=N; for(int i=N;i>0;i--,pX--) (*pX)-=avg;
	}

public:
	//
	void PCG( const VECTOR& B,
                    VECTOR& X, Preconditioner prec, 
					           bool   with_initial_guess = false,
							   bool   nonzero_kernel = true,
                               int    max_it = 2000, double residual_norm = 1E-6 )
    {
        assert(N==B.size);
		// preconditioner
		VECTOR D(N); 
		     if(prec== ILU) calculate_ILU (D);
		else if(prec==MILU) calculate_MILU(D);
		else                assert(false);
		// 
		if(X.size!=N){ X.resize(N); X=0; with_initial_guess=false; }
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
		if(nonzero_kernel) projection_average_zero(R);
		memcpy(pZ,pR,sizeof(double)*N); apply_M(D,Z); // Z=invert(M)*R
		memcpy(pP,pZ,sizeof(double)*N);               // P=Z;  
		double  rz = inner_product( R,Z); 
		//---------------------------------------------------------------------
		// iteration
		//---------------------------------------------------------------------
		int iteration = 0; 
		while( iteration++ < max_it ) {
			// a
			apply_A(P,AP); if(nonzero_kernel) projection_average_zero(AP);
			double app = inner_product(AP,P);
			double a   = rz/app;
			// check convergence
			if(ABS(rz)<residual_norm) {
				if(prec== ILU) printf( "ICCG  converged in %d iterations\n", iteration );
				if(prec==MILU) printf( "MICCG converged in %d iterations\n", iteration );
				return;
			}
			// x & r
			for(int i=N;i>0;i--,pX++,pP++,pR++,pAP++){
				*pX+= a*(*pP );
				*pR-= a*(*pAP);
			}
			pX-=N; pP -=N; pR-=N; pAP-=N;
			// b
			memcpy(pZ,pR,sizeof(double)*N);
			apply_M(D,Z);
			double rz_new = inner_product(R,Z);
			double b      = rz_new/rz; rz=rz_new;
			for(int i=N;i>0;i--,pZ++,pP++)
				*pP = (*pZ)+b*(*pP);
			pP-=N; pZ-=N;
		}
		if(prec== ILU) printf( "ICCG  failed to converge in maximum iteration %d\n", iteration );
		if(prec==MILU) printf( "MICCG failed to converge in maximum iteration %d\n", iteration );
	}

	//---------------------------------------------------------------------
	// print
	//---------------------------------------------------------------------
	void print(const char* name) const {
		printf("-------------------SPARSEMATRIX : %s--------------\n",name);
		for(int i=0;i<N;i++) {
			printf("i:%3d | ",i);
			int index_i  =IA[i  ];
			int index_ip1=IA[i+1];
			for(int a=index_i;a<index_ip1;a++){
				int    j = JA[a];
				double v = AA[a];
				printf("(%3d,% 3.5f), ", j,v);}
			printf("\n");}}

	void print_matlab(const char* filename) const {
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
};
#endif
