#ifndef _ICCG_H
#define _ICCG_H

#include "SPARSEMATRIX.h"

//---------------------------------------------------------------------
// The Incomplete Cholesky decomposition of the sparsematrix
//---------------------------------------------------------------------
class IC : public SPARSEMATRIX
{
public:
    void create_preconditioner( const SPARSEMATRIX& A )
    {
        // copy the elements of A to this elements
        clear_all();
        Nmax=N=A.number_of_rows    (); IA=new int   [N+1]; memcpy((void*)IA,(const void*)A.get_IA(),(N+1)*sizeof(int   ));
        Mmax=M=A.number_of_elements(); JA=new int   [M  ]; memcpy((void*)JA,(const void*)A.get_JA(),(M  )*sizeof(int   ));
                                       AA=new double[M  ]; memcpy((void*)AA,(const void*)A.get_AA(),(M  )*sizeof(double));
    
        // IC decomposition : IKJ version
		for(int i=0;i<N;i++)
		{
			int ai   = IA[i  ];
			int aip1 = IA[i+1];

			double aii;
					
			for(int k_=ai; k_<aip1; k_++)
			{
				int k = JA[k_];
				
				// aik /= akk
				if(k<i)
				{
					double aik = AA[k_]/get_val(k,k);
					AA[k_] = aik;
					
					// aij -= aik*akj
					for(int j_=k_+1;j_<aip1;j_++)
					{
						int j = JA[j_];
						AA[j_]-=aik*get_val(k,j);
					}
				}
				// aik /= sqrt(aii)
				else
				{
					if(k==i) aii = sqrt(AA[k_]);
					if(aii<1E-12 && aii>-1E-12)
						assert(false);
					AA[k_]/=aii;
				}
			}
		}
	}
    
    // return U_inverse * L_inverse * x
    void inverse( VECTOR& X ) const
    {
        double* pX = (double*) X; // 0-based
        		
        // L_inverse * x
		for(int i=0;i<N;i++)
        {
            double sum = pX[i];
            
            int index_i   = IA[i  ];
            int index_ip1 = IA[i+1];
            
            for(int j_=index_i; j_<index_ip1;j_++)
            {
                int j = JA[j_];
                
                     if(j< i) sum -= AA[j_]*pX[j];
                else if(j==i) sum /= AA[j_];
                else          break;
            }
            
            pX[i] = sum;
        }
        
        // U_inverse * x
        for(int i=N-1;i>=0;i--)
        {
            double sum = pX[i];
            
            int index_i   = IA[i  ];
            int index_ip1 = IA[i+1];
            
            for(int j_=index_ip1-1; j_>=index_i;j_--)
            {
                int j = JA[j_];
                
                     if(j> i) sum -= AA[j_]*pX[j];
                else if(j==i) sum /= AA[j_];
                else          break;
            }
            
            pX[i] = sum;
        }
    }
};


class ICCG
{
public:
    static void solve( const SPARSEMATRIX  & A,
                       const VECTOR& B,
                             VECTOR& X, bool   with_initial_guess = false,
                                                int    max_it             = 2000,
                                                double residual_norm      = 1E-10 )
    {
        //---------------------------------------------------------------------
		// memory allocation
		//---------------------------------------------------------------------
		assert(A.dim()==B.size);
		//StopWatch sw; sw.start();
        IC M; M.create_preconditioner(A);
		//M.print("M");
		//sw.stop(); sw.read_duration("create IC preconditioner");
		//sw.start(); 
     	int N = B.size;
		
		if(X.size!=N){ X.resize(N); X=0; }
		
						const double* pB  = (const double*)B;
		                      double* pX  = (      double*)X;
		VECTOR  R(N); double* pR  = (      double*)R;
		VECTOR  Z(N); double* pZ  = (      double*)Z;
		VECTOR  P(N); double* pP  = (      double*)P;
		VECTOR AP(N); double* pAP = (      double*)AP;
		
		//---------------------------------------------------------------------
		// initilialize R,P,R0
		//---------------------------------------------------------------------
		if(with_initial_guess==false) memset((void*)pX,0,sizeof(double)*N);
		
		A.transform(X,R); for(int i=0;i<N;i++) pR[i] = pB[i] - pR[i];  // R=B-AX
		memcpy(pZ,pR,sizeof(double)*N);M.inverse(Z);//Z=R;  // Z=invert(M)*R
		memcpy(pP,pZ,sizeof(double)*N);//P=Z;  // P=Z
		
		//---------------------------------------------------------------------
		// iteration
		//---------------------------------------------------------------------
		int iteration = 0; 
		while( iteration++ < max_it )
		{
			// a
			A.transform(P,AP); 
			double  rz = inner_product( R,Z); 
			double app = inner_product(AP,P);
			double a   = rz/app;
			
			// check convergence
			if(sqrt(rz)<residual_norm)
			{
				printf( "ICCG converged in %d iterations\n", iteration );
				//sw.stop(); sw.read_duration("ICCG");
				return;
			}
			
			// x & r
			for(int i=0; i<N; i++ ) pX[i]+= a* pP[i];
			for(int i=0; i<N; i++ ) pR[i]-= a*pAP[i];
			// b
			memcpy(pZ,pR,sizeof(double)*N);
			M.inverse(Z);
			double rz_new = inner_product(R,Z);
			double b      = rz_new/rz;
			for(int i=0; i<N; i++ ) pP[i] = pZ[i] + b*pP[i];
		}
		printf( "ICCG failed to converge in maximum iteration %d\n", iteration );
	}
};


/*---------------------------------------------------------------------------
 * solve AX=B with constraint 1*X=0
 *--------------------------------------------------------------------------*/
class ICCG_aug
{
public:
	// projection P=I-11^T/(1^T1)
	static void Project( VECTOR& X ) {
		double avg=0;
		for(int i=0;i<X.size;i++) avg+=X[i]; avg/=X.size;
		for(int i=0;i<X.size;i++) X[i]-=avg;}
	// ICCG : PAPX=PB
	static void solve( const SPARSEMATRIX  & A,
                       const VECTOR& B,
                             VECTOR& X, bool   with_initial_guess = false,
                                                int    max_it             = 2000,
                                                double residual_norm      = 1E-10 )
    {
        //---------------------------------------------------------------------
		// memory allocation
		//---------------------------------------------------------------------
		assert(A.dim()==B.size);
		//StopWatch sw; sw.start();
        IC M; M.create_preconditioner(A);
		//M.print("M");
		//sw.stop(); sw.read_duration("create IC preconditioner");
		//sw.start(); 
     	int N = B.size;
		
		if(X.size!=N){ X.resize(N); X=0; }
		
						const double* pB  = (const double*)B;
		                      double* pX  = (      double*)X;
		VECTOR  R(N); double* pR  = (      double*)R;
		VECTOR  Z(N); double* pZ  = (      double*)Z;
		VECTOR  P(N); double* pP  = (      double*)P;
		VECTOR AP(N); double* pAP = (      double*)AP;
		
		//---------------------------------------------------------------------
		// initilialize R,P,R0
		//---------------------------------------------------------------------
		if(with_initial_guess==false) memset((void*)pX,0,sizeof(double)*N);
		
		Project(X);A.transform(X,R);
		for(int i=0;i<N;i++) pR[i] = pB[i] - pR[i];  // R=B-APX
		Project(R);                                  // R=PB-PAPX
		memcpy(pZ,pR,sizeof(double)*N);M.inverse(Z);//Z=R;  // Z=invert(M)*R
		memcpy(pP,pZ,sizeof(double)*N);//P=Z;  // P=Z
		
		//---------------------------------------------------------------------
		// iteration
		//---------------------------------------------------------------------
		int iteration = 0; 
		while( iteration++ < max_it )
		{
			// a
			Project(P); A.transform(P,AP); Project(AP);
			double  rz = inner_product( R,Z); 
			double app = inner_product(AP,P);
			double a   = rz/app;
			
			// check convergence
			if(sqrt(rz)<residual_norm)
			{
				printf( "ICCG_aug converged in %d iterations\n", iteration );
				//sw.stop(); sw.read_duration("ICCG");
				return;
			}
			
			// x & r
			for(int i=0; i<N; i++ ) pX[i]+= a* pP[i];
			for(int i=0; i<N; i++ ) pR[i]-= a*pAP[i];
			// b
			memcpy(pZ,pR,sizeof(double)*N);
			M.inverse(Z);
			double rz_new = inner_product(R,Z);
			double b      = rz_new/rz;
			for(int i=0; i<N; i++ ) pP[i] = pZ[i] + b*pP[i];
		}
		printf( "ICCG_aug failed to converge in maximum iteration %d\n", iteration );
	}
};

#endif
