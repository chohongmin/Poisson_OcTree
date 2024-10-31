#ifndef LINEAR_PROGRAMMING_H
#define LINEAR_PROGRAMMING_H

/*--------------------------------------------------------*\
 * maximize C*x where X>=0 and  A   X<=B, or
 * minimize B*Y where Y>=0 and (A^T)Y>=C
 * by the Simplex method with two phases
 * first phase searches a basic feasible solution
 * second phase runs primal or dual simplex on the feasible
 * solution
\*--------------------------------------------------------*/
class LinearProgramming
{
public:
	double** A;
	double*  C;
	double*  B;
	double   v;
	int    m,n;
protected:
	int   *  J;
	int   *  I;
	//
	void pivot(int i, int j)
	{
		// update I and J
		int save=J[j]; J[j]=I[i]; I[i]=save;
		// update Aij
		double Aij=1./(A[i][j]); A[i][j]=Aij;
		// update v
		v+=C[j]*B[i]*Aij;
		// update C and B
		C[j]*=-Aij; for(int k=0;k<n;k++) if(k!=j) C[k]+=C[j]*A[i][k];
		B[i]*= Aij; for(int h=0;h<m;h++) if(h!=i) B[h]-=B[i]*A[h][j];
		// update A[i][k]
		for(int k=0;k<n;k++) if(k!=j) A[i][k]*= Aij;
		// update A[h][k]
		double* Ai=A[i];
		for(int h=0;h<m;h++)
		{
			if(h!=i)
			{
				double* Ah=A[h];
				for(int k=0;k<n;k++)
					if(k!=j)
						Ah[k]-=Ah[j]*Ai[k];
			}
		}
		// updae A[h][j]
		for(int h=0;h<m;h++) if(h!=i) A[h][j]*=-Aij;
	}

	enum STATE
	{
		DONE,
		PRIMAL_INFEASIBLE,
		DUAL_INFEASIBLE,
		STILL_DOING
	};
	//
	STATE primal_simplex()
	{
		// choose j0 with C[j0]>0, "the most positive"
		int j0=-1; double Cj0=0;
		for(int j=0;j<n;j++)
			if(C[j]>Cj0)
			{ j0=j; Cj0=C[j]; }
		// 
		if(j0==-1) return DONE;
		else
		{
			// choose i0 with smallest B[i]/A[i][j0] with A[i][j0]>0
			int i0=-1; double ratio_min=0;
			for(int i=0;i<m;i++)
				if(A[i][j0]>0)
				{
					double ratio = B[i]/A[i][j0];
					if(i0==-1 || ratio<=ratio_min)
					{
						i0=i; ratio_min=ratio;
					}
				}
			//
			if(i0==-1) return DUAL_INFEASIBLE;
			else
			{
				pivot(i0,j0);
				return STILL_DOING;
			}
		}
	}

	STATE dual_simplex()
	{
		// choose i0 with B[i0]<0, "the most negative"
		int i0=-1; double Bi0=0;
		for(int i=0;i<m;i++)
			if(B[i]<Bi0)
			{ i0=i; Bi0=B[i]; }
		// 
		if(i0==-1) return DONE;
		else
		{
			// choose j0 with largest  C[j]/A[i0][j] with A[i0][j]<0
			int j0=-1; double ratio_max=0;
			for(int j=0;j<n;j++)
				if(A[i0][j]<0)
				{
					double ratio = C[j]/A[i0][j];
					if(j0==-1 || ratio>=ratio_max)
					{
						j0=j; ratio_max=ratio;
					}
				}
			//
			if(j0==-1) return PRIMAL_INFEASIBLE;
			else
			{
				pivot(i0,j0);
				return STILL_DOING;
			}
		}
	}

	STATE dual_simplex_first_phase() // pretend C=0
	{
		// choose i0 with B[i0]<0, "the most negative"
		int i0=-1; double Bi0=0;
		for(int i=0;i<m;i++)
			if(B[i]<Bi0)
			{ i0=i; Bi0=B[i]; }
		// 
		if(i0==-1) return DONE;
		else
		{
			// choose j0 with largest  C[j]/A[i0][j] with A[i0][j]<0
			int j0=-1; double ratio_max=0;
			for(int j=0;j<n;j++)
				if(A[i0][j]<0)
				{
					double ratio = 0;
					if(j0==-1 || ratio>=ratio_max)
					{
						j0=j; ratio_max=ratio;
					}
				}
			//
			if(j0==-1) return PRIMAL_INFEASIBLE;
			else
			{
				pivot(i0,j0);
				return STILL_DOING;
			}
		}
	}

	STATE solve()
	{
		bool all_C_negative = true; for(int j=0;j<n;j++) if(C[j]>0){ all_C_negative=false; break; } 
		bool all_B_positive = true; for(int i=0;i<m;i++) if(B[i]<0){ all_B_positive=false; break; }

		STATE st;

		// DONE
		if(all_C_negative && all_B_positive) st=DONE;
		// dual simplex
		else if(all_C_negative)
		{
			do{st=dual_simplex();} while(st==STILL_DOING);
		}
		// primal simplex
		else if(all_B_positive)
		{
			do{st=primal_simplex();} while(st==STILL_DOING);
		}
		// 
		else
		{
			// 1st phase
			do{st=dual_simplex_first_phase();} while(st==STILL_DOING);
			if(st==DONE)
			{
				// 2nd phase
				do{st=primal_simplex();} while(st==STILL_DOING);
			}
		}
		
		return st;
	}

public:
	 LinearProgramming(){m=n=0;A=0;C=0;B=0;J=0;I=0;}
	~LinearProgramming()
	{
		if(m!=0 || n!=0)
		{
			free(C); free(B); free(J); free(I);
			for(int i=0;i<m;i++) free(A[i]); free(A);
		}
	}

	void initialize(int m_, int n_)
	{
		m=m_; n=n_;
		A=(double**)malloc(sizeof(double*)*m);
		C=(double* )malloc(sizeof(double )*n);
		B=(double* )malloc(sizeof(double )*m);
		J=(int   * )malloc(sizeof(int    )*n);
		I=(int   * )malloc(sizeof(int    )*m);
		for(int i=0;i<m;i++)
		{
			A[i]=(double* )malloc(sizeof(double )*n);
			I[i]=n+i;
		}
		for(int j=0;j<n;j++)
			J[j]=j;
		v=0;
	}

	//
	void print_tableau() const
	{
		printf("  |"      ); for(int j=0;j<n;j++) printf("% 2d    ",J[j]   ); printf("|  \n"); for(int i=0;i<m;i++){
		printf("%2d|",I[i]); for(int j=0;j<n;j++) printf("% 02.2f ",A[i][j]); printf("|%02.2f\n",B[i]);}
		printf("  |"      ); for(int j=0;j<n;j++) printf("% 02.2f ",C[j]   ); printf("|%02.2f\n\n",v);
	}

	//
	STATE solve_primal(double* X, double& maximum)
	{
		STATE st=solve();
		for(int j=0;j<n;j++) if(J[j]<n) X[J[j]]=0;
		for(int i=0;i<m;i++) if(I[i]<n) X[I[i]]=B[i];
		maximum=v;
		return st;
	}

	STATE solve_dual(double* Y, double& minimum)
	{
		STATE st=solve();
		for(int j=0;j<n;j++) if(J[j]>=n) Y[J[j]-n]=-C[j];
		for(int i=0;i<m;i++) if(I[i]>=n) Y[I[i]-n]=0;
		minimum=v;
		return st;
	}
};

#endif
