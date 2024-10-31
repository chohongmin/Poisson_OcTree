#ifndef _LAPLACETYPEMATRIX3D_H_
#define _LAPLACETYPEMATRIX3D_H_

class LaplaceTypeMatrix3D {
public:
	ARRAY_3D<double> Adiag, Aiph, Ajph, Akph, precond;
	
	double dotproduct( const ARRAY_3D<double>& X,
		               const ARRAY_3D<double>& Y)
	{
		return 0;
	}

	void applyA( const ARRAY_3D<double>&  X,
		               ARRAY_3D<double>& AX )
	{
		for(int i=1;i<X.isize-1;i++)
		for(int j=1;j<X.jsize-1;j++)
		for(int k=1;k<X.ksize-1;k++)
		{
			AX(i,j,k) = Adiag(i  ,j  ,k  )*X(i  ,j  ,k  )
				      + Aiph (i  ,j  ,k  )*X(i+1,j  ,k  )
					  + Aiph (i-1,j  ,k  )*X(i-1,j  ,k  )
					  + Ajph (i  ,j  ,k  )*X(i  ,j+1,k  )
					  + Ajph (i  ,j-1,k  )*X(i  ,j-1,k  )
					  + Akph (i  ,j  ,k  )*X(i  ,j  ,k+1)
					  + Akph (i  ,j  ,k-1)*X(i  ,j  ,k-1);
		}
	}

	void calculate_MIC0()
	{
		double t=0.97;
		double s=0.25;
		for(int i=0;i<isize;i++)
		for(int j=0;j<jsize;j++)
		for(int k=0;k<ksize;k++)
		{
			double e = Adiag(i,j,k)- SQR(Aiph(i-1,j,k)*precond(i-1,j,k))
				                   - SQR(Ajph(i,j-1,k)*precond(i,j-1,k))
								   - SQR(Akph(i,j,k-1)*precond(i,j,k-1))
				     - t*(Aiph(i-1,j,k)*(Ajph(i-1,j,k)+Akph(i-1,j,k))*SQR(precond(i-1,j,k))
					     +Ajph(i,j-1,k)*(Aiph(i,j-1,k)+Akph(i,j-1,k))*SQR(precond(i,j-1,k))
						 +Akph(i,j,k-1)*(Aiph(i,j,k-1)+Ajph(i,j,k-1))*SQR(precond(i,j,k-1)));
			if(e<s*Adiag(i,j,k)) e=Adiag(i,j,k);
			precond(i,j,k) = 1/SQRT(e);
		}
	}

	void solve_MIC0_PCG( ARRAY_3D<double>& B )
	{
		ARRAY_3D<double>& P=B;
		ARRAY_3D<double> R; R=B;
		ARRAY_3D<double> Z; Z=R; applyPreconditioner(Z);
		ARRAY_3D<double> S; S=Z;
		double sigma = dotproduct(Z,R);
		for(int i=0;i<2000;i++)
		{
			 Z=S; applyA(Z);
			double alpha = sigma/dotproduct(Z,S);
			P+=alpha*P;
			R-=alpha*Z;
		}
	}
};
#endif
