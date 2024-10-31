#ifndef SPARSE_SYMM_MATRIX_H
#define SPARSE_SYMM_MATRIX_H

#include "../ARRAY.h"
#include "VECTOR.h"

//----------------------------------------------------------------------------
// Sparse symmetric matrix : store only upper triangular part
//----------------------------------------------------------------------------
class SparseSymmMatrix
{
protected:
	ARRAY<   int> JA; //0~N  (      IA), N+1~N+M(JA                 )
	ARRAY<double> AA; //0~N-1(diagonal), N+1~N+M(AA:upper triangular)
	VECTOR D;
	int N,M,Mmax,imax;
	int find_index( int i, int j ) const;
public:
	SparseSymmMatrix(int N=0);
	void set_dimension (int N);
	double get_val(int i,int j) const;
	void add_element(int i,int j, double v);
	void apply_A( const double* X, double* Y ) const;
	void apply_M( double* X ) const;
	static double inner_product( const double* X, const double* Y, int N );
	void PCG( const VECTOR& B, VECTOR& X, bool with_initial_guess=false,
						 bool nonzero_kernel    =true ,
                         int  max_it            =2000 , 
						 double tolerance = 1E-6  );
	void print       (const char*     name) const;
	void print_matlab(const char* filename) const;
	void projection_sum_zero( double* X ) const;
	void calculate_MILU() ;
};

#endif
