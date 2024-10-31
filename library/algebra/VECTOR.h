#ifndef VECTOR_H
#define VECTOR_H

#include <malloc.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include "../ARRAY.h"

//------------------------------------------------------------------------
// VECTOR : 0-based indexing
//------------------------------------------------------------------------
class VECTOR : public ARRAY<double>
{
public:
	VECTOR(int dim=1) : ARRAY<double>(dim){}

	/*-------------------------------------------------------------*\
	 * array operations
	\*-------------------------------------------------------------*/
	void operator= (double s){for(int n=0;n<size;n++)data[n] =s;}
	void operator*=(double s){for(int n=0;n<size;n++)data[n]*=s;}
	void operator/=(double s){for(int n=0;n<size;n++)data[n]/=s;}
	void operator+=(double s){for(int n=0;n<size;n++)data[n]+=s;}
	void operator-=(double s){for(int n=0;n<size;n++)data[n]-=s;}

	void operator+=(const VECTOR& V){for(int n=0;n<size;n++)data[n]+=V.data[n];}
	void operator-=(const VECTOR& V){for(int n=0;n<size;n++)data[n]-=V.data[n];}
	             inline VECTOR operator*(double s                    ){VECTOR copy(*this); copy*=s; return copy;}
	friend       inline VECTOR operator*(double s, const VECTOR& that){VECTOR copy( that); copy*=s; return copy;}
	inline VECTOR operator+(const VECTOR& that){VECTOR copy(*this); copy+=that;return copy;}
	inline VECTOR operator-(const VECTOR& that){VECTOR copy(*this); copy-=that;return copy;}

	operator const double*() const {return data;}
	operator       double*()       {return data;}

	//-------------------------------------------------------------------------
	// accessors
	//-------------------------------------------------------------------------
	double  operator[](int r) const {return data[r];}
	double& operator[](int r)       {return data[r];}
	
	void print(const char* name) const
	{
		printf("--------------------%s--------------------\n", name);
		for(int j=0;j<size;j++)
			printf("% .4f    ", data[j]);
		printf("\n");
	}
};

static double inner_product( const VECTOR& V1, const VECTOR& V2)
{
	assert(V1.size==V2.size);
	int    dim=V1.size;
	double sum=0;
	const double* p1 = (const double*)V1;
	const double* p2 = (const double*)V2;
	while(dim>0)
	{
		sum += (*p1)*(*p2);
		p1++;
		p2++;
		dim--;
	}
	return sum;
}

#endif

