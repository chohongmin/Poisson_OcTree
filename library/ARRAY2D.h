#ifndef ARRAY2D_H__
#define ARRAY2D_H__

#include <stdio.h>
#include <assert.h>
#include <memory.h>
#include <iostream>

template<class T>
struct ARRAY2D
{
	int isize; 
	int jsize;
	T**  data_2d;
	T*   data;
	//-------------------------------------------------------------
	// memory related functions
	//-------------------------------------------------------------
	ARRAY2D    (int isize=1,int jsize=1);
	void resize(int isize  ,int jsize  );
	ARRAY2D       (const ARRAY2D& V);
	void operator=(const ARRAY2D& V);
	~ARRAY2D();
	//-------------------------------------------------------------
	// array operations
	//-------------------------------------------------------------
	void operator= (T s);
	void operator*=(T s);
	void operator/=(T s);
	void operator+=(T s);
	void operator-=(T s);
	void operator+=(const ARRAY2D& V);
	void operator-=(const ARRAY2D& V);
	//-------------------------------------------------------------------------
	// accessors
	//-------------------------------------------------------------------------
	operator const T*()const{ return data;}
	operator       T*()     { return data;}
	inline virtual const T* operator[](int i) const{ return data_2d[i]; }
	inline virtual       T* operator[](int i)      { return data_2d[i]; }
	inline T iph   (int i,int j)const{ return .5*(data_2d[i][j] + data_2d[i+1][j]);}
	inline T jph   (int i,int j)const{ return .5*(data_2d[i][j] + data_2d[i][j+1]);}
	inline T iphjph(int i,int j)const{ return .5*(iph(i,j) + iph(i,j+1));}
	inline T Dx       (int i,int j,double dx)const{ return (data_2d[i+1][j]-data_2d[i-1][j])/(2*dx);}
	inline T Dy       (int i,int j,double dy)const{ return (data_2d[i][j+1]-data_2d[i][j-1])/(2*dy);}
	inline T Dx_iph   (int i,int j,double dx)const{ return (data_2d[i+1][j]-data_2d[i][j])/dx;}
	inline T Dy_jph   (int i,int j,double dy)const{ return (data_2d[i][j+1]-data_2d[i][j])/dy;}
//	inline T Dx_jph   (int i,int j,double dx)const{ return (iphjph(i,j)-iphjph(i-1,j))/dx;}
//	inline T Dy_iph   (int i,int j,double dy)const{ return (iphjph(i,j)-iphjph(i,j-1))/dy;}
	void print( const char* name, int step=1 ) const;
	void save ( const char* filename=0 ) const;
	void load ( const char* filename   )      ;
	T    max_norm() const;
};

typedef ARRAY2D<double> A2D;

//
template<class T>
ARRAY2D<T>::ARRAY2D(int isize_, int jsize_)
{
	isize=isize_;
	jsize=jsize_;
	data    = new T [isize*jsize];
	data_2d = new T*[isize];
	for(int i=0;i<isize;i++) data_2d[i]=data+i*jsize;
}
//
template<class T>
void ARRAY2D<T>::resize(int isize_, int jsize_)
{
	if(isize!=isize_ || jsize!=jsize_)
	{
		delete[] data; delete[] data_2d;
		isize=isize_;
		jsize=jsize_;
		data    = new T [isize*jsize];
		data_2d = new T*[isize];
		for(int i=0;i<isize;i++) data_2d[i]=data+i*jsize;
	}
}
//
template<class T>
ARRAY2D<T>::ARRAY2D(const ARRAY2D& V)
{	isize=V.isize;
	jsize=V.jsize;
	data    = new T [isize*jsize];
	data_2d = new T*[isize];
	for(int i=0;i<isize;i++) data_2d[i]=data+i*jsize;
	     if(sizeof(T)==sizeof(int   )) memcpy(data,V.data,sizeof(int   )*isize*jsize);
	else if(sizeof(T)==sizeof(double)) memcpy(data,V.data,sizeof(double)*isize*jsize);
	else
	{
		T*       dest =   data;
		const T* sour = V.data;
		for(int n=isize*jsize;n>0;n--,dest++,sour++) *dest=*sour;
	}
}
//
template<class T>
ARRAY2D<T>::~ARRAY2D(){ delete[] data; delete[] data_2d; }
//
template<class T>
void ARRAY2D<T>::operator=(const ARRAY2D& V)
{	
	resize(V.isize,V.jsize);
	     if(sizeof(T)==sizeof(int   )) memcpy(data,V.data,sizeof(int   )*isize*jsize);
	else if(sizeof(T)==sizeof(double)) memcpy(data,V.data,sizeof(double)*isize*jsize);
	else
	{
		T*       dest =   data;
		const T* sour = V.data;
		for(int n=isize*jsize;n>0;n--,dest++,sour++) *dest=*sour;
	}
}
template<class T>void ARRAY2D<T>::operator= (T s){T* p=data;for(int n=isize*jsize;n>0;n--,p++) *p =s;}
template<class T>void ARRAY2D<T>::operator*=(T s){T* p=data;for(int n=isize*jsize;n>0;n--,p++) *p*=s;}
template<class T>void ARRAY2D<T>::operator/=(T s){T* p=data;for(int n=isize*jsize;n>0;n--,p++) *p/=s;}
template<class T>void ARRAY2D<T>::operator+=(T s){T* p=data;for(int n=isize*jsize;n>0;n--,p++) *p+=s;}
template<class T>void ARRAY2D<T>::operator-=(T s){T* p=data;for(int n=isize*jsize;n>0;n--,p++) *p-=s;}
template<class T>void ARRAY2D<T>::operator+=(const ARRAY2D& V){T* p=data;const T* q=V.data;for(int n=isize*jsize;n>0;n--,p++,q++) *p+=*q;}
template<class T>void ARRAY2D<T>::operator-=(const ARRAY2D& V){T* p=data;const T* q=V.data;for(int n=isize*jsize;n>0;n--,p++,q++) *p-=*q;}

template<class T>
void ARRAY2D<T>::print(const char* name, int step) const{
	printf("--------------------%s--------------------\n", name);
	for(int j=jsize-1;j>=0;j-=step) {
		for(int i=0;i<isize;i+=step) {
			if(sizeof(T)==sizeof(double)) printf("% .3f", data_2d[i][j]);
			if(sizeof(T)==sizeof(int   )) printf("% d " , data_2d[i][j]); }
		printf("\n"); }
	printf("\n"); }
template<class T>
T ARRAY2D<T>::max_norm() const{
	T max=0;
	T* p=data;
	for(int n=isize*jsize;n>0;n--,p++) 
	{
		T  datan = *p;
		if(datan<  0) datan=-datan;
		if(datan>max)   max= datan;
	}
	return max;
}
//-------------------------------------------------------------------------
// File I/O
//-------------------------------------------------------------------------
template<class T>
void ARRAY2D<T>::save( const char* name ) const{
	char filename[100]; static int count = 0;
	if(name==0) sprintf( filename, "function_%04d.df2", count++ );
	else        sprintf( filename, "%s", name );
	std::cout << "saving " << filename << std::endl;
	std::ofstream os(filename); os<<isize<<' '<<jsize<<std::endl;
	for(int i=0;i<isize;i++)
	for(int j=0;j<jsize;j++) os<<data_2d[i][j]<<std::endl;}

template<class T>
void ARRAY2D<T>::load( const char* filename ) {
	std::ifstream is(filename); int isize_; is>>isize_;
	                            int jsize_; is>>jsize_; resize(isize_,jsize_);
	for(int i=0;i<isize;i++)
		for(int j=0;j<jsize;j++) is>>data_2d[i][j];}


//-------------------------------------------------------------------------
//
//-------------------------------------------------------------------------
static double get_f_cons_ext   ( const ARRAY2D<double>& F, int i, int j )
{
	i=MAX(i,0); i=MIN(i,F.isize-1);
	j=MAX(j,0); j=MIN(j,F.jsize-1);
	return F[i][j];
}
static double get_fiph_quad   ( const ARRAY2D<double>& F, int i, int j )
{
	double f0 = get_f_cons_ext(F,i-1,j);
	double f1 = get_f_cons_ext(F,i  ,j), f01=f1-f0;
	double f2 = get_f_cons_ext(F,i+1,j), f12=f2-f1, f012=f12-f01;
	double f3 = get_f_cons_ext(F,i+2,j), f23=f3-f2, f123=f23-f12;
	return .5*(f1+f2)-.125*MINMOD(f012,f123);
}
static double get_fjph_quad   ( const ARRAY2D<double>& F, int i, int j )
{
	double f0 = get_f_cons_ext(F,i,j-1);
	double f1 = get_f_cons_ext(F,i,j  ), f01=f1-f0;
	double f2 = get_f_cons_ext(F,i,j+1), f12=f2-f1, f012=f12-f01;
	double f3 = get_f_cons_ext(F,i,j+2), f23=f3-f2, f123=f23-f12;
	return .5*(f1+f2)-.125*MINMOD(f012,f123);
}
static double get_fiphjph_quad( const ARRAY2D<double>& F, int i, int j )
{
	double fxx_1;
	{ double f0 = get_f_cons_ext(F,i-1,j  );
	  double f1 = get_f_cons_ext(F,i  ,j  ), f01=f1-f0;
	  double f2 = get_f_cons_ext(F,i+1,j  ), f12=f2-f1, f012=f12-f01;
	  double f3 = get_f_cons_ext(F,i+2,j  ), f23=f3-f2, f123=f23-f12; fxx_1=MINMOD(f012,f123);
	}
	double fxx_2;
	{ double f0 = get_f_cons_ext(F,i-1,j+1);
	  double f1 = get_f_cons_ext(F,i  ,j+1), f01=f1-f0;
	  double f2 = get_f_cons_ext(F,i+1,j+1), f12=f2-f1, f012=f12-f01;
	  double f3 = get_f_cons_ext(F,i+2,j+1), f23=f3-f2, f123=f23-f12; fxx_2=MINMOD(f012,f123);
	}
	double fyy_1;
	{ double f0 = get_f_cons_ext(F,i  ,j-1);
	  double f1 = get_f_cons_ext(F,i  ,j  ), f01=f1-f0;
	  double f2 = get_f_cons_ext(F,i  ,j+1), f12=f2-f1, f012=f12-f01;
	  double f3 = get_f_cons_ext(F,i  ,j+2), f23=f3-f2, f123=f23-f12; fyy_1=MINMOD(f012,f123);
	}
	double fyy_2;
	{ double f0 = get_f_cons_ext(F,i+1,j-1);
	  double f1 = get_f_cons_ext(F,i+1,j  ), f01=f1-f0;
	  double f2 = get_f_cons_ext(F,i+1,j+1), f12=f2-f1, f012=f12-f01;
	  double f3 = get_f_cons_ext(F,i+1,j+2), f23=f3-f2, f123=f23-f12; fyy_2=MINMOD(f012,f123);
	}
	double f00 = get_f_cons_ext(F,i  ,j  );
	double f01 = get_f_cons_ext(F,i  ,j+1);
	double f10 = get_f_cons_ext(F,i+1,j  );
	double f11 = get_f_cons_ext(F,i+1,j+1);
	return (f00+f01+f10+f11)*.25 - .125*.5*(fxx_1+fxx_2)
		                         - .125*.5*(fyy_1+fyy_2);
}

#endif

