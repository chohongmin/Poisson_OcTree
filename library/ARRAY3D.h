#ifndef ARRAY3D_H__
#define ARRAY3D_H__

#include <stdio.h>
#include <assert.h>
#include <memory.h>
#include <iostream>
#include <fstream>

template<class T>
struct ARRAY3D
{
	int isize; 
	int jsize;
	int ksize;
	T*** data_3d;
	T*   data   ;
	//-------------------------------------------------------------
	// memory related functions
	//-------------------------------------------------------------
	ARRAY3D    (int isize=1,int jsize=1,int ksize=1);
	void resize(int isize  ,int jsize  ,int ksize  );
	ARRAY3D       (const ARRAY3D& V);
	void operator=(const ARRAY3D& V);
	~ARRAY3D();
	//-------------------------------------------------------------
	// array operations
	//-------------------------------------------------------------
	void operator= (T s);
	void operator*=(T s);
	void operator/=(T s);
	void operator+=(T s);
	void operator-=(T s);
	void operator+=(const ARRAY3D& V);
	void operator-=(const ARRAY3D& V);
	//-------------------------------------------------------------------------
	// accessors
	//-------------------------------------------------------------------------
	operator const T*()const{ return data;}
	operator       T*()     { return data;}
	inline virtual T*const* operator[](int i) const{ return data_3d[i]; }
	inline virtual      T** operator[](int i)      { return data_3d[i]; }
	inline T iph      (int i,int j,int k )const{ return .5*(data_3d[i][j][k] + data_3d[i+1][j][k]);}
	inline T jph      (int i,int j,int k )const{ return .5*(data_3d[i][j][k] + data_3d[i][j+1][k]);}
	inline T kph      (int i,int j,int k )const{ return .5*(data_3d[i][j][k] + data_3d[i][j][k+1]);}
	inline T jphkph   (int i,int j,int k )const{ return .5*(jph(i,j,k) + jph(i,j,k+1));}
	inline T iphkph   (int i,int j,int k )const{ return .5*(kph(i,j,k) + kph(i+1,j,k));}
	inline T iphjph   (int i,int j,int k )const{ return .5*(iph(i,j,k) + iph(i,j+1,k));}
	inline T iphjphkph(int i,int j,int k )const{ return .5*(iphjph(i,j,k) + iphjph(i,j,k+1));}
	inline T Dx       (int i,int j,int k,double dx)const{ return (data_3d[i+1][j][k]-data_3d[i-1][j][k])/(2*dx);}
	inline T Dy       (int i,int j,int k,double dy)const{ return (data_3d[i][j+1][k]-data_3d[i][j-1][k])/(2*dy);}
	inline T Dz       (int i,int j,int k,double dz)const{ return (data_3d[i][j][k+1]-data_3d[i][j][k-1])/(2*dz);}
	inline T Dx_iph   (int i,int j,int k,double dx)const{ return (data_3d[i+1][j][k]-data_3d[i][j][k])/dx;}
	inline T Dy_jph   (int i,int j,int k,double dy)const{ return (data_3d[i][j+1][k]-data_3d[i][j][k])/dy;}
	inline T Dz_kph   (int i,int j,int k,double dz)const{ return (data_3d[i][j][k+1]-data_3d[i][j][k])/dz;}
	//inline T Dx_jph   (int i,int j,int k,double dx)const{ return (iphjph(i,j,k)-iphjph(i-1,j,k))/dx;}
	//inline T Dx_kph   (int i,int j,int k,double dx)const{ return (iphkph(i,j,k)-iphkph(i-1,j,k))/dx;}
	//inline T Dy_iph   (int i,int j,int k,double dy)const{ return (iphjph(i,j,k)-iphjph(i,j-1,k))/dy;}
	//inline T Dy_kph   (int i,int j,int k,double dy)const{ return (jphkph(i,j,k)-jphkph(i,j-1,k))/dy;}
	//inline T Dz_iph   (int i,int j,int k,double dz)const{ return (iphkph(i,j,k)-iphkph(i,j,k-1))/dz;}
	//inline T Dz_jph   (int i,int j,int k,double dz)const{ return (jphkph(i,j,k)-jphkph(i,j,k-1))/dz;}
	void print_j(const char* name, int j, int step=1) const;
	void print_k(const char* name, int k, int step=1) const;
	void save ( const char* filename=0 ) const;
	void load ( const char* filename   )      ;
	T    max_norm() const;
};

typedef ARRAY3D<double> A3D;


//
template<class T>
ARRAY3D<T>::ARRAY3D(int isize_, int jsize_, int ksize_)
{
	isize=isize_;
	jsize=jsize_;
	ksize=ksize_;
	                          data         = new T  [isize*jsize*ksize];
	                          data_3d      = new T**[isize];
	for(int i=0;i<isize;i++){ data_3d[i]   = new T* [jsize];
	for(int j=0;j<jsize;j++)  data_3d[i][j]=data+(i*jsize+j)*ksize;}
}
//
template<class T>
void ARRAY3D<T>::resize(int isize_, int jsize_, int ksize_)
{
	if( isize!=isize_ || 
		jsize!=jsize_ ||
		ksize!=ksize_ )
	{
		for(int i=0;i<isize;i++) delete[] data_3d[i];
		                         delete[] data_3d   ;
								 delete[] data      ;
		isize=isize_;
		jsize=jsize_;
		ksize=ksize_;
		                          data         = new T  [isize*jsize*ksize];
								  data_3d      = new T**[isize];
		for(int i=0;i<isize;i++){ data_3d[i]   = new T* [jsize];
		for(int j=0;j<jsize;j++)  data_3d[i][j]=data+(i*jsize+j)*ksize;}
	}
}
//
template<class T>
ARRAY3D<T>::ARRAY3D(const ARRAY3D& V)
{	isize=V.isize;
	jsize=V.jsize;
	ksize=V.ksize;
	                          data         = new T  [isize*jsize*ksize];
	                          data_3d      = new T**[isize];
	for(int i=0;i<isize;i++){ data_3d[i]   = new T* [jsize];
	for(int j=0;j<jsize;j++)  data_3d[i][j]=data+(i*jsize+j)*ksize;}
	     if(sizeof(T)==sizeof(int   )) memcpy(data,V.data,sizeof(int   )*isize*jsize*ksize);
	else if(sizeof(T)==sizeof(double)) memcpy(data,V.data,sizeof(double)*isize*jsize*ksize);
	else
	{
		T*       dest =   data;
		const T* sour = V.data;
		for(int n=isize*jsize*ksize;n>0;n--,dest++,sour++) *dest=*sour;
	}
}
//
template<class T>
ARRAY3D<T>::~ARRAY3D(){
	for(int i=0;i<isize;i++) delete[] data_3d[i];
		                     delete[] data_3d   ;
							 delete[] data      ;
}
//
template<class T>
void ARRAY3D<T>::operator=(const ARRAY3D& V)
{	
	resize(V.isize,V.jsize,V.ksize);
	     if(sizeof(T)==sizeof(int   )) memcpy(data,V.data,sizeof(int   )*isize*jsize*ksize);
	else if(sizeof(T)==sizeof(double)) memcpy(data,V.data,sizeof(double)*isize*jsize*ksize);
	else
	{
		T*       dest =   data;
		const T* sour = V.data;
		for(int n=isize*jsize*ksize;n>0;n--,dest++,sour++) *dest=*sour;
	}
}
template<class T>void ARRAY3D<T>::operator= (T s){T* p=data;for(int n=isize*jsize*ksize;n>0;n--,p++) *p =s;}
template<class T>void ARRAY3D<T>::operator*=(T s){T* p=data;for(int n=isize*jsize*ksize;n>0;n--,p++) *p*=s;}
template<class T>void ARRAY3D<T>::operator/=(T s){T* p=data;for(int n=isize*jsize*ksize;n>0;n--,p++) *p/=s;}
template<class T>void ARRAY3D<T>::operator+=(T s){T* p=data;for(int n=isize*jsize*ksize;n>0;n--,p++) *p+=s;}
template<class T>void ARRAY3D<T>::operator-=(T s){T* p=data;for(int n=isize*jsize*ksize;n>0;n--,p++) *p-=s;}
template<class T>void ARRAY3D<T>::operator+=(const ARRAY3D& V){T* p=data;const T* q=V.data;for(int n=isize*jsize*ksize;n>0;n--,p++,q++) *p+=*q;}
template<class T>void ARRAY3D<T>::operator-=(const ARRAY3D& V){T* p=data;const T* q=V.data;for(int n=isize*jsize*ksize;n>0;n--,p++,q++) *p-=*q;}

//
template<class T>
void ARRAY3D<T>::print_j(const char* name, int j, int step) const {
printf("--------------------%s :(j=%d/%d)--------------------\n", name,j,jsize);
for(int k=ksize-1;k>=0;k-=step) {
	for(int i=0;i<isize;i+=step) {
		if(sizeof(T)==sizeof(double)) printf("% .3f", data_3d[i][j][k]);
		if(sizeof(T)==sizeof(int   )) printf("% d " , data_3d[i][j][k]); }
	printf("\n"); }
printf("\n"); }
//
template<class T>
void ARRAY3D<T>::print_k(const char* name, int k, int step) const {
	printf("--------------------%s :(k=%d/%d)--------------------\n", name,k,ksize);
	for(int j=jsize-1;j>=0;j-=step) {
		for(int i=0;i<isize;i+=step) {
			if(sizeof(T)==sizeof(double)) printf("% .3f", data_3d[i][j][k]);
			if(sizeof(T)==sizeof(int   )) printf("% d " , data_3d[i][j][k]); }
		printf("\n"); }
	printf("\n"); }
//
template<class T>
T ARRAY3D<T>::max_norm() const{
	T max=0;
	T* p=data;
	for(int n=isize*jsize*ksize;n>0;n--,p++) 
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
void ARRAY3D<T>::save( const char* name ) const{
	char filename[100]; static int count = 0;
		if(name==0) sprintf( filename, "function_%04d.df3", count++ );
		else        sprintf( filename, "%s", name );
		std::cout << "saving " << filename << std::endl;
		std::ofstream os(filename); os<<isize<<' '
			                          <<jsize<<' '
									  <<ksize<<std::endl;
		for(int i=0;i<isize;i++)
		for(int j=0;j<jsize;j++)
		for(int k=0;k<ksize;k++) os<<data_3d[i][j][k]<<std::endl;}

template<class T>
void ARRAY3D<T>::load( const char* filename ) {
	std::ifstream is(filename); int isize_; is>>isize_;
		                            int jsize_; is>>jsize_;
									int ksize_; is>>ksize_; resize(isize_,jsize_,ksize_);
		for(int i=0;i<isize;i++)
		for(int j=0;j<jsize;j++)
		for(int k=0;k<ksize;k++) is>>data_3d[i][j][k];}
#endif

