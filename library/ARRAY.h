#ifndef ARRAY_H__
#define ARRAY_H__

#include <stdio.h>
#include <assert.h>
#include <memory.h>
#include <fstream>
#include "Macros.h"

template<class T>
struct ARRAY
{
	int size; 
	T*  data;
	//-------------------------------------------------------------
	// memory related functions
	//-------------------------------------------------------------
	ARRAY      (int size=1);
	void resize(int size  );
	void resize_with_existing_data_kept(int size  );
	ARRAY         (const ARRAY& V);
	void operator=(const ARRAY& V);
	~ARRAY();
	//-------------------------------------------------------------
	// array operations
	//-------------------------------------------------------------
	void operator= (T s);
	void operator*=(double s);
	void operator/=(double s);
	void operator+=(T s);
	void operator-=(T s);
	void operator+=(const ARRAY& V);
	void operator-=(const ARRAY& V);
	//-------------------------------------------------------------------------
	// accessors
	//-------------------------------------------------------------------------
	operator const T*()const{ return data;}
	operator       T*()     { return data;}
	inline virtual T  operator[](int i) const{ return data[i]; }
	inline virtual T& operator[](int i)      { return data[i]; }
	void print( const char* name       ) const;
	void save ( const char* filename=0 ) const;
	void load ( const char* filename   )      ;
	T    max_norm() const;
};



//
template<class T>
ARRAY<T>::ARRAY(int size_){
	size=size_;
	data=new T[size];}
//
template<class T>
void ARRAY<T>::resize(int size_){
	if(size!=size_){
		delete[] data;
		size=size_;
		data=new T[size];}}
//
template<class T>
void ARRAY<T>::resize_with_existing_data_kept(int size_){
	if(size!=size_){
		T* data_new = new T[size_]; int size_copy=MIN(size,size_);
		memcpy(data_new,data,sizeof(T)*size_copy);
		delete[] data;size=size_;data=data_new;}}
//
template<class T>
ARRAY<T>::ARRAY(const ARRAY& V){	size=V.size;
	data=new T[size];
	     if(sizeof(T)==sizeof(int   )) memcpy(data,V.data,sizeof(int   )*size);
	else if(sizeof(T)==sizeof(double)) memcpy(data,V.data,sizeof(double)*size);
	else{
		T*       dest =   data;
		const T* sour = V.data;
		for(int i=size;i>0;i--,dest++,sour++)
			*dest=*sour;}}
//
template<class T>
ARRAY<T>::~ARRAY(){ delete[] data; }
//
template<class T>
void ARRAY<T>::operator=(const ARRAY& V)
{	
	resize(V.size);
	     if(sizeof(T)==sizeof(int   )) memcpy(data,V.data,sizeof(int   )*size);
	else if(sizeof(T)==sizeof(double)) memcpy(data,V.data,sizeof(double)*size);
	else
	{
		T*       dest =   data;
		const T* sour = V.data;
		for(int i=size;i>0;i--,dest++,sour++)
			*dest=*sour;
	}
}
template<class T>void ARRAY<T>::operator= (T s){for(int n=0;n<size;n++)data[n] =s;}
template<class T>void ARRAY<T>::operator*=(double s){for(int n=0;n<size;n++)data[n]*=s;}
template<class T>void ARRAY<T>::operator/=(double s){for(int n=0;n<size;n++)data[n]/=s;}
template<class T>void ARRAY<T>::operator+=(T s){for(int n=0;n<size;n++)data[n]+=s;}
template<class T>void ARRAY<T>::operator-=(T s){for(int n=0;n<size;n++)data[n]-=s;}
template<class T>void ARRAY<T>::operator+=(const ARRAY& V){for(int n=0;n<size;n++)data[n]+=V.data[n];}
template<class T>void ARRAY<T>::operator-=(const ARRAY& V){for(int n=0;n<size;n++)data[n]-=V.data[n];}

template<class T>
void ARRAY<T>::print(const char* name) const {
	printf("----------------------%s--------------------\n", name);
	for(int j=0;j<size;j++) {
		if(sizeof(T)==sizeof(double)) printf("% .4f ", data[j]);
		if(sizeof(T)==sizeof(int   )) printf("% d "  , data[j]); }
	printf("\n");}
template<class T>
T ARRAY<T>::max_norm() const{
	T max=0;
	for(int i=0;i<size;i++)
	{
		T datai = data[i];
		if(datai<0) datai=-datai;
		if(datai>max) max=datai;
	}
	return max;
}
//-------------------------------------------------------------------------
// File I/O
//-------------------------------------------------------------------------
template<class T>
void ARRAY<T>::save( const char* name ) const{
	char filename[100]; static int count = 0;
	if(name==0) sprintf( filename, "function_%04d.df1", count++ );
	else        sprintf( filename, "%s", name );
	std::ofstream os(filename); os<<size<<std::endl;
	for(int i=0;i<size;i++) os<<data[i]<<std::endl;}

template<class T>
void ARRAY<T>::load( const char* filename ) {
	std::ifstream is(filename); int size_; is>>size_; resize(size_);
	for(int i=0;i<size;i++) is>>data[i];}



#endif

