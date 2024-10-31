#ifndef _QUATERNION_ALGEBRA__
#define _QUATERNION_ALGEBRA__

#include <math.h>
#include "../geometry/POINT3.h"
#include "MATRIX.h"

//---------------------------------------------------------------------
// QUATERNION Algebra 
//---------------------------------------------------------------------
struct QUATERNION
{
	double a,u1,u2,u3;

	//---------------------------------------------------------------------------
	// 
	//---------------------------------------------------------------------------
	QUATERNION( const POINT3& P ){ a=0; u1=P.x; u2=P.y; u3=P.z; }
	QUATERNION(){a=1;u1=0;u2=0;u3=0;}
	QUATERNION( double a_, double u1_, double u2_, double u3_ ){a=a_;u1=u1_;u2=u2_;u3=u3_;}

	//---------------------------------------------------------------------------
	//
	//---------------------------------------------------------------------------
	void set_rotation( double theta, double u1_, 
		                             double u2_, 
									 double u3_)
	{
		double abs = sqrt(u1_*u1_ + u2_*u2_ + u3_*u3_);
		
		double cos_theta = cos(theta*.5);
		double sin_theta = sin(theta*.5);

		a = cos_theta;
		u1= sin_theta*u1_/abs;
		u2= sin_theta*u2_/abs;
		u3= sin_theta*u3_/abs;
	}

	//---------------------------------------------------------------------
	// QUATERNION Algebra
	//---------------------------------------------------------------------
	QUATERNION operator+( const QUATERNION& Q)const{ return QUATERNION(a+Q.a,u1+Q.u1,u2+Q.u2,u3+Q.u3);}
	QUATERNION operator-( const QUATERNION& Q)const{ return QUATERNION(a-Q.a,u1-Q.u1,u2-Q.u2,u3-Q.u3);}
	QUATERNION operator*( double s           )const{ return QUATERNION(s*a,s*u1,s*u2,s*u3);}
	QUATERNION operator/( double s           )const{ return QUATERNION(a/s,u1/s,u2/s,u3/s);}
	QUATERNION operator*( const QUATERNION& Q)const
	{
		QUATERNION Q_;
		Q_.a  = a*Q.a  - u1*Q.u1 - u2*Q.u2 - u3*Q.u3;
		Q_.u2 = a*Q.u2 - u1*Q.u3 + u2*Q.a  + u3*Q.u1;
		Q_.u3 = a*Q.u3 + u1*Q.u2 - u2*Q.u1 + u3*Q.a;
		Q_.u1 = a*Q.u1 + u1*Q.a  + u2*Q.u3 - u3*Q.u2;
		return Q_;
	}
	
	QUATERNION friend operator*(const POINT3& P, const QUATERNION& Q)
	{
		QUATERNION Q_;
		Q_.a  = - P.x*Q.u1 - P.y*Q.u2 - P.z*Q.u3;
		Q_.u2 = - P.x*Q.u3 + P.y*Q.a  + P.z*Q.u1;
		Q_.u3 = + P.x*Q.u2 - P.y*Q.u1 + P.z*Q.a;
		Q_.u1 = + P.x*Q.a  + P.y*Q.u3 - P.z*Q.u2;
		return Q_;
	}
	
	
	QUATERNION inverse() const
	{
		double sqr = a*a+u1*u1+u2*u2+u3*u3;

		QUATERNION Q_;
		Q_.a  =   a/sqr;
		Q_.u1 = -u1/sqr;
		Q_.u2 = -u2/sqr;
		Q_.u3 = -u3/sqr;
		return Q_;
	}

	void operator+=(const QUATERNION& Q){a+=Q.a;u1+=Q.u1;u2+=Q.u2;u3+=Q.u3;}
	void operator-=(const QUATERNION& Q){a-=Q.a;u1-=Q.u1;u2-=Q.u2;u3-=Q.u3;}
	void operator*=(double            s){a*=  s;u1*=  s ;u2*=  s ;u3*=  s ;}
	void operator/=(double            s){a/=  s;u1/=  s ;u2/=  s ;u3/=  s ;}
	//---------------------------------------------------------------------
	//
	//---------------------------------------------------------------------
	void normalize()
	{
		double norm = sqrt(a*a+u1*u1+u2*u2+u3*u3);
		a /=norm;
		u1/=norm;
		u2/=norm;
		u3/=norm;
	}
	
	//---------------------------------------------------------------------
	//
	//---------------------------------------------------------------------
	void apply_rotation( POINT3& P ) const
	{
		double& v1 = P.x;
		double& v2 = P.y;
		double& v3 = P.z;
		
		double  b =      - u1*v1 - u2*v2 - u3*v3;
		double w1 = a*v1         + u2*v3 - u3*v2;
		double w2 = a*v2 - u1*v3         + u3*v1;
		double w3 = a*v3 + u1*v2 - u2*v1             ;

		v1 =-b*u1 + w1*a  - w2*u3 + w3*u2;
		v2 =-b*u2 + w1*u3 + w2*a  - w3*u1;
		v3 =-b*u3 - w1*u2 + w2*u1 + w3*a;
	}
	
	//---------------------------------------------------------------------
	//
	//---------------------------------------------------------------------
	MATRIX QUATERNION_to_matrix() const
	{
		MATRIX R(3,3);
		
		R(0,0)=1-2*u2*u2-2*u3*u3; R(0,1)=2*u1*u2-2*a*u3; R(0,2)=2*u1*u3+2*a*u2;
		R(1,1)=1-2*u1*u1-2*u3*u3; R(1,0)=2*u1*u2+2*a*u3; R(1,2)=2*u2*u3-2*a*u1;
		R(2,2)=1-2*u1*u1-2*u2*u2; R(2,0)=2*u1*u3-2*a*u2; R(2,1)=2*u2*u3+2*a*u1;
		
		return R;
	}
};

#endif


