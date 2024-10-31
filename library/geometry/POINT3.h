#ifndef POINT3_
#define POINT3_

#include "POINT2.h"

//-----------------------------------------------
// 
//  POINT3 class
//
//-----------------------------------------------
class POINT3
{
public:
	double x, y, z;

	double& operator()(int);
	POINT3() ;
	POINT3(double c1, double c2, double c3);
	POINT3 (const POINT2& pt);
	POINT3 operator+( const POINT3& r ) const { POINT3 out; out = *this; out += r; return out; }
	POINT3 operator-( const POINT3& r ) const { POINT3 out; out = *this; out -= r; return out; }
	POINT3 cross    ( const POINT3& r ) const { POINT3 out; cross(*this,r,out);    return out; }
	POINT3 operator*( double s          ) const { POINT3 out; out = *this; out *= s; return out; }
	POINT3 operator/( double s          ) const { POINT3 out; out = *this; out /= s; return out; }

	void operator =( const POINT3& P );
	void operator-=( const POINT3& r );
	void operator+=( const POINT3& r );
	void operator/=( double r );
	void operator*=( double r );
	void print() const;
	double abs() const;
	double dot(const POINT3& pt) const;

	static POINT3 cross( const POINT3& l, 
		                 const POINT3& r );
	static void   cross( const POINT3& l, 
		                 const POINT3& r,
						       POINT3& out);

	static double area( const POINT3& P0,
		                const POINT3& P1,
						const POINT3& P2 )
	{
		double Ux = P1.x-P0.x; double Uy  = P1.y-P0.y; double Uz = P1.z-P0.z;
		double Vx = P2.x-P0.x; double Vy  = P2.y-P0.y; double Vz = P2.z-P0.z;

		double Dx = Uy*Vz-Vy*Uz;
		double Dy = Vx*Uz-Ux*Vz;
		double Dz = Ux*Vy-Vx*Uy;

		return .5*sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
	}

	static double volume( const POINT3& P0,
		                  const POINT3& P1,
						  const POINT3& P2,
						  const POINT3& P3 )
	{
		double a11 = P1.x-P0.x; double a12 = P1.y-P0.y; double a13 = P1.z-P0.z;
		double a21 = P2.x-P0.x; double a22 = P2.y-P0.y; double a23 = P2.z-P0.z;
		double a31 = P3.x-P0.x; double a32 = P3.y-P0.y; double a33 = P3.z-P0.z;

		double vol = a11*(a22*a33-a23*a32)
                   + a21*(a32*a13-a12*a33)
                   + a31*(a12*a23-a22*a13);

		if(vol>0) return vol/6; else return -vol/6;
	}

	// binary operations
	friend double operator*(const POINT3& P1, const POINT3& P2){ return P1.x*P2.x + P1.y*P2.y + P1.z*P2.z; }
	friend POINT3 operator*(double         s, const POINT3& P ){ return POINT3(s*P.x,s*P.y,s*P.z);}
};

#endif
