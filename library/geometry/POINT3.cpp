#include <math.h>
#include <stdio.h>
#include "POINT3.h"

//----------------------------------------------
//
//	POINT3 impl
//
//----------------------------------------------
POINT3::POINT3 () { x=y=z=0; }
POINT3::POINT3 (double c1, double c2, double c3){x = c1;y = c2;z = c3;}
POINT3::POINT3 (const POINT2& pt)
{
	x = pt.x;
	y = pt.y;
	z = 0;
}

double& POINT3::operator()( int i )
{
	if( i==0 ) return x;
	if( i==1 ) return y;
	else       return z;
}

void POINT3::operator=( const POINT3& P )
{
	x = P.x;
	y = P.y;
	z = P.z;
}
void POINT3::operator-=( const POINT3& r )
{
	x -= r.x;
	y -= r.y;
	z -= r.z;
}

void POINT3::operator+=( const POINT3& r )
{
	x += r.x;
	y += r.y;
	z += r.z;
}

void POINT3::operator/=( double r )
{
	x /= r;
	y /= r;
	z /= r;
}

void POINT3::operator*=( double r )
{
	x *= r;
	y *= r;
	z *= r;
}

double POINT3::abs() const
{
	return sqrt( x*x + y*y + z*z );
}

double POINT3::dot( const POINT3& pt ) const
{
	return x*pt.x + y*pt.y + z*pt.z;
}

POINT3 POINT3::cross( const POINT3& l, const POINT3& r )
{
	POINT3 out;

	out.x = l.y * r.z - l.z * r.y;
	out.y =-l.x * r.z + l.z * r.x;
	out.z = l.x * r.y - l.y * r.x;

	return out;
}

void POINT3::cross(  const POINT3& l, 
				     const POINT3& r,
					       POINT3& out )
{
	out.x = l.y * r.z - l.z * r.y;
	out.y =-l.x * r.z + l.z * r.x;
	out.z = l.x * r.y - l.y * r.x;
}

void POINT3::print()const{ printf("% 3.2f, % 3.2f, % 3.2f\n",x,y,z); }
