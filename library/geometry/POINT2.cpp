#include <math.h>
#include "POINT2.h"

//----------------------------------------------
//
//	POINT2 impl
//
//----------------------------------------------
POINT2::POINT2 () { x=y=0; }
POINT2::POINT2 (double c1, double c2)
{
	x = c1;		y = c2;		
}
POINT2& POINT2::operator=(const POINT2& pt)
{
	x = pt.x;	y = pt.y;
	return *this;
}
POINT2::POINT2 (const POINT2& pt)
{
	*this = pt;
}

double& POINT2::operator()( int i )
{
	if( i==0 ) return x;
	else       return y;
}

void POINT2::operator-=( const POINT2& r )
{
	x -= r.x;
	y -= r.y;
}

void POINT2::operator+=( const POINT2& r )
{
	x += r.x;
	y += r.y;
}

void POINT2::operator/=( double r )
{
	x /= r;
	y /= r;
}

void POINT2::operator*=( double r )
{
	x *= r;
	y *= r;
}

double POINT2::abs() const
{
	return sqrt( x*x + y*y );
}

double POINT2::dot( const POINT2& pt ) const
{
	return x*pt.x + y*pt.y;
}

