#ifndef POINT2_
#define POINT2_

#include <math.h>

//-----------------------------------------------
// 
//  POINT2 class
//
//-----------------------------------------------
class POINT2
{
public:
	double x, y;

	double& operator()(int);
	POINT2() ;
	POINT2(double c1, double c2);
	POINT2& operator=(const POINT2& pt);
	POINT2 (const POINT2& pt);
	POINT2 operator+( const POINT2& r ) const { POINT2 out; out = *this; out += r; return out; }
	POINT2 operator-( const POINT2& r ) const { POINT2 out; out = *this; out -= r; return out; }
	POINT2 operator*( double s          ) const { POINT2 out; out = *this; out *= s; return out; }
	POINT2 operator/( double s          ) const { POINT2 out; out = *this; out /= s; return out; }
	void operator-=( const POINT2& r );
	void operator+=( const POINT2& r );
	void operator/=( double r );
	void operator*=( double r );
	double abs() const;
	double dot(const POINT2& pt) const;
	double sqr()const{ return x*x+y*y; }

	static POINT2 curl( double w, const POINT2& P)
	{
		POINT2 out;
		out.x = -w*P.y;
		out.y =  w*P.x; return out;
	}

	static POINT2 curl( const POINT2& P, double w)
	{
		POINT2 out;
		out.x = w*P.y;
		out.y =-w*P.x; return out;
	}

	static double curl( const POINT2& P1, const POINT2& P2){return P1.x*P2.y-P1.y*P2.x;}

	static double area( const POINT2& P1,
		                const POINT2& P2,
						const POINT2& P3 )
	{
		double sum = .5*((P2.x-P1.x)*(P3.y-P1.y) - (P3.x-P1.x)*(P2.y-P1.y));
		if(sum>0) return sum; else return -sum;
	}

	static double length(const POINT2& P1,
		                 const POINT2& P2)
	{
		return sqrt((P2.x-P1.x)*(P2.x-P1.x)+
			        (P2.y-P1.y)*(P2.y-P1.y));
	}

	static double interpolation(const POINT2& P1, double u1,
								const POINT2& P2, double u2,
								const POINT2& P3, double u3, 
								const POINT2& P )
	{
		double area_123 = area(P1,P2,P3);
		double area_023 = area(P ,P2,P3);
		double area_103 = area(P1,P ,P3);
		double area_120 = area(P1,P2,P );

		return(u1*area_023 +
			   u2*area_103 +
			   u3*area_120 )/area_123;
	}

	static POINT2 cross(double l, const POINT2&        r){return POINT2(-l*r.x, l*r.y); }
	static POINT2 cross(const POINT2& l, double        r){return POINT2( r*l.y,-r*l.x); }
	static double cross(const POINT2& l, const POINT2& r){return l.x*r.y-l.y*r.x; }
};
#endif
