#ifndef MY_DEFINITIONS
#define MY_DEFINITIONS

#include <time.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>

static int BC_NEUMANN   = 1;
static int BC_CONSTANT  = 2;
static int BC_LINEAR    = 3;
static int BC_QUADRATIC = 4;
static int BC_PERIODIC  = 5;

template<class T>
T ABS(T a){ return (a)>0 ? (a) : -(a); }

#ifndef SGN
#define SGN(a) ((a)>0 ? (1) : ((a)<0 ? -1 : 0))
#endif
static int    MIN( int    a, int    b){ return (a<b)?a:b;}
static double MIN( double a, double b){ return (a<b)?a:b;}
static int    MAX( int    a, int    b){ return (a>b)?a:b;}
static double MAX( double a, double b){ return (a>b)?a:b;}
static double MIN( double a, double b, double c )
{
	if(a<b) return (a<c)?a:c;
	else	return (b<c)?b:c;
}
static int MIN( int a, int b, int c )
{
	if(a<b) return (a<c)?a:c;
	else	return (b<c)?b:c;
}
static double MAX( double a, double b, double c )
{
	if(a>b) return (a>c)?a:c;
	else	return (b>c)?b:c;
}
static int MAX( int a, int b, int c )
{
	if(a>b) return (a>c)?a:c;
	else	return (b>c)?b:c;
}
static double MAX( double a, double b, double c, double d ){
	return MAX(MAX(a,b),MAX(c,d));}

static double CUBE(double a){ return a*a*a; }
static double SQR (double a){ return a*a  ; }
static int    SQR (int    a){ return a*a  ; }
static int    CUBE(int    a){ return a*a*a; }

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif

#define NATURAL_NUMBER  2.71828182845904523536

template<class T>
T POS(T a){ return ((a)>0 ? (a) : 0); }

template<class T>
T NEG(T a){ return ((a)<0 ? (a) : 0);}

static double MINMOD(double a, double b){ return ((a)*(b)<0 ? 0 : (ABS(a)<ABS(b) ? (a) : (b)));}

const double eps_machine = 1E-12;

class Function_1D{public:virtual double operator()(double x) const =0;};
class Function_2D{public:virtual double operator()(double x, double y) const =0;};
class Function_3D{public:virtual double operator()(double x, double y, double z) const =0;};
class Function_4D{public:virtual double operator()(double x, double y, double z, double w) const =0;};

class Function_1D_Formula : public Function_1D{double (*formula)(double                     );public:Function_1D_Formula(double (*formula_)(double                     )){formula=formula_;}double operator()(double x                            ) const{return (*formula)(x      );}};
class Function_2D_Formula : public Function_2D{double (*formula)(double,double              );public:Function_2D_Formula(double (*formula_)(double,double              )){formula=formula_;}double operator()(double x,double y                   ) const{return (*formula)(x,y    );}};
class Function_3D_Formula : public Function_3D{double (*formula)(double,double,double       );public:Function_3D_Formula(double (*formula_)(double,double,double       )){formula=formula_;}double operator()(double x,double y, double z         ) const{return (*formula)(x,y,z  );}};
class Function_4D_Formula : public Function_4D{double (*formula)(double,double,double,double);public:Function_4D_Formula(double (*formula_)(double,double,double,double)){formula=formula_;}double operator()(double x,double y, double z,double w) const{return (*formula)(x,y,z,w);}};

class SubFunction2D : public Function_2D
{ 
	const Function_3D* pF; int sub; double subval; 
	public:
	SubFunction2D(const Function_3D& F, int sub_, double subval_){pF=&F;sub=sub_;subval=subval_;} 
	double operator()(double x, double y)const
	{
		     if(sub==1)return (*pF)(subval,x,y);
		else if(sub==2)return (*pF)(x,subval,y);
		else           return (*pF)(x,y,subval);
	}
};

class SubFunction3D : public Function_3D
{ 
	const Function_4D* pF; int sub; double subval; 
	public:
	SubFunction3D(const Function_4D& F, int sub_, double subval_){pF=&F;sub=sub_;subval=subval_;} 
	double operator()(double x, double y, double z)const
	{
		     if(sub==1)return (*pF)(subval,x,y,z);
		else if(sub==2)return (*pF)(x,subval,y,z);
		else if(sub==3)return (*pF)(x,y,subval,z);
		else           return (*pF)(x,y,z,subval);
	}
};

//----------------------------------------------------------------------------------
//
//----------------------------------------------------------------------------------
struct Statistics {
	int number_of_data;
	double sum_of_data;
	double sum_of_data_square;
	double min_of_data;
	double max_of_data;
	Statistics() { 
		number_of_data=0; 
		sum_of_data=0; 
		sum_of_data_square=0;
		min_of_data = 1E8;
		max_of_data = -1E8; }
	void put_data(double data){
		number_of_data++;
		sum_of_data += data;
		sum_of_data_square += data*data;
		if(min_of_data>data) min_of_data=data;
		if(max_of_data<data) max_of_data=data; }
	double average_of_data(){ return sum_of_data/number_of_data; }
	double stdev_of_data(){ return sqrt(sum_of_data_square/number_of_data - average_of_data()*average_of_data());}
	void print_data( const char* name ) {
		if(number_of_data==0) {printf("No data\n"); return; }
		printf("Statistics of %s\n",name);
		printf("maximum is %3.2e\n",max_of_data);
		//printf("minimum is %3.2e\n",min_of_data);
		printf("L1 norm is %3.2e\n",sum_of_data/number_of_data);
		//printf("L2 norm is %3.2e\n",sqrt(sum_of_data_square/number_of_data));
	}};
//-----------------------------------------------------------------------
static double determinant( double a11, double a12, double a13,
	                       double a21, double a22, double a23,
						   double a31, double a32, double a33){
	return a11*(a22*a33-a23*a32)
		  -a12*(a21*a33-a23*a31)
		  +a13*(a21*a32-a22*a31); }
//---------------------------------------------------------------------------
static void solve_3_by_3_linear_system
        (double a11, double a12, double a13, double b1, double& x1,
         double a21, double a22, double a23, double b2, double& x2,
         double a31, double a32, double a33, double b3, double& x3) {
    // Kramer's rule
    double det = a11*(a22*a33-a23*a32)
               - a21*(a12*a33-a13*a32)
               + a31*(a12*a23-a13*a22);
    x1 = (  b1*(a22*a33-a23*a32)
          - b2*(a12*a33-a13*a32)
          + b3*(a12*a23-a13*a22))/det;
    x2 = (  a11*(b2*a33-a23*b3)
          - a21*(b1*a33-a13*b3)
          + a31*(b1*a23-a13*b2))/det;
    x3 = (  a11*(a22*b3-b2*a32)
          - a21*(a12*b3-b1*a32)
          + a31*(a12*b2-b1*a22))/det; }

//---------------------------------------------------------------------
// Stopwatch
//---------------------------------------------------------------------
class StopWatch {
private: 
	clock_t startTime; // time the stop watch was started
	clock_t  stopTime; // stop the stop watch
	int total_ticks;
public:
	StopWatch() {total_ticks=0;}
	void stop (){ stopTime =clock(); } // stop the stopwatch
	void start(){ startTime=clock(); } // start the stopwatch
	void read_duration( const char* name_Of_work) 	{ 
		double duration = (double)(stopTime - startTime) / CLOCKS_PER_SEC;
		printf( "%2.3f seconds in %s\n", duration, name_Of_work );	}
	void read_ticks( const char* name_Of_work) 	{ 
		int number_of_ticks = (int)(stopTime - startTime);
		printf( "%d ticks in %s\n", number_of_ticks, name_Of_work );	}
	void add_into_total_ticks()	{total_ticks += (int)(stopTime - startTime);	}
	void read_total_ticks( const char* name_Of_work) { printf( "%d ticks in %s\n", total_ticks, name_Of_work );}};

//----------------------------------------------------------------------------------------------------
// SQRT
//----------------------------------------------------------------------------------------------------
static double SQRT(double a) {
	double xn   = .5*a+.5;
	double xnp1 = xn/2 + a/(2*xn);
	while ( xn - xnp1 > 1E-14) {
		xn   = xnp1;
		xnp1 = xn/2 + a/(2*xn); }
	return xnp1;}

//----------------------------------------------------------------------------------------------------
// EXPONENTIAL
//----------------------------------------------------------------------------------------------------
static double EXPONENTIAL_01( double x) {
	return  1+x 
		  *(1+x/2 
		  *(1+x/3 
		  *(1+x/4 
		  *(1+x/5 
		  *(1+x/6 
		  *(1+x/7 
		  *(1+x/8 
		  *(1+x/9 
		  *(1+x/10
		  *(1+x/11
		  *(1+x/12
		  *(1+x/13
		  *(1+x/14
		  *(1+x/15)))))))))))))); }

static double EXPONENTIAL_POSITIVE( double x)
{
	double product=1;
	while(x>1)
	{
		product*=NATURAL_NUMBER;
		x      -=1;
	}
	return product*EXPONENTIAL_01(x);
}

static double EXPONENTIAL(double x)
{
	if(x>=0)	return   EXPONENTIAL_POSITIVE( x);
	else		return 1/EXPONENTIAL_POSITIVE(-x);
}

//----------------------------------------------------------------------------------------------------
// SINE & COSECANT
//----------------------------------------------------------------------------------------------------
static double SINE(double x)
{   
	//	move the x in (-PI,PI)
	while(x>PI) x=x-2*PI;
	while(x<PI) x=x+2*PI;

	// in (-PI/2,PI/2)
	if(x> PI/2)	x =  PI-x;
	if(x<-PI/2)	x = -PI-x;
	
	// Taylor expansion with error bound 1E(-15)
	double x2 = x*x;
	return x*(1-x2/  6*
	         (1-x2/ 20*
 	         (1-x2/ 42*
 	         (1-x2/ 72*
 	         (1-x2/110*
 	         (1-x2/156*
 	         (1-x2/210*
 	         (1-x2/272*
   			 (1-x2/342*
			 (1-x2/420*
			 (1-x2/506)))))))))));
}

static double COSECANT(double x){	return 1/SINE(x);	}

//----------------------------------------------------------------------------------------------------
// COSINE & SECANT
//----------------------------------------------------------------------------------------------------
static double COSINE(double x){	return SINE(x+PI/2);	}
static double SECANT(double x){	return 1/SINE(x+PI/2);	}

//----------------------------------------------------------------------------------------------------
// TANGENT & COTANGENT
//----------------------------------------------------------------------------------------------------
static double TANGENT  (double x){	return SINE(x)/COSINE(x);	}
static double COTANGENT(double x){	return COSINE(x)/SINE(x);	}

//----------------------------------------------------------------------------------------------------
// ARCTANGENT & ARCCOTANGENT
//----------------------------------------------------------------------------------------------------
static double ARCTANGENT_ROOT5M2(double x)
{	
	double x2=x*x;
	
	return x*(1-x2*(0.3333333333333333333 -x2
		          *(0.2000000000000000000 -x2
				  *(0.1428571428571428571 -x2
				  *(0.1111111111111111111 -x2
				  *(0.09090909090909090909-x2
				  *(0.07692307692307692308-x2
				  *(0.06666666666666666667-x2
				  *(0.05882352941176470588-x2
				  *(0.05263157894736842105-x2
				  * 0.04761904761904761905 ))))))))));
}

static double ARCTANGENT_ROOT2M1(double x)
{	
	const double ROOT_5 = 2.236067977499789696;
	if ( x < ROOT_5-2)	return      ARCTANGENT_ROOT5M2(x);
	else			    return PI/6-ARCTANGENT_ROOT5M2((1-2*x)/(2+x));
}

static double ARCTANGENT_01(double x)
{
	const double ROOT_2 = 1.414213562373095049;
	if ( x < ROOT_2-1)	return        ARCTANGENT_ROOT2M1(x);
	else			    return .25*PI-ARCTANGENT_ROOT2M1((1-x)/(1+x));
}

static double ARCTANGENT_POSITIVE(double x)
{
	if ( x <= 1)	return       ARCTANGENT_01(  x);
	else			return .5*PI-ARCTANGENT_01(1/x);
}

static double ARCTANGENT(double x)
{
	if ( x >=0 )	return  ARCTANGENT_POSITIVE( x);
	else			return -ARCTANGENT_POSITIVE(-x);
}

static double ARCCOTANGENT(double x){ return  PI/2-ARCTANGENT(x);	}

//----------------------------------------------------------------------------------------------------
// ARCSINE & ARCCOSECANT
//----------------------------------------------------------------------------------------------------
static double ARCSINE    (double x){return 2*ARCTANGENT(x/(1+SQRT(1-x*x)));	}
static double ARCCOSECANT(double x){return 2*ARCTANGENT(1/(x+SQRT(x*x-1)));	}

//----------------------------------------------------------------------------------------------------
// ARCCOSINE & ARCSECANT
//----------------------------------------------------------------------------------------------------
static double ARCCOSINE(double x){return 2*ARCTANGENT(SQRT(1-x*x)/(1+x));	}
static double ARCSECANT(double x){return 2*ARCTANGENT(SQRT(x*x-1)/(1+x));	}

//----------------------------------------------------------------------------------------------------
// subcell resolution
// input : phi_i, phi_ip1, phixx
// output : xGamma \in (0,1)
//----------------------------------------------------------------------------------------------------
static double subcell_resolution( double phi_i, double phi_ip1, double phixx ){
	if(ABS(phixx) < 1E-10)
		return phi_i/(phi_i-phi_ip1);
	else {
		double D = SQR(.5*phixx-phi_i-phi_ip1)-4*phi_i*phi_ip1;
		return .5+(phi_i-phi_ip1+SGN(phi_ip1-phi_i)*SQRT(D))/phixx;}}

//----------------------------------------------------------------------------------------------------
// load distribution for parallel computing
// input : NT(number of threads), i_nodes[]
// output : ibeg_T[], iend_T[]
//----------------------------------------------------------------------------------------------------
static void load_distribution( const int NT, const int NI, const int i_nodes[],
					    int ibeg_T[], int iend_T[] )
{
	// count the total number of nodes
	int NN=0; for(int i=0;i<NI;i++) NN+=i_nodes[i];
	//
	int avg = int(NN/NT)+1;
	//
	int i=0;
	for(int T=0;T<NT;T++)
	{
		ibeg_T[T]=i; int nodes_T=0;
		while(nodes_T<avg && i<NI)
		{   nodes_T+=i_nodes[i]; 
			i++; 
			if(nodes_T==0) ibeg_T[T]=i;
		}
		iend_T[T]=i-1;
	}
}

//---------------------------------------------------------------------
// ENO-2 interpolation
//---------------------------------------------------------------------
static double ENO2_interpolation( double f0, double f1, double f2, double f3, double i )
{
	double f01=(f1-f0);
	double f12=(f2-f1),f012=(f12-f01)*.5;
	double f23=(f3-f2),f123=(f23-f12)*.5;
	     if(f012*f123<0        ){i-=1; return f1*(1-i)+f2*(i);      }
	else if(ABS(f012)<ABS(f123)){      return f0+f01*i+f012*i*(i-1);}
	else						{i-=1; return f1+f12*i+f123*i*(i-1);}
}

static double ENO2_interpolation( double f[4][4], double i, double j )
{
	double fi[4];
	double fj[4];
	for(int a=0;a<4;a++) fi[a]=ENO2_interpolation(f[0][a],f[1][a],f[2][a],f[3][a],i);
	for(int a=0;a<4;a++) fj[a]=ENO2_interpolation(f[a][0],f[a][1],f[a][2],f[a][3],j);
	double fij1 = ENO2_interpolation(fi[0],fi[1],fi[2],fi[3],j);
	double fij2 = ENO2_interpolation(fj[0],fj[1],fj[2],fj[3],i);
	return .5*(fij1+fij2);
}

static double ENO2_interpolation( double f[4][4][4], double i, double j, double k )
{
	double fjk[4];
	double fki[4];
	double fij[4]; double temp[4][4];
	for(int a=0;a<4;a++) fjk[a]=ENO2_interpolation(f[a],j,k);
	for(int b=0;b<4;b++)
	{
		for(int a=0;a<4;a++)
		for(int c=0;c<4;c++)
			temp[a][c]=f[a][b][c];
		fki[b]=ENO2_interpolation(temp,i,k);
	}
	for(int c=0;c<4;c++)
	{
		for(int a=0;a<4;a++)
		for(int b=0;b<4;b++)
			temp[a][b]=f[a][b][c];
		fij[c]=ENO2_interpolation(temp,i,j);
	}
	double fijk1 = ENO2_interpolation(fjk[0],fjk[1],fjk[2],fjk[3],i);
	double fijk2 = ENO2_interpolation(fki[0],fki[1],fki[2],fki[3],j);
	double fijk3 = ENO2_interpolation(fij[0],fij[1],fij[2],fij[3],k);
	return (fijk1+fijk2+fijk3)/3.;
}


#endif
