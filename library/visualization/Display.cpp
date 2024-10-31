#include <math.h>
#include <stdio.h>
#include <iostream>
#include "Display.h"
#include "../STACK.h"
#include "../geometry/MarchingCubes.h"
#include "../uniform_grid/GRID2D.h"
#include "../uniform_grid/SimplicialIsosurfacing.h"
#include "../uniform_grid/TrilinearInterpolation.h"


int RGB_COMBINE( unsigned char r, 
				 unsigned char g, 
				 unsigned char b){ return (r<<16)+(g<<8)+b;}

void RGB_SPLIT(int rgb,  unsigned char& r,
				         unsigned char& g,
				         unsigned char& b)
{
	b = rgb & 0x0000ff; rgb>>=8;
	g = rgb & 0x0000ff; rgb>>=8;
	r = rgb           ;
}

//
typedef unsigned short WORD ; // 16 bit unsigned integer
typedef int            LONG ; // 32 bit integer
typedef unsigned int   DWORD; // 32 bit unsigned integer
   
struct BITMAPINFOHEADER_LINUX{
  DWORD  biSize; 
  LONG   biWidth; 
  LONG   biHeight; 
  WORD   biPlanes; 
  WORD   biBitCount; 
  DWORD  biCompression; 
  DWORD  biSizeImage; 
  LONG   biXPelsPerMeter; 
  LONG   biYPelsPerMeter; 
  DWORD  biClrUsed; 
  DWORD  biClrImportant; 
}; 

struct BITMAPFILEHEADER_LINUX { 
  WORD    bfType; 
  DWORD   bfSize; 
  WORD    bfReserved1; 
  WORD    bfReserved2; 
  DWORD   bfOffBits; 
}; 

//
void Display::save_bmp(const char* file_name) 
{
	// create file name
	char file_name_new[100];

	if(file_name==0)
	{
		static int count = 0;
		sprintf( file_name_new, "image_%04d.bmp", count++ );

		POINT3 pm; from_logi_to_phys_coords(m_imax-150,m_jmax-50,0,pm);
		POINT3 pM; from_logi_to_phys_coords(m_imax- 10,m_jmax-10,0,pM);
		
		// draw_number(count,pm.x,pM.x,pm.y,pM.y,0);
	}
	else
		sprintf( file_name_new, "%s", file_name );

	printf("saving %s\n",file_name_new);

    
	//
    BITMAPFILEHEADER_LINUX bmFileHeader;
    BITMAPINFOHEADER_LINUX bmInfoHeader;

    bmFileHeader.bfSize  = 14;
    bmFileHeader.bfType  = 0x4D42;
    bmFileHeader.bfOffBits = 14 + 40;

    bmInfoHeader.biSize    = 40;
    bmInfoHeader.biWidth    = m_imax+1;
    bmInfoHeader.biHeight   = m_jmax+1;
    bmInfoHeader.biPlanes   = 1;
    bmInfoHeader.biBitCount   = 24; 
    bmInfoHeader.biCompression = 0;
    bmInfoHeader.biSizeImage = bmInfoHeader.biWidth 
                             * bmInfoHeader.biHeight 
                             * (bmInfoHeader.biBitCount/8);
    bmInfoHeader.biXPelsPerMeter = 0;
    bmInfoHeader.biYPelsPerMeter = 0;
    bmInfoHeader.biClrUsed = 0;
    bmInfoHeader.biClrImportant = 0;

	//
	FILE* fp = fopen(file_name_new,"wb");

   long size_w = fwrite((const void*)&(bmFileHeader.bfType     ),1,2,fp);
		size_w = fwrite((const void*)&(bmFileHeader.bfSize     ),1,4,fp);
		size_w = fwrite((const void*)&(bmFileHeader.bfReserved1),1,2,fp);
		size_w = fwrite((const void*)&(bmFileHeader.bfReserved2),1,2,fp);
		size_w = fwrite((const void*)&(bmFileHeader.bfOffBits  ),1,4,fp);

        size_w = fwrite((const void*)&(bmInfoHeader.biSize    	   ),1,4,fp);
		size_w = fwrite((const void*)&(bmInfoHeader.biWidth        ),1,4,fp);
		size_w = fwrite((const void*)&(bmInfoHeader.biHeight       ),1,4,fp);
		size_w = fwrite((const void*)&(bmInfoHeader.biPlanes       ),1,2,fp);
		size_w = fwrite((const void*)&(bmInfoHeader.biBitCount     ),1,2,fp);
		size_w = fwrite((const void*)&(bmInfoHeader.biCompression  ),1,4,fp);
		size_w = fwrite((const void*)&(bmInfoHeader.biSizeImage    ),1,4,fp);
		size_w = fwrite((const void*)&(bmInfoHeader.biXPelsPerMeter),1,4,fp);
		size_w = fwrite((const void*)&(bmInfoHeader.biYPelsPerMeter),1,4,fp);
		size_w = fwrite((const void*)&(bmInfoHeader.biClrUsed      ),1,4,fp);
		size_w = fwrite((const void*)&(bmInfoHeader.biClrImportant ),1,4,fp);

	unsigned char r=0,g=0,b=0; 

    for(int j=0;j<=m_jmax;j++)
    for(int i=0;i<=m_imax;i++)
    {
		RGB_SPLIT(m_rgb[i][j],r,g,b);

		size_w = fwrite((const void*)&b,1,1,fp);
		size_w = fwrite((const void*)&g,1,1,fp);
		size_w = fwrite((const void*)&r,1,1,fp);
    }
                                                                      
    fclose(fp);
}

/*----------------------------------------------------------------*\
 * mapping of logical and physical coordinates
\*----------------------------------------------------------------*/
//
void Display::set_physical_coordinates(double xmin, double xmax, 
									   double ymin, double ymax,
									   double zmin, double zmax)
{
	double xcen = (xmax+xmin)/2, xlen = (xmax-xmin)/2;
	double ycen = (ymax+ymin)/2, ylen = (ymax-ymin)/2;
	double zcen = (zmax+zmin)/2, zlen = (zmax-zmin)/2;

	m_translation.x = -xcen;
	m_translation.y = -ycen;
	m_translation.z = -zcen;
	
	m_scaling.x = 1./(xlen*1.733);
	m_scaling.y = 1./(ylen*1.733);
	m_scaling.z = 1./(zlen*1.733);
}

//
void Display::from_phys_to_logi_coords(const POINT3& P_,double& i, 
													    double& j, 
													    double& z) const
{
	// translation
	POINT3 P(P_);
	P.x += m_translation.x;
	P.y += m_translation.y;
	P.z += m_translation.z;
	
	// scaling
	P.x *= m_scaling.x;
	P.y *= m_scaling.y;
	P.z *= m_scaling.z;
	
	// rotation 
	m_rotation.apply_rotation(P);

	// projection onto xy plane
	z = P.z;

	int max_ = MAX(m_imax,m_jmax);
	i = (int)(.5*(P.x)*max_ + .5*m_imax + .5);
	j = (int)(.5*(P.y)*max_ + .5*m_jmax + .5);
}

//
void Display::from_logi_to_phys_coords( double i, 
										double j, 
										double z, POINT3& P) const
{

	// projection onto xy plane
	P.z=z;
	int max_ = MAX(m_imax,m_jmax);
	P.x = (i-.5-.5*m_imax)*2/max_;
	P.y = (j-.5-.5*m_jmax)*2/max_;

	// rotation
	QUATERNION rotation_inv = m_rotation.inverse();
	rotation_inv.apply_rotation(P);

	// scaling
	P.x /= m_scaling.x;
	P.y /= m_scaling.y;
	P.z /= m_scaling.z;
	
	// translation
	P.x -= m_translation.x;
	P.y -= m_translation.y;
	P.z -= m_translation.z;
}

/*----------------------------------------------------------------*\
 * basic drawings
\*----------------------------------------------------------------*/
void Display::draw_pixel( const POINT3& P, int rgb, int w)
{
	double id,jd,z; from_phys_to_logi_coords(P,id,jd,z);
	
	int i = (int)id;
	int j = (int)jd;

	for(int a=i-w;a<=i+w;a++)
	for(int b=j-w;b<=j+w;b++)
	if( a>=0 && a<=m_imax &&
		b>=0 && b<=m_jmax && z>=m_z[a][b])
	{
		m_rgb[a][b] = rgb;
		m_z  [a][b] = z  ;
	}
}

void Display::draw_line( double x0, double y0,
						 double x1, double y1, int rgb )
{
	draw_line(POINT3(x0,y0,0),POINT3(x1,y1,0),rgb);
}

void Display::POINT( int i, int j, double z, int rgb )
{	
	if(i>=0 && i<=m_imax &&
	   j>=0 && j<=m_jmax && z>=m_z[i][j])
	{
		m_rgb[i][j]=rgb;
		m_z  [i][j]=z;
	}
}

void Display::draw_line( const POINT3& P0,
	                     const POINT3& P1, int rgb )
{
	// the line segment is in rectangle [imin,imax]x[jmin,jmax]
	double i0,j0,z0; from_phys_to_logi_coords(P0,i0,j0,z0);
	double i1,j1,z1; from_phys_to_logi_coords(P1,i1,j1,z1);
	
	// clipping
	if( i0<0 || i0>m_imax || i1<0 || i1>m_imax ||
		j0<0 || j0>m_jmax || j1<0 || j1>m_jmax ) return;

	// singular case
	int di = (int)(i1-i0);
	int dj = (int)(j1-j0);
	if( di==0 && dj==0 ) 
	{
		if(z0<z1) POINT((int)(i0),(int)(j0),z0,rgb);
		else      POINT((int)(i1),(int)(j1),z1,rgb);
		return;
	}

	// x-y symmetry
	if( ABS(di)>=ABS(dj) )
	{
		double sj = ((double)j1-j0)/(i1-i0); double j=j0;
		double sz = ((double)z1-z0)/(i1-i0); double z=z0;

		if( di >= 0 )
			for( int i=(int)(i0); i<=i1; i++, j+=sj, z+=sz)
				POINT( i,int(j+.5),z, rgb);
		else
			for( int i=(int)(i0); i>=i1; i--, j-=sj, z-=sz)
				POINT( i,int(j+.5),z, rgb);
	}
	else
	{

		double si = ((double)i1-i0)/(j1-j0); double i=i0;
		double sz = ((double)z1-z0)/(j1-j0); double z=z0;

		if( dj >= 0 )
			for( int j=(int)(j0); j<=j1; j++, i+=si, z+=sz)
				POINT( int(i+.5),j,z, rgb);
		else
			for( int j=(int)(j0); j>=j1; j--, i-=si, z-=sz)
				POINT( int(i+.5),j,z, rgb);	
	}
}

//
void Display::draw_triangle(const POINT3& P1,
							const POINT3& P2,
							const POINT3& P3, int rgb)
{
	double i1,j1,z1; from_phys_to_logi_coords(P1,i1,j1,z1);
	double i2,j2,z2; from_phys_to_logi_coords(P2,i2,j2,z2);
	double i3,j3,z3; from_phys_to_logi_coords(P3,i3,j3,z3);

    int imin=(int)(MIN(MIN(i1,i2),i3))-1; int imax=(int)(MAX(MAX(i1,i2),i3))+1;
	int jmin=(int)(MIN(MIN(j1,j2),j3))-1; int jmax=(int)(MAX(MAX(j1,j2),j3))+1;
	
	if(imin<0) imin=0; if(imax>m_imax) imax=m_imax;
	if(jmin<0) jmin=0; if(jmax>m_jmax) jmax=m_jmax;

    double det = i2*j1 - i1*j2 + i1*j3 - j1*i3 - i2*j3 + i3*j2;
    
    if(det<1E-13 && det>-1E-13) return; // singular case

    for(int i=imin;i<=imax;i++)
    for(int j=jmin;j<=jmax;j++)
    {
        double l1 = (j*i2 - i*j2 + i*j3 - j*i3 - i2*j3 + i3*j2)/det; if(l1<-1E-5) continue;
        double l2 = (i*j1 - j*i1 - i*j3 + j*i3 + i1*j3 - j1*i3)/det; if(l2<-1E-5) continue;
        double l3 = (j*i1 - i*j1 + i*j2 - j*i2 - i1*j2 + i2*j1)/det; if(l3<-1E-5) continue;

        double z = l1*z1 + l2*z2 + l3*z3;

        if(z>=m_z[i][j])
        {
            m_rgb[i][j] = rgb;
            m_z  [i][j] = z;
        }
    }
}


//
void Display::draw_curve( double x0, double y0,double a0, double b0,
			              double x1, double y1,double a1, double b1, int rgb)
{
    double c3 = 2*(x0-x1)+(a0+a1);
    double c2 = 3*(x1-x0)-(2*a0+a1);
    double c1 = a0;
    double c0 = x0;

    double d3 = 2*(y0-y1)+(b0+b1);
    double d2 = 3*(y1-y0)-(2*b0+b1);
    double d1 = b0;
    double d0 = y0;

	int N=15;

    for(int n=0;n<N;n++)
    {
        double tmin = ((double)n   )/N;
        double tmax = ((double)n+1.)/N;

        double xmin = ((c3*tmin+c2)*tmin+c1)*tmin+c0;
        double ymin = ((d3*tmin+d2)*tmin+d1)*tmin+d0;

        double xmax = ((c3*tmax+c2)*tmax+c1)*tmax+c0;
        double ymax = ((d3*tmax+d2)*tmax+d1)*tmax+d0;

        draw_line(xmin,ymin,xmax,ymax,rgb);
    }
}

//
void Display::draw_number_fr_0_to_9( int number, double xmin, double xmax,
								                 double ymin, double ymax, int rgb)
{
	// scale
    double rect_w = xmax-xmin; double scale_w = rect_w/0.7;
    double rect_h = ymax-ymin; double scale_h = rect_h/1.0;
    double scale = (scale_w>scale_h) ? scale_h : scale_w;

    // translation
    double trans_x = (xmax+xmin)/2-.35*scale;
    double trans_y = (ymax+ymin)/2-.5 *scale;

    // Data of Bezier curves for drawing numbers 0 to 9
    static const double c0_1[8]={0.359,0.847,-1.079, 0.0  ,0.359,0.133, 1.079, 0.0  };
    static const double c0_2[8]={0.359,0.133, 1.062, 0.0  ,0.359,0.847,-1.062, 0.0  };
    static const double c1_1[8]={0.190,0.706, 0.300, 0.048,0.356,0.898, 0.083, 0.400};
    static const double c1_2[8]={0.356,0.898, 0.000,-0.400,0.348,0.138, 0.000,-0.105};
    static const double c2_1[8]={0.067,0.095, 0.678, 0.080,0.631,0.089, 0.767, 0.246};
    static const double c2_2[8]={0.067,0.828, 1.316, 0.466,0.067,0.095,-2.118,-1.316};
    static const double c3_1[8]={0.146,0.741, 0.881, 1.102,0.307,0.517,-1.774,-0.197};
    static const double c3_2[8]={0.307,0.517, 1.130,-0.163,0.179,0.116,-1.915, 0.031};
    static const double c4_1[8]={0.198,0.879,-0.375,-1.208,0.580,0.424, 1.497, 0.220};
    static const double c4_2[8]={0.433,0.894,-0.163,-0.704,0.430,0.157,-0.197,-0.532};
    static const double c5_1[8]={0.561,0.803,-0.065,-0.071,0.166,0.809,-0.246,-0.057};
    static const double c5_2[8]={0.166,0.809,-0.057,-0.294,0.154,0.520,-0.114,-0.180};
    static const double c5_3[8]={0.154,0.520, 1.571,-0.088,0.173,0.111,-1.703,-0.105};
    static const double c6_1[8]={0.356,0.100, 0.752,-0.137,0.236,0.539,-1.585,-0.303};
    static const double c6_2[8]={0.446,0.869,-0.251,-0.687,0.356,0.100, 1.445,-0.466};
    static const double c7_1[8]={0.133,0.831, 0.463, 0.048,0.555,0.836, 0.349, 0.148};
    static const double c7_2[8]={0.555,0.836,-0.349,-0.417,0.359,0.149, 0.237,-0.695};
    static const double c8_1[8]={0.343,0.831,-1.660,-0.440,0.337,0.166,-1.683,-0.589};
    static const double c8_2[8]={0.337,0.166,-1.797, 0.449,0.343,0.831,-1.740, 0.572};
    static const double c9_1[8]={0.173,0.787,-0.535,-0.400,0.531,0.550, 0.924, 0.744};
    static const double c9_2[8]={0.173,0.787, 1.047, 0.646,0.288,0.138,-1.245,-0.638};

	static const double ca_1[8]={0.4303, 0.8397, -1.3626, -0.3865, 0.5229, 0.2863,  1.6317,  0.9160};
	static const double ca_2[8]={0.4303, 0.8397,  0.8101,  0.2691, 0.6594, 0.1851,  0.8588,  0.2948};
	static const double cb_1[8]={0.1469, 0.8721,  0.0086, -0.6126, 0.1632, 0.1164, -0.0315,  0.4122};
	static const double cb_2[8]={0.1718, 0.0954,  1.5372, -0.1660, 0.1574, 0.4962, -1.2767, -0.4408};
	static const double cc_1[8]={0.5067, 0.7987, -1.3569,  0.0000, 0.4790, 0.1441,  1.7548,  0.0315};
	static const double cd_1[8]={0.4303, 0.5200, -1.2968,  0.0172, 0.4332, 0.1469,  1.5029,  0.0344};
	static const double cd_2[8]={0.4408, 0.8989,  0.0258, -0.8645, 0.6460, 0.1059,  0.8674,  0.4380};
	static const double ce_1[8]={0.5611, 0.4933, -0.0172,  2.4132, 0.5639, 0.2882,  3.2118,  0.6126};
	static const double ce_2[8]={0.5506, 0.4981, -0.5410, -0.1374, 0.1002, 0.4905, -0.5010,  0.0401};
	static const double cf_1[8]={0.1632, 0.5258,  0.2462,  0.0172, 0.4198, 0.5172,  0.1718, -0.0515};
	static const double cf_2[8]={0.4437, 0.6489, -0.1231,  0.7185, 0.2939, 0.1498,  0.2691, -2.4447};

	//
    const double* curves_to_draw[3];
    if(number==0){curves_to_draw[0]=c0_1; curves_to_draw[1]=c0_2; curves_to_draw[2]=0;}
    if(number==1){curves_to_draw[0]=c1_1; curves_to_draw[1]=c1_2; curves_to_draw[2]=0;}
    if(number==2){curves_to_draw[0]=c2_1; curves_to_draw[1]=c2_2; curves_to_draw[2]=0;}
    if(number==3){curves_to_draw[0]=c3_1; curves_to_draw[1]=c3_2; curves_to_draw[2]=0;}
    if(number==4){curves_to_draw[0]=c4_1; curves_to_draw[1]=c4_2; curves_to_draw[2]=0;}
    if(number==5){curves_to_draw[0]=c5_1; curves_to_draw[1]=c5_2; curves_to_draw[2]=c5_3;}
    if(number==6){curves_to_draw[0]=c6_1; curves_to_draw[1]=c6_2; curves_to_draw[2]=0;}
    if(number==7){curves_to_draw[0]=c7_1; curves_to_draw[1]=c7_2; curves_to_draw[2]=0;}
    if(number==8){curves_to_draw[0]=c8_1; curves_to_draw[1]=c8_2; curves_to_draw[2]=0;}
    if(number==9){curves_to_draw[0]=c9_1; curves_to_draw[1]=c9_2; curves_to_draw[2]=0;}

	if(number==10){curves_to_draw[0]=ca_1; curves_to_draw[1]=ca_2; curves_to_draw[2]=0;}
	if(number==11){curves_to_draw[0]=cb_1; curves_to_draw[1]=cb_2; curves_to_draw[2]=0;}
	if(number==12){curves_to_draw[0]=cc_1; curves_to_draw[1]=   0; curves_to_draw[2]=0;}
	if(number==13){curves_to_draw[0]=cd_1; curves_to_draw[1]=cd_2; curves_to_draw[2]=0;}
	if(number==14){curves_to_draw[0]=ce_1; curves_to_draw[1]=ce_2; curves_to_draw[2]=0;}
	if(number==15){curves_to_draw[0]=cf_1; curves_to_draw[1]=cf_2; curves_to_draw[2]=0;}
    // 
    for(int n=0;n<3;n++)
    {
        const double* curve_to_draw = curves_to_draw[n];

        if(curve_to_draw!=0)
        {
            draw_curve(curve_to_draw[0]*scale + trans_x,
                       curve_to_draw[1]*scale + trans_y,
                       curve_to_draw[2]*scale          ,
                       curve_to_draw[3]*scale          ,
                       curve_to_draw[4]*scale + trans_x,
                       curve_to_draw[5]*scale + trans_y,
                       curve_to_draw[6]*scale          ,
                       curve_to_draw[7]*scale          ,rgb);
        }
    }
}

//
void Display::draw_signum( int signum, double xmin, double xmax,
						               double ymin, double ymax, int rgb)
{
	// scale
    double rect_w = xmax-xmin; double scale_w = rect_w/0.7;
    double rect_h = ymax-ymin; double scale_h = rect_h/1.0;
    double scale = (scale_w>scale_h) ? scale_h : scale_w;

    // translation
    double trans_x = (xmax+xmin)/2-.35*scale;
    double trans_y = (ymax+ymin)/2-.5 *scale;

    // Data of Bezier curves 
    static const double minus[8]={0.1,0.5, 1.0,0.0, 0.6 ,0.5,-1.0, 0.0};
    static const double plus [8]={0.5,0.1, 0.0,1.0, 0.35,0.9, 0.0,-1.0};

	draw_curve(	minus[0]*scale + trans_x,
				minus[1]*scale + trans_y,
				minus[2]*scale          ,
				minus[3]*scale          ,
				minus[4]*scale + trans_x,
				minus[5]*scale + trans_y,
				minus[6]*scale          ,
				minus[7]*scale          ,rgb);

	if(signum>0)
	{
		draw_curve(	plus[0]*scale + trans_x,
					plus[1]*scale + trans_y,
					plus[2]*scale          ,
					plus[3]*scale          ,
					plus[4]*scale + trans_x,
					plus[5]*scale + trans_y,
					plus[6]*scale          ,
					plus[7]*scale          ,rgb);
	}
}

//
void Display::draw_number( int integer, double xmin, double xmax,
					                    double ymin, double ymax, int rgb)
{
	// calculation of scale and center
	bool b_negative; if(integer<0){ b_negative=true; integer=-integer; }
    int num_digits=0;
	
	int integer_bak = integer;
	while(integer_bak>0)
	{
		num_digits++;
		integer_bak/=10;
	}

	if(integer==0) num_digits=1;
	if(b_negative) num_digits++;

    double center_x = (xmin+xmax)/2; double width  = xmax-xmin;
    double center_y = (ymin+ymax)/2; double height = ymax-ymin;

    double scale_x = width /(num_digits*0.7);
    double scale_y = height/(1.0);

    if(scale_x>scale_y){ scale_x=scale_y; width =scale_x*num_digits*0.7;}
    else               { scale_y=scale_x; height=scale_y*           1.0;}

    // drawing each digit
    for(int d=0;d<num_digits;d++)
    {
        double ymin_d = center_y - height/2;
        double ymax_d = center_y + height/2;

        double xmin_d = center_x - width/2 + scale_x*0.7*(num_digits-1-d);
        double xmax_d = center_x - width/2 + scale_x*0.7*(num_digits  -d);

		if(b_negative && d==num_digits-1)
			draw_signum(-1, xmin_d,xmax_d,ymin_d,ymax_d,rgb);
		else
		{
			int integer_d = integer%10;
        	    integer   = integer/10;

			draw_number_fr_0_to_9(integer_d,xmin_d,xmax_d,
											ymin_d,ymax_d,rgb);
		}
    }
}

//
void Display::draw_number_hexa( int integer, double xmin, double xmax,
					                         double ymin, double ymax, int rgb)
{
	// calculation of scale and center
	bool b_negative; if(integer<0){ b_negative=true; integer=-integer; }
    int num_digits=0;
	
	int integer_bak = integer;
	while(integer_bak>0)
	{
		num_digits++;
		integer_bak/=16;
	}

	if(integer==0) num_digits=1;
	if(b_negative) num_digits++;

    double center_x = (xmin+xmax)/2; double width  = xmax-xmin;
    double center_y = (ymin+ymax)/2; double height = ymax-ymin;

    double scale_x = width /(num_digits*0.7);
    double scale_y = height/(1.0);

    if(scale_x>scale_y){ scale_x=scale_y; width =scale_x*num_digits*0.7;}
    else               { scale_y=scale_x; height=scale_y*           1.0;}

    // drawing each digit
    for(int d=0;d<num_digits;d++)
    {
        double ymin_d = center_y - height/2;
        double ymax_d = center_y + height/2;

        double xmin_d = center_x - width/2 + scale_x*0.7*(num_digits-1-d);
        double xmax_d = center_x - width/2 + scale_x*0.7*(num_digits  -d);

    	if(b_negative && d==num_digits-1)
			draw_signum(-1, xmin_d,xmax_d,ymin_d,ymax_d,rgb);
		else
		{
			int integer_d = integer%16;
        	    integer   = integer/16;

			draw_number_fr_0_to_9(integer_d,xmin_d,xmax_d,
											ymin_d,ymax_d,rgb);
		}
    }
}

//
void Display::draw_number( double integer, double xmin, double xmax,
					                       double ymin, double ymax, int rgb)
{
	int i1 = (int)integer;
	int i2 = (int)((integer - i1)*100);
	
	double xcen =   .5*(xmin+xmax);
	double ycen =   .5*(ymin+ymax);
	double  len = .009*(xmax-xmin);
	
	draw_number(i1,xmin,xcen,ymin,ymax,rgb);
	draw_number(i2,xcen,xmax,ymin,ymax,rgb);
	box(xcen-len,xcen+len,
		ycen-len,ycen+len,
		    -len,     len);
}

//
int Display::Phong_Shading( POINT3 N, int rgb ) const
{
	// default setting
	POINT3 L(-1,1,1.8); L/=L.abs();
	                    N/=N.abs();
	POINT3 V(0,0,1);
	
	m_rotation.apply_rotation(N);

	if(rgb==0) rgb = 0xFF0000;

	double Cd_b = (double)(rgb & 0x0000ff); double Cs_b=Cd_b; double ks_b = 60; rgb>>=8;
	double Cd_g = (double)(rgb & 0x0000ff); double Cs_g=Cd_g; double ks_g = 60; rgb>>=8;
	double Cd_r = (double)(rgb & 0x0000ff); double Cs_r=Cd_r; double ks_r = 60;

	// reflection vector
	double NL = N*L;
	POINT3 R = (2*NL)*N - L;
	
	// Phong shading colour
	double RV = R*V;

    if(NL<0) NL=-NL; // no negative colour
    if(RV<0) RV=-RV; // no negative colour
    
    RV=0;
    
    int r = (int)(Cd_r*NL + Cs_r*pow(RV,ks_r));
    int g = (int)(Cd_g*NL + Cs_g*pow(RV,ks_g));
    int b = (int)(Cd_b*NL + Cs_b*pow(RV,ks_b));

    if(r>255) r=255; // maximum intensity is 255
    if(g>255) g=255;
    if(b>255) b=255;

    return (r<<16)+(g<<8)+(b);
}


//
void Display::box(double xmin, double xmax,
				  double ymin, double ymax,
				  double zmin, double zmax, int rgb)
{
	draw_line(POINT3(xmin,ymin,zmin), POINT3(xmin,ymin,zmax),rgb);
	draw_line(POINT3(xmin,ymin,zmin), POINT3(xmin,ymax,zmin),rgb);
	draw_line(POINT3(xmin,ymin,zmin), POINT3(xmax,ymin,zmin),rgb);
	draw_line(POINT3(xmax,ymin,zmin), POINT3(xmax,ymax,zmin),rgb);
	draw_line(POINT3(xmax,ymin,zmin), POINT3(xmax,ymin,zmax),rgb);
	draw_line(POINT3(xmin,ymax,zmin), POINT3(xmax,ymax,zmin),rgb);
	draw_line(POINT3(xmin,ymax,zmin), POINT3(xmin,ymax,zmax),rgb);
	draw_line(POINT3(xmin,ymin,zmax), POINT3(xmax,ymin,zmax),rgb);
	draw_line(POINT3(xmin,ymin,zmax), POINT3(xmin,ymax,zmax),rgb);
	draw_line(POINT3(xmax,ymin,zmax), POINT3(xmax,ymax,zmax),rgb);
	draw_line(POINT3(xmin,ymax,zmax), POINT3(xmax,ymax,zmax),rgb);
	draw_line(POINT3(xmax,ymax,zmin), POINT3(xmax,ymax,zmax),rgb);
}

/*----------------------------------------------------------------*\
 * Drawings of functions
\*----------------------------------------------------------------*/
void Display::graph( const ARRAY<double>& F, int rgb)
{
	// finding fmin and famx
	double fmin=1E10, fmax=-1E10;
	for(int i=0;i<F.size;i++){ double f=F[i];
		if(f<fmin) fmin=f;
		if(f>fmax) fmax=f;}
	set_physical_coordinates(0,F.size-1,fmin,fmax);
	box                     (0,F.size-1,fmin,fmax,-1,1);
	// draw graph
	for(int i=0;i<F.size-1;i++)
	{
		//
		double l=i-.005/m_scaling.x, b=F[i]-.005/m_scaling.y;
		double r=i+.005/m_scaling.x, t=F[i]+.005/m_scaling.y;
		draw_line(l,b,l,t,rgb);
		draw_line(l,b,r,b,rgb);
		draw_line(r,t,l,t,rgb);
		draw_line(r,t,r,b,rgb);
		//
		draw_line(i,F[i],i+1,F[i+1],rgb);
	}
}


//
void Display::graph( const ARRAY2D<double>& F, int rgb, bool wire_frame) {
	//
	double dx=1.0/F.isize;
	double dy=1.0/F.jsize;
	// finding fmin and famx
	double fmin=1E10, fmax=-1E10;
	for(int i=0;i<F.isize;i++)
	for(int j=0;j<F.jsize;j++){ double f=F[i][j];
		if(f<fmin) fmin=f;
		if(f>fmax) fmax=f;}
	// draw box
	set_physical_coordinates(0,F.isize-1,0,F.jsize-1,fmin,fmax);
	box                     (0,F.isize-1,0,F.jsize-1,fmin,fmax);
	// draw graph
	for(int i=0;i<F.isize-1;i++){
	for(int j=0;j<F.jsize-1;j++){
		double f00=F[i  ][j  ]; POINT3 P00(i  ,j  ,f00);
		double f01=F[i  ][j+1]; POINT3 P01(i  ,j+1,f01);
		double f10=F[i+1][j  ]; POINT3 P10(i+1,j  ,f10);
		double f11=F[i+1][j+1]; POINT3 P11(i+1,j+1,f11);
		if(wire_frame) {
			draw_line(P00,P01,rgb);
			draw_line(P00,P10,rgb);
			draw_line(P11,P01,rgb);
			draw_line(P11,P10,rgb); }
		else{
			double fhh=.25*(f00+f01+f10+f11); POINT3 Phh(i+.5,j+.5,fhh);
			POINT3 Nh0(1,-(f10-f00)/dx,-(2*fhh-f00-f10)/dy); int rgb_h0=Phong_Shading(Nh0,rgb);
			POINT3 Nh1(1,-(f11-f01)/dx,-(f01+f11-2*fhh)/dy); int rgb_h1=Phong_Shading(Nh1,rgb);
			POINT3 N0h(1,-(2*fhh-f00-f01)/dx,-(f01-f00)/dy); int rgb_0h=Phong_Shading(N0h,rgb);
			POINT3 N1h(1,-(f10+f11-2*fhh)/dx,-(f11-f10)/dy); int rgb_1h=Phong_Shading(N1h,rgb);
			draw_triangle(Phh,P00,P10,rgb_h0);
			draw_triangle(Phh,P10,P11,rgb_1h);
			draw_triangle(Phh,P11,P01,rgb_h1);
			draw_triangle(Phh,P01,P00,rgb_0h);}}}
}

//
void Display::graph( const ARRAY3D<double>& F ) {
	//
	double dx=1.0/F.isize;
	double dy=1.0/F.jsize;
	double dz=1.0/F.ksize;
	// draw box
	set_physical_coordinates(0,F.isize-1,0,F.jsize-1,0,F.ksize-1);
	int isize = F.isize;
	int jsize = F.jsize;
	int ksize = F.ksize;
	int num_of_segments=50; double EPS=1E-7;

	GRID3D grid; grid.set_grid(0,F.isize-1,F.isize,
		                       0,F.jsize-1,F.jsize,
							   0,F.ksize-1,F.ksize);
	TrilinearInterpolation Fcont(grid,F);
	for(int a=0;a<=m_imax;a++)
	for(int b=0;b<=m_jmax;b++)
	{
		POINT3 X0; from_logi_to_phys_coords(a,b,-1,X0);
		POINT3 X1; from_logi_to_phys_coords(a,b, 1,X1);
		// line : (1-t)X0 + t*X1
		double i0=X0.x, i1=X1.x; if(i1==i0) i1=i0+EPS;
		double j0=X0.y, j1=X1.y; if(j1==j0) j1=j0+EPS;
		double k0=X0.z, k1=X1.z; if(k1==k0) k1=k0+EPS;

		double tmin = MAX(MIN(-i0/(i1-i0),(isize-1-i0)/(i1-i0)),
						  MIN(-j0/(j1-j0),(jsize-1-j0)/(j1-j0)),
						  MIN(-k0/(k1-k0),(ksize-1-k0)/(k1-k0)));
		double tmax = MIN(MAX(-i0/(i1-i0),(isize-1-i0)/(i1-i0)),
						  MAX(-j0/(j1-j0),(jsize-1-j0)/(j1-j0)),
						  MAX(-k0/(k1-k0),(ksize-1-k0)/(k1-k0)));
		double sum=0;
		if(tmin<tmax)
		{
			double  t = tmin;
			double dt = (tmax-tmin)/num_of_segments;
			double  f = Fcont(	(1-t)*i0+t*i1,
								(1-t)*j0+t*j1,
								(1-t)*k0+t*k1);
			for(int n=1;n<num_of_segments;n++)
			{
				t+=dt;
				double f_ = Fcont(	(1-t)*i0+t*i1,
									(1-t)*j0+t*j1,
									(1-t)*k0+t*k1);
				sum += .5*(f+f_)*dt; f=f_;
			}
		}
		{ 
			if(sum>0)
				const char* test="test";
			char red,gre,blu; red=gre=blu= 1.5*255*sum;	
			m_rgb[a][b] = (red<<16)+(gre<<8)+blu;
			m_z  [a][b] = -3;
		}
	}
	box(0,F.isize-1,0,F.jsize-1,0,F.ksize-1,0xffffff);
}


//
void Display::contour( const ARRAY2D<double>& f, int rgb, double level )
{
	// draw box
	set_physical_coordinates(0,f.isize-1,0,f.jsize-1,-1,1);
	
	m_scaling.x = MIN(m_scaling.x,m_scaling.y);
	m_scaling.y =     m_scaling.x;

	// box(0,f.isize-1,0,f.jsize-1,-1,1);

	// drawing
	GRID2D grid; grid.set_grid(0,f.isize-1,f.isize,
		                       0,f.jsize-1,f.jsize);
	ARRAY<POINT2> Vs,Ns;
	SimplicialIsosurfacing::isosurface( grid, f, Vs,Ns,level);
	
	for(int n=0;n<Vs.size/2;n++)
	{
		POINT2 P0=Vs[2*n  ];
		POINT2 P1=Vs[2*n+1];
		draw_line(P0.x,P0.y,P1.x,P1.y,rgb); 
	}
}

//
void Display::graph_xconst( const ARRAY3D<double>& F,int rgb,bool wire_frame)
{
	//
	double dz=1.0/F.ksize;
	double dy=1.0/F.jsize;
	// finding fmin and famx
	int i=F.isize/2;
	double fmin=1E10, fmax=-1E10;
	for(int k=0;k<F.ksize;k++)
	for(int j=0;j<F.jsize;j++){ double f=F[i][j][k];
		if(f<fmin) fmin=f;
		if(f>fmax) fmax=f;}
	// draw box
	set_physical_coordinates(fmin,fmax,0,F.jsize-1,0,F.ksize-1);
	box                     (fmin,fmax,0,F.jsize-1,0,F.ksize-1);
	for(int j=0;j<F.jsize-1;j++){ 
	for(int k=0;k<F.ksize-1;k++){ 
		double f00=F[i][j  ][k  ]; POINT3 P00(f00,j  ,k  );
		double f01=F[i][j  ][k+1]; POINT3 P01(f01,j  ,k+1);
		double f10=F[i][j+1][k  ]; POINT3 P10(f10,j+1,k  );
		double f11=F[i][j+1][k+1]; POINT3 P11(f11,j+1,k+1);
		if(wire_frame) {
			draw_line(P00,P01,rgb);
			draw_line(P00,P10,rgb);
			draw_line(P11,P01,rgb);
			draw_line(P11,P10,rgb); }
		else{
			double fhh=.25*(f00+f01+f10+f11); POINT3 Phh(fhh,j+.5,k+.5);
			POINT3 Nh0(1,-(f10-f00)/dy      ,-(2*fhh-f00-f10)/dz); int rgb_h0=Phong_Shading(Nh0,rgb);
			POINT3 Nh1(1,-(f11-f01)/dy      ,-(f01+f11-2*fhh)/dz); int rgb_h1=Phong_Shading(Nh1,rgb);
			POINT3 N0h(1,-(2*fhh-f00-f01)/dy,-(f01-f00)/dz      ); int rgb_0h=Phong_Shading(N0h,rgb);
			POINT3 N1h(1,-(f10+f11-2*fhh)/dy,-(f11-f10)/dz      ); int rgb_1h=Phong_Shading(N1h,rgb);
			draw_triangle(Phh,P00,P10,rgb_h0);
			draw_triangle(Phh,P10,P11,rgb_1h);
			draw_triangle(Phh,P11,P01,rgb_h1);
			draw_triangle(Phh,P01,P00,rgb_0h);}}}
}

//
void Display::graph_yconst( const ARRAY3D<double>& F, int rgb, bool wire_frame)
{
	//
	double dz=1.0/F.ksize;
	double dx=1.0/F.isize;
	// finding fmin and famx
	int j=F.jsize/2;
	double fmin=1E10, fmax=-1E10;
	for(int k=0;k<F.ksize;k++)
	for(int i=0;i<F.isize;i++){ double f=F[i][j][k];
		if(f<fmin) fmin=f;
		if(f>fmax) fmax=f;}
	// draw box
	set_physical_coordinates(0,F.isize-1,fmin,fmax,0,F.ksize-1);
	box                     (0,F.isize-1,fmin,fmax,0,F.ksize-1);
	for(int i=0;i<F.isize-1;i++){ 
	for(int k=0;k<F.ksize-1;k++){ 
		double f00=F[i  ][j][k  ]; POINT3 P00(i  ,f00,k  );
		double f01=F[i  ][j][k+1]; POINT3 P01(i  ,f01,k+1);
		double f10=F[i+1][j][k  ]; POINT3 P10(i+1,f10,k  );
		double f11=F[i+1][j][k+1]; POINT3 P11(i+1,f11,k+1);
		if(wire_frame) {
			draw_line(P00,P01,rgb);
			draw_line(P00,P10,rgb);
			draw_line(P11,P01,rgb);
			draw_line(P11,P10,rgb); }
		else{
			double fhh=.25*(f00+f01+f10+f11); POINT3 Phh(i+.5,fhh,k+.5);
			POINT3 Nh0(-(f10-f00)/dx      ,1,-(2*fhh-f00-f10)/dz); int rgb_h0=Phong_Shading(Nh0,rgb);
			POINT3 Nh1(-(f11-f01)/dx      ,1,-(f01+f11-2*fhh)/dz); int rgb_h1=Phong_Shading(Nh1,rgb);
			POINT3 N0h(-(2*fhh-f00-f01)/dx,1,-(f01-f00)/dz      ); int rgb_0h=Phong_Shading(N0h,rgb);
			POINT3 N1h(-(f10+f11-2*fhh)/dx,1,-(f11-f10)/dz      ); int rgb_1h=Phong_Shading(N1h,rgb);
			draw_triangle(Phh,P00,P10,rgb_h0);
			draw_triangle(Phh,P10,P11,rgb_1h);
			draw_triangle(Phh,P11,P01,rgb_h1);
			draw_triangle(Phh,P01,P00,rgb_0h);}}}
}


//
void Display::graph_zconst( const ARRAY3D<double>& F, int rgb, bool wire_frame)
{
	//
	double dy=1.0/F.jsize;
	double dx=1.0/F.isize;
	// finding fmin and famx
	int k=F.ksize/2;
	double fmin=1E10, fmax=-1E10;
	for(int j=0;j<F.jsize;j++)
	for(int i=0;i<F.isize;i++){ double f=F[i][j][k];
		if(f<fmin) fmin=f;
		if(f>fmax) fmax=f;}
	// draw box
	set_physical_coordinates(0,F.isize-1,0,F.jsize-1,fmin,fmax);
	box                     (0,F.isize-1,0,F.jsize-1,fmin,fmax);
	for(int i=0;i<F.isize-1;i++){ 
	for(int j=0;j<F.jsize-1;j++){ 
		double f00=F[i  ][j  ][k]; POINT3 P00(i  ,j  ,f00);
		double f01=F[i  ][j+1][k]; POINT3 P01(i  ,j+1,f01);
		double f10=F[i+1][j  ][k]; POINT3 P10(i+1,j  ,f10);
		double f11=F[i+1][j+1][k]; POINT3 P11(i+1,j+1,f11);
		if(wire_frame) {
			draw_line(P00,P01,rgb);
			draw_line(P00,P10,rgb);
			draw_line(P11,P01,rgb);
			draw_line(P11,P10,rgb); }
		else{
			double fhh=.25*(f00+f01+f10+f11); POINT3 Phh(i+.5,j+.5,fhh);
			POINT3 Nh0(-(f10-f00)/dx      ,-(2*fhh-f00-f10)/dy,1); int rgb_h0=Phong_Shading(Nh0,rgb);
			POINT3 Nh1(-(f11-f01)/dx      ,-(f01+f11-2*fhh)/dy,1); int rgb_h1=Phong_Shading(Nh1,rgb);
			POINT3 N0h(-(2*fhh-f00-f01)/dx,-(f01-f00)/dy      ,1); int rgb_0h=Phong_Shading(N0h,rgb);
			POINT3 N1h(-(f10+f11-2*fhh)/dx,-(f11-f10)/dy      ,1); int rgb_1h=Phong_Shading(N1h,rgb);
			draw_triangle(Phh,P00,P10,rgb_h0);
			draw_triangle(Phh,P10,P11,rgb_1h);
			draw_triangle(Phh,P11,P01,rgb_h1);
			draw_triangle(Phh,P01,P00,rgb_0h);}}}
}


//
//void Display::contour( const Grid_2D& grid, const QuadTree& Tr, const Array<double>& f, int rgb, double level)
//{
//	int imax = grid.imax;
//	int jmax = grid.jmax;
//
//	// finding zmin and zmax
//	double zmin =  100000;
//	double zmax = -100000;
//	int num_of_nodes = Tr.num_of_nodes();
//
//	for(int n=0;n<num_of_nodes;n++)
//	{
//		double z = f(n);
//
//		zmin = MIN(z,zmin);
//		zmax = MAX(z,zmax);
//	}
//
//	// display region a little larger than the physical region
//	double xmin = grid.xmin, xmax = grid.xmax;
//	double ymin = grid.ymin, ymax = grid.ymax;
//
//	double xcen = (xmax+xmin)/2, xlen = (xmax-xmin)/2;
//	double ycen = (ymax+ymin)/2, ylen = (ymax-ymin)/2;
//	double zcen = (zmax+zmin)/2, zlen = (zmax-zmin)/2;
//
//	double len_max = MAX(MAX(xlen,ylen),zlen);
//
//	set_physical_coordinates(xcen-xlen*1.2,xcen+xlen*1.2,
//		                     ycen-ylen*1.2,ycen+ylen*1.2,
//							 zcen-zlen*1.2,zcen+zlen*1.2);
//
//	box(xcen-xlen,xcen+xlen,
//        ycen-ylen,ycen+ylen,
//        zcen-zlen,zcen+zlen);
//
//	// drawing
//	STACK<const QuadCell*> cells; cells.push(Tr.root);
//	STACK<LineSegment2> list_linesegments;
//
//	while(!cells.is_empty())
//	{
//		const QuadCell* cell = cells.pop();
//
//		if(cell->is_leaf())
//		{
//			const QuadLeaf* leaf = (const QuadLeaf*)cell;
//
//			POINT2 P00(grid.x_fr_i(leaf->imin  ),grid.y_fr_j(leaf->jmin  ));
//			POINT2 P01(grid.x_fr_i(leaf->imin  ),grid.y_fr_j(leaf->jmax()));
//			POINT2 P10(grid.x_fr_i(leaf->imax()),grid.y_fr_j(leaf->jmin  ));
//			POINT2 P11(grid.x_fr_i(leaf->imax()),grid.y_fr_j(leaf->jmax()));
//
//			double f00 = f(leaf->nodes[0][0]) - level;
//			double f01 = f(leaf->nodes[0][1]) - level;
//			double f10 = f(leaf->nodes[1][0]) - level;
//			double f11 = f(leaf->nodes[1][1]) - level;
//
//			Triangulation::do_in_2D(P00,f00,
//									P01,f01,
//									P10,f10,
//									P11,f11,list_linesegments);
//		}
//		else
//		{
//			const QuadBranch* branch = (const QuadBranch*)cell;
//
//			for(int a=0;a<2;a++)
//			for(int b=0;b<2;b++)
//				cells.push(branch->children[a][b]);
//		}
//	}
//
//
//	for(int n=0;n<list_linesegments.size();n++)
//	{
//		LineSegment2 lsg = list_linesegments(n);
//
//		draw_line(lsg.P1,lsg.P2,rgb);
//	}
//}

//
void Display::isosurface( const ARRAY3D<double>& f, int rgb, double level, bool wire_frame)
{
	set_physical_coordinates(0,f.isize-1,0,f.jsize-1,0,f.ksize-1);
	m_scaling.x = MIN(m_scaling.x,
		              m_scaling.y,
					  m_scaling.z);
	m_scaling.y = m_scaling.x;
	m_scaling.z = m_scaling.x;
	//box(0,f.isize-1,0,f.jsize-1,0,f.ksize-1);
	box(-.5,f.isize-.5,
		-.5,f.jsize-.5,
		-.5,f.ksize-.5);

	//	Isosurfacing based on TetraCubes
//	STACK<Triangle3> Ts;
//	STACK<   POINT3> Ns;
//	Triangulation::do_in_3D(grid,f,Ts,Ns,level);
//
//	for(int n=0;n<Ts.size();n++)
//	{
//		const Triangle3& T = Ts(n);
//		const    POINT3& N = Ns(n);
//
//		if(wire_frame)
//		{
//			draw_line(T.P1,T.P2,rgb);
//			draw_line(T.P1,T.P3,rgb);
//			draw_line(T.P2,T.P3,rgb);
//		}
//		else
//		{
//			int rgb_shading = Phong_Shading(N,rgb);
//			draw_triangle(T.P1,T.P2,T.P3,rgb_shading);
//		}
//	}
//	
	// Modified Marching Cubes
	int i,j,k;
//#ifndef WIN32
//#pragma omp parallel for private(i,j,k)
//#endif
	for(i=0;i<f.isize-1;i++)
	for(j=0;j<f.jsize-1;j++)
	for(k=0;k<f.ksize-1;k++)
	{
		// test를 해보자
		//if(i==3 && j==2 && k==3)
		{
			double phi[2][2][2];
			
			for(int a=0;a<2;a++)
			for(int b=0;b<2;b++)
			for(int c=0;c<2;c++)
				phi[a][b][c] = f[i+a][j+b][k+c] - level;
			
			double xmin = i; double xmax = i+1;
			double ymin = j; double ymax = j+1;
			double zmin = k; double zmax = k+1;
			
			STACK<POINT3> pts;
			
			MarchingCubes::isosurface( phi,
						            xmin, xmax,
						            ymin, ymax,
						            zmin, zmax, pts);
			
			int num_trs = pts.size()/3;
			
			for(int n=0;n<num_trs;n++)
			{
				// DEBUG
				static int jmax=0;
				if(j>jmax)
					jmax=j;
				const POINT3& P1 = pts(3*n  );
				const POINT3& P2 = pts(3*n+1);
				const POINT3& P3 = pts(3*n+2);
				
				POINT3 N; POINT3::cross(P2-P1,P3-P1,N);
				N/=N.abs();
			
				int rgb_shading = Phong_Shading(N,rgb);
				draw_triangle(P1,P2,P3,rgb_shading);
			}
		}
	}
}



/*
////
//void Display::graph( const Grid_2D& grid, const QuadTree& Tr, const Array<double>& F, int rgb, bool wire_frame)
//{
//	make_it_white();
//
//	int imax = grid.imax;
//	int jmax = grid.jmax;
//
//	// finding zmin and zmax
//	double zmin =  100000;
//	double zmax = -100000;
//
//	for(int n=0;n<F.size();n++)
//	{
//		double z = F(n);
//
//		zmin = MIN(z,zmin);
//		zmax = MAX(z,zmax);
//	}
//
//	// display region a little larger than the physical region
//	double xmin = grid.xmin, xmax = grid.xmax;
//	double ymin = grid.ymin, ymax = grid.ymax;
//
//	double xcen = (xmax+xmin)/2, xlen = (xmax-xmin)/2;
//	double ycen = (ymax+ymin)/2, ylen = (ymax-ymin)/2;
//	double zcen = (zmax+zmin)/2, zlen = (zmax-zmin)/2;
//
//	double len_max = MAX(MAX(xlen,ylen),zlen);
//
//	set_physical_coordinates(xcen-xlen*1.2,xcen+xlen*1.2,
//		                     ycen-ylen*1.2,ycen+ylen*1.2,
//							 zcen-zlen*1.2,zcen+zlen*1.2);
//
//	box(xcen-xlen,xcen+xlen,
//        ycen-ylen,ycen+ylen,
//        zcen-zlen,zcen+zlen);
//
//	// drawing
//	STACK<const QuadCell*> STACK;
//	STACK.push(Tr.root);
//
//	while(STACK.size()>0)
//	{
//		const QuadCell* cell = STACK.pop();
//
//		if(cell->is_leaf())
//		{
//			const QuadLeaf* leaf = (const QuadLeaf*)cell;
//
//			POINT3 P00(grid.x_fr_i(leaf->imin  ), grid.y_fr_j(leaf->jmin  ), F(leaf->nodes[0][0]));
//			POINT3 P01(grid.x_fr_i(leaf->imin  ), grid.y_fr_j(leaf->jmax()), F(leaf->nodes[0][1]));
//			POINT3 P10(grid.x_fr_i(leaf->imax()), grid.y_fr_j(leaf->jmin  ), F(leaf->nodes[1][0]));
//			POINT3 P11(grid.x_fr_i(leaf->imax()), grid.y_fr_j(leaf->jmax()), F(leaf->nodes[1][1]));
//
//			if(wire_frame)
//			{
//				draw_line(P00,P01,rgb);
//				draw_line(P00,P10,rgb);
//				draw_line(P11,P01,rgb);
//				draw_line(P11,P10,rgb);
//			}
//			else
//			{
//				POINT3 N00_10_11(-(F(leaf->nodes[1][0])-F(leaf->nodes[0][0]))/(grid.dx*leaf->size),
//					             -(F(leaf->nodes[1][1])-F(leaf->nodes[1][0]))/(grid.dy*leaf->size),1);
//				POINT3 N00_01_11(-(F(leaf->nodes[1][1])-F(leaf->nodes[0][1]))/(grid.dx*leaf->size),
//					             -(F(leaf->nodes[0][1])-F(leaf->nodes[0][0]))/(grid.dy*leaf->size),1);
//
//				int rgb_00_10_11 = Phong_Shading(N00_10_11,rgb);
//				int rgb_00_01_11 = Phong_Shading(N00_01_11,rgb);
//
//				draw_triangle(P00,P10,P11,rgb_00_10_11);
//				draw_triangle(P00,P01,P11,rgb_00_01_11);
//			}	
//		}
//		else
//		{
//			const QuadBranch* branch = (const QuadBranch*)cell;
//
//			for(int a=0;a<2;a++)
//			for(int b=0;b<2;b++)
//				STACK.push(branch->children[a][b]);
//		}
//	}
//}

void Display::vf2( const GRID_2D& grid,
				   const ARRAY2D<double>& U,
				   const ARRAY2D<double>& V, int rgb)
{
	int imax = grid.isize-1;
	int jmax = grid.jsize-1;

	// display region a little larger than the physical region
	double xmin = grid.xmin, xmax = grid.xmax; double dx=grid.dx;
	double ymin = grid.ymin, ymax = grid.ymax; double dy=grid.dy;

	double xcen = (xmax+xmin)/2, xlen = (xmax-xmin)/2;
	double ycen = (ymax+ymin)/2, ylen = (ymax-ymin)/2;

	double len_max = MAX(xlen,ylen);

	set_physical_coordinates(xcen-xlen*1.2,xcen+xlen*1.2,
		                     ycen-ylen*1.2,ycen+ylen*1.2,
							 ycen-ylen*1.2,ycen+ylen*1.2);
	//
	double max_magnitude = 1E-9;

	for(int i=0;i<grid.isize;i++)
	for(int j=0;j<grid.jsize;j++)
	{
		double u = U[i][j];
		double v = V[i][j];
		double abs = sqrt(u*u+v*v); if(abs>max_magnitude)max_magnitude = abs;
	}

	double scale = MIN(grid.dx,grid.dy)/max_magnitude;

	for(int i=0;i<grid.isize;i++)
	for(int j=0;j<grid.jsize;j++)
	{
		double x = grid.x_fr_i(i);
		double y = grid.y_fr_j(j);
		double u = U[i][j];
		double v = V[i][j];
		draw_line(x,y,x+u*scale,y+v*scale,rgb);

		if(i<grid.isize               ) draw_line(x   ,y   ,x+dx,y   ,0);
		if(j<grid.jsize               ) draw_line(x   ,y   ,x   ,y+dy,0);
		if(i<grid.isize&& j<grid.jsize) draw_line(x+dx,y+dy,x+dx,y   ,0);
		if(i<grid.isize&& j<grid.jsize) draw_line(x+dx,y+dy,x+dx,y+dy,0);
	}
}

/*
void Display::vf3( const GRID_3D& grid,
				   const ARRAY3D<double>& U,
				   const ARRAY3D<double>& V,
				   const ARRAY3D<double>& W, int rgb)
{
	int imax = grid.imax;
	int jmax = grid.jmax;
	int kmax = grid.kmax;

	//
	double max_magnitude = 1E-9;

	for(int i=0;i<=grid.imax;i++)
	for(int j=0;j<=grid.jmax;j++)
	for(int k=0;k<=grid.kmax;k++)
	{
		double u = U(i,j,k);
		double v = V(i,j,k);
		double w = W(i,j,k);
		double abs = sqrt(u*u+v*v+w*w); if(abs>max_magnitude)max_magnitude = abs;
	}

	double scale = MIN(MIN(grid.dx,grid.dy),grid.dz)/max_magnitude;

	for(int i=0;i<=grid.imax;i++)
	for(int j=0;j<=grid.jmax;j++)
	for(int k=0;k<=grid.kmax;k++)
	{
		double x = grid.x_fr_i(i);
		double y = grid.y_fr_j(j);
		double z = grid.z_fr_k(k);
		double u = U(i,j,k);
		double v = V(i,j,k);
		double w = W(i,j,k);
		draw_line( POINT3(x,y,z),
				   POINT3( x+u*scale,
						   y+v*scale,
						   z+w*scale), rgb);
	}
}
*/
