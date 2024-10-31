#ifndef DISPLAY_H
#define DISPLAY_H

#include "../ARRAY.h"
#include "../ARRAY2D.h"
#include "../ARRAY3D.h"
#include "../algebra/QUATERNION.h"


/*----------------------------------------------------------------*\
 *
 * Display structure, where we can display functions and their
 * contours. Chohong Min, 2008, chohong@khu.ac.kr
 *
\*----------------------------------------------------------------*/
class Display
{
	public:
	/*----------------------------------------------------------------*\
	 * colours and z-buffers
	\*----------------------------------------------------------------*/
	int m_imax; // i=0...imax
	int m_jmax; // j=0...jmax

	ARRAY2D< int    > m_rgb;
	ARRAY2D< double > m_z;

	//
	Display(int isize=600, int jsize=600)
	{
		m_imax = isize-1;
		m_jmax = jsize-1;
		
		m_rgb.resize(isize,jsize);
		m_z  .resize(isize,jsize); 
		m_rotation.set_rotation(0,1,0,0);
	}
	//
	void make_it_white()
	{
		m_rgb = 0xffffff;
		m_z =-1000;
	}
	//
	void save_bmp(const char* file_name=0) ;
	
	/*----------------------------------------------------------------*\
	 * mapping of logical and physical coordinates
	 * given (x,y,z),
	 *  	apply translation, scaling, rotation,
	 *      then  draw [-1,1]^3
	\*----------------------------------------------------------------*/
	POINT3      m_translation;
	POINT3      m_scaling;
	QUATERNION m_rotation;
	
	//
	void set_physical_coordinates(double xmin=-1, double xmax=1, 
                                  double ymin=-1, double ymax=1,
								  double zmin=-1, double zmax=1);

	//
    void from_phys_to_logi_coords(const POINT3& P, double& i, 
		                                           double& j, 
												   double& z) const;    

    //
    void from_logi_to_phys_coords(double i, 
		                          double j, 
								  double z, POINT3& P) const;
    
    //
    void rotate_to_x_plus (){ QUATERNION Q; Q.set_rotation( 0.11,1,0,0); m_rotation = Q*m_rotation; }
    void rotate_to_x_minus(){ QUATERNION Q; Q.set_rotation(-0.11,1,0,0); m_rotation = Q*m_rotation; }
    void rotate_to_y_plus (){ QUATERNION Q; Q.set_rotation( 0.11,0,1,0); m_rotation = Q*m_rotation; }
    void rotate_to_y_minus(){ QUATERNION Q; Q.set_rotation(-0.11,0,1,0); m_rotation = Q*m_rotation; }
    void rotate_to_z_plus (){ QUATERNION Q; Q.set_rotation( 0.11,0,0,1); m_rotation = Q*m_rotation; }
    void rotate_to_z_minus(){ QUATERNION Q; Q.set_rotation(-0.11,0,0,1); m_rotation = Q*m_rotation; }

	/*----------------------------------------------------------------*\
	 * basic drawings
	\*----------------------------------------------------------------*/
	void POINT( int i, int j, double z, int rgb );
    void draw_line( double x0, double y0,
    				double x1, double y1, int rgb );

    void draw_pixel( const POINT3& P, int rgb, int w=1);
    
	void draw_line( const POINT3& P1,
		            const POINT3& P2, int rgb );

	void draw_triangle( const POINT3& P1,
		                const POINT3& P2,
						const POINT3& P3, int rgb);

	void draw_curve(  double x0, double y0,double a0, double b0,
					  double x1, double y1,double a1, double b1, int rgb);

	void draw_number_fr_0_to_9( int number, double xmin, double xmax,
									        double ymin, double ymax, int rgb);

	void draw_signum( int sgn, double xmin, double xmax,
							   double ymin, double ymax, int rgb);

	void draw_number(  int number, double xmin, double xmax,
						           double ymin, double ymax, int rgb);
	
	void draw_number( double number, double xmin, double xmax,
			                         double ymin, double ymax, int rgb);

	void draw_number_hexa(  int number, double xmin, double xmax,
						                double ymin, double ymax, int rgb);

	int Phong_Shading( POINT3 N, int rgb ) const;

	/*----------------------------------------------------------------*\
	 * Drawings of functions
	\*----------------------------------------------------------------*/
	void graph       ( const ARRAY  <double>& f, int rgb=0);
	void graph       ( const ARRAY2D<double>& f, int rgb=0,bool wire_frame=false );
	void graph       ( const ARRAY3D<double>& f );
	void graph_xconst( const ARRAY3D<double>& f, int rgb=0,bool wire_frame=true );
	void graph_yconst( const ARRAY3D<double>& f, int rgb=0,bool wire_frame=true );
	void graph_zconst( const ARRAY3D<double>& f, int rgb=0,bool wire_frame=true );
	void contour     ( const ARRAY2D<double>& f, int rgb=0,double level=0);
	void isosurface  ( const ARRAY3D<double>& f, int rgb=0xff0000,double level=0,bool wire_frame=false );

//	void contour     ( const Grid_2D& grid, const QuadTree& Tr, const Array<double>& f, int rgb=0, double level=0);
//	void graph       ( const Grid_2D& grid, const QuadTree& Tr, const Array<double>& f, int rgb=0, bool wire_frame=true );
	void box(double xmin, double xmax,
		     double ymin, double ymax,
			 double zmin, double zmax, int rgb=0x000000);
//	void vf2( const GRID_2D& grid, const ARRAY2D<double>& U, 
//			                       const ARRAY2D<double>& V, int rgb=0);
	//void vf3( const GRID_3D& grid, const ARRAY3D<double>& U, 
	//		                       const ARRAY3D<double>& V,
	//							   const ARRAY3D<double>& W, int rgb=0);							   
};

#endif
