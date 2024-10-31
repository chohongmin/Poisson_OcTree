#ifndef __MY_WINDOW__
#define __MY_WINDOW__

#include <windows.h>
#include "../visualization/Display.h"
#include "../ARRAY.h"
#include "../ARRAY2D.h"
#include "../ARRAY3D.h"

/*-----------------------------------------------------------------------------
 *
 *  GTK Visualization Window
 *
 * File Extensions : "*.df2" for 2D discrete functions
 *                   "*.df3" for 3D discrete functions
 *                   "*.vf3" for 3D vector fields
 *
 *---------------------------------------------------------------------------*/
class MyWindow 
{
public:
	static void unselect_all_drawing_modes();
	static void _2d_contour();
	static void _2d_graph  ();
	static void _3d_isosurface_zero();
	static void _3d_isosurface_plus_dx ();
	static void _3d_isosurface_minus_dx();
	static void _3d_graph();
	static void _3d_graph_xcenter();
	static void _3d_graph_ycenter();
	static void _3d_graph_zcenter();
	//static void _3d_vector_field();
	
	static bool _mode_df1 ;
	static bool _mode_df2 ;
	static bool _mode_df3 ;
	static bool _mode_vf3 ;
	static bool	_mode_df2_contour;
	static bool	_mode_df2_graph  ;
	static bool	_mode_df3_isosurface_zero ;
	static bool	_mode_df3_isosurface_plus_dx ;
	static bool	_mode_df3_isosurface_minus_dx;
	static bool	_mode_df3_graph  ;
	static bool	_mode_df3_graph_xcenter ;
	static bool	_mode_df3_graph_ycenter ;
	static bool	_mode_df3_graph_zcenter ;
	static ARRAY  <double>    _f_1d;
	static ARRAY2D<double>    _f_2d;
	static ARRAY3D<double>    _f_3d;
	
public:
	//static Point2 mouse_click; 
	static Display m_display;
	static HWND    m_hwd;
	static unsigned char* m_rgb_buf;
	
	static void display();
	static void render();
};

#endif
