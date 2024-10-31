#ifndef __MY_WINDOW__
#define __MY_WINDOW__

#include <windows.h>
#include "Display.h"
#include "../library.h"

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
	static bool	_mode_df3_graph_xcenter ;
	static bool	_mode_df3_graph_ycenter ;
	static bool	_mode_df3_graph_zcenter ;
	//static  Grid_1D         _grid_1d;
	static  GRID_2D         _grid_2d;
	static  GRID_3D         _grid_3d;
	//static Array   <double>    _f_1d;
	static ARRAY_2D<double>    _f_2d;
	static ARRAY_3D<double>    _f_3d;
	//static Array_3D<double>    _vf_3d_u;
	//static Array_3D<double>    _vf_3d_v;
	//static Array_3D<double>    _vf_3d_w;
	
public:
	//static Point2 mouse_click; 
	static Display m_display;
	static HWND    m_hwd;
	static unsigned char* m_rgb_buf;
	
	static void display();
	static void render();
};

#endif
