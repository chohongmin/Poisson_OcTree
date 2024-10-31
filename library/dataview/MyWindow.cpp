#include <math.h>
#include <string.h>
#include "MyWindow.h"
#include "resource.h"

void _RGB_SPLIT(int rgb, unsigned char& r,
				         unsigned char& g,
				         unsigned char& b)
{
	b = rgb & 0x0000ff; rgb>>=8;
	g = rgb & 0x0000ff; rgb>>=8;
	r = rgb           ;
}

// static members
Display        MyWindow::m_display;
HWND           MyWindow::m_hwd;
unsigned char* MyWindow::m_rgb_buf;

bool MyWindow::_mode_df1 = false;
bool MyWindow::_mode_df2 = false;
bool MyWindow::_mode_df3 = false;
bool MyWindow::_mode_vf3 = false;
bool MyWindow::_mode_df2_contour = false;
bool MyWindow::_mode_df2_graph   = false;
bool MyWindow::_mode_df3_graph   = false;
bool MyWindow::_mode_df3_isosurface_zero = false;
bool MyWindow::_mode_df3_isosurface_plus_dx = false;
bool MyWindow::_mode_df3_isosurface_minus_dx= false;
bool MyWindow::_mode_df3_graph_xcenter = false;
bool MyWindow::_mode_df3_graph_ycenter = false;
bool MyWindow::_mode_df3_graph_zcenter = false;


ARRAY  <double>    MyWindow::_f_1d;
ARRAY2D<double>    MyWindow::_f_2d;
ARRAY3D<double>    MyWindow::_f_3d;

// select drawing mode
void MyWindow::unselect_all_drawing_modes()
{
	 _mode_df1 = false;
	 _mode_df2 = false;
	 _mode_df3 = false;
	 _mode_df2_contour = false;
	 _mode_df2_graph   = false;
	 _mode_df3_isosurface_zero = false;
	 _mode_df3_isosurface_plus_dx = false;
	 _mode_df3_isosurface_minus_dx= false;
	 _mode_df3_graph         = false;
	 _mode_df3_graph_xcenter = false;
	 _mode_df3_graph_ycenter = false;
	 _mode_df3_graph_zcenter = false;
	 _mode_vf3 = false;
}

void MyWindow::_2d_contour          () {unselect_all_drawing_modes(); _mode_df2 = true; _mode_df2_contour=true; display();}
void MyWindow::_2d_graph            () {unselect_all_drawing_modes(); _mode_df2 = true; _mode_df2_graph  =true; display();}
void MyWindow::_3d_graph            () {unselect_all_drawing_modes(); _mode_df3 = true; _mode_df3_graph  =true; display();}
void MyWindow::_3d_isosurface_zero()         {unselect_all_drawing_modes(); _mode_df3 = true; _mode_df3_isosurface_zero=true    ; display();}
void MyWindow::_3d_isosurface_plus_dx ()     {unselect_all_drawing_modes(); _mode_df3 = true; _mode_df3_isosurface_plus_dx=true ; display();}
void MyWindow::_3d_isosurface_minus_dx()     {unselect_all_drawing_modes(); _mode_df3 = true; _mode_df3_isosurface_minus_dx=true; display();}
void MyWindow::_3d_graph_xcenter(){unselect_all_drawing_modes(); _mode_df3 = true; _mode_df3_graph_xcenter=true; display();}
void MyWindow::_3d_graph_ycenter(){unselect_all_drawing_modes(); _mode_df3 = true; _mode_df3_graph_ycenter=true; display();}
void MyWindow::_3d_graph_zcenter(){unselect_all_drawing_modes(); _mode_df3 = true; _mode_df3_graph_zcenter=true; display();}
//void MyWindow::_3d_vector_field (){unselect_all_drawing_modes(); _mode_vf3 = true; display();}


void MyWindow::render()
{
	m_display.make_it_white();
	
	if(_mode_df1                                 ) m_display.graph       (_f_1d);
	if(_mode_df2 && _mode_df2_contour            ) m_display.contour     (_f_2d,0xff0000);
	if(_mode_df2 && _mode_df2_graph              ) m_display.graph       (_f_2d,0x0000ff,false);
	if(_mode_df3 && _mode_df3_isosurface_zero    ) m_display.isosurface  (_f_3d,0xff0000);
	if(_mode_df3 && _mode_df3_isosurface_plus_dx ) m_display.isosurface  (_f_3d,0xff0000, .01);
	if(_mode_df3 && _mode_df3_isosurface_minus_dx) m_display.isosurface  (_f_3d,0xff0000,-.01);
	if(_mode_df3 && _mode_df3_graph_xcenter      ) m_display.graph_xconst(_f_3d,0xff0000,false);
	if(_mode_df3 && _mode_df3_graph_ycenter      ) m_display.graph_yconst(_f_3d,0xff0000,false);
	if(_mode_df3 && _mode_df3_graph_zcenter      ) m_display.graph_zconst(_f_3d,0xff0000,false);
	if(_mode_df3 && _mode_df3_graph              ) m_display.graph       (_f_3d );
	//if(_mode_vf3                                 ) m_display.vf3         (_grid_3d,_vf_3d_u,
	//	                                                                           _vf_3d_v,
	//																	           _vf_3d_w);
}


void MyWindow::display()
{
	render();
	
	//
	int w = m_display.m_imax+1;
	int h = m_display.m_jmax+1;
	
	//
	BITMAPINFO bminfo;
	BITMAPINFOHEADER& bmInfoHeader = bminfo.bmiHeader;

    bmInfoHeader.biSize    = 40;
	bmInfoHeader.biWidth    = m_display.m_imax+1;
    bmInfoHeader.biHeight   = m_display.m_jmax+1;
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
	unsigned char* buf = m_rgb_buf;
	for (int j=0; j < h; j++)
	for (int i=0; i < w; i++)
	{
		int rgb = m_display.m_rgb[i][j];
		unsigned char r,g,b; _RGB_SPLIT(rgb,r,g,b);
		*buf = b; buf++;
		*buf = g; buf++;
		*buf = r; buf++;
	}

	//
	HDC hdc = GetDC(m_hwd);
	SetDIBitsToDevice(hdc,0,0,w,h,0,0,0,h,m_rgb_buf,&bminfo,DIB_RGB_COLORS);
	ReleaseDC(m_hwd,hdc);
}

LRESULT CALLBACK WndProc(HWND,UINT,WPARAM,LPARAM);

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance,
				   PSTR szCmdLine, int iCmdShow)
{
	static TCHAR szAppName[] = TEXT("DataView");
	HWND hwnd;
	MSG  msg;
	WNDCLASS wndclass;

	wndclass.style = CS_HREDRAW | CS_VREDRAW;
	wndclass.lpfnWndProc = WndProc;
	wndclass.cbClsExtra =0;
	wndclass.cbWndExtra = 0;
	wndclass.hInstance = hInstance;
	wndclass.hIcon = LoadIcon(NULL,IDI_APPLICATION);
	wndclass.hCursor= LoadCursor(NULL,IDC_ARROW);
	wndclass.hbrBackground = (HBRUSH) GetStockObject(WHITE_BRUSH);
	wndclass.lpszMenuName = szAppName;
	wndclass.lpszClassName = szAppName;

	RegisterClass(&wndclass);
	hwnd = CreateWindow(szAppName,szAppName,WS_OVERLAPPEDWINDOW,
		CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT,
		NULL,NULL,hInstance, NULL);
	MyWindow::m_hwd = hwnd;
	ShowWindow(hwnd,iCmdShow);
	UpdateWindow(hwnd);

	HACCEL hAccel = LoadAccelerators(hInstance,szAppName);

	while(GetMessage(&msg,NULL,0,0))
	{
		if(!TranslateAccelerator(hwnd,hAccel,&msg))
		{
			TranslateMessage(&msg);
			DispatchMessage (&msg);
		}
	}

	return msg.wParam;
}

//
LRESULT CALLBACK WndProc( HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam )
{
	HDC			hdc;
	PAINTSTRUCT ps;
	
	switch(message)
	{
	case WM_CREATE:
		{
			int w = MyWindow::m_display.m_imax+1;
			int h = MyWindow::m_display.m_jmax+1;
			MyWindow::m_rgb_buf = new unsigned char[w*h*3];
			return 0;
		}

	case WM_DESTROY:
		delete[] MyWindow::m_rgb_buf;
		PostQuitMessage(0);
		return 0;
			
	case WM_COMMAND:
		{
			HMENU hMenu = GetMenu(hwnd);

			switch( LOWORD(wParam))
			{
			case IDM_FILE_OPEN :
				{
					TCHAR szFilter[] =	  TEXT("Data Files\0*.df1;*.df2;*.df3;*.vf3\0");

					TCHAR FileName[500];
					TCHAR TitleName[500];

					OPENFILENAME ofn;
					ZeroMemory(&ofn, sizeof(ofn));
					ofn.lStructSize = sizeof(ofn);
					ofn.hwndOwner = hwnd;
					ofn.lpstrFile = FileName;
					ofn.lpstrFile[0] = '\0';
					ofn.nMaxFile = sizeof(FileName);
					ofn.lpstrFilter = szFilter;
					ofn.nFilterIndex = 1;
					ofn.lpstrFileTitle = TitleName;
					ofn.nMaxFileTitle = sizeof(TitleName);
					ofn.lpstrInitialDir = NULL;
					ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

					if(GetOpenFileName(&ofn))
					{
						char filename [500];
						char titlename[500];
						WideCharToMultiByte(CP_ACP, 0,  FileName, -1,  filename, 499, NULL, NULL);
						WideCharToMultiByte(CP_ACP, 0, TitleName, -1, titlename, 499, NULL, NULL);

						SetWindowText(hwnd, TitleName);
						
						int len=strlen(titlename);
						char ext[4]; ext[0] = titlename[len-3];
									 ext[1] = titlename[len-2];
									 ext[2] = titlename[len-1];
									 ext[3] = 0;

						MyWindow::unselect_all_drawing_modes();

						// vf3 
						//if(strcmp(ext,"vf3")==0)
						//{
						//	MyWindow::_mode_vf3 = true;
						//	FileIO::load(MyWindow::_grid_3d,MyWindow::_vf_3d_u,
						//									MyWindow::_vf_3d_v,
						//									MyWindow::_vf_3d_w, filename);
						//	MyWindow::m_display.set_physical_coordinates(MyWindow::_grid_3d.xmin, MyWindow::_grid_3d.xmax,
						//												 MyWindow::_grid_3d.ymin, MyWindow::_grid_3d.ymax,
						//												 MyWindow::_grid_3d.zmin, MyWindow::_grid_3d.zmax);
						//}

						// df3 
						if(strcmp(ext,"df3")==0)
						{
							MyWindow::_mode_df3 = true;
							MyWindow::_mode_df3_isosurface_zero = true;
							MyWindow::_f_3d.load(filename);
							MyWindow::display();
						}

						// df2 
						if(strcmp(ext,"df2")==0)
						{
							MyWindow::_mode_df2 = true;
							MyWindow::_mode_df2_graph = true;
							MyWindow::_f_2d.load(filename);
							MyWindow::display();
						}

						// df1 
						if(strcmp(ext,"df1")==0)
						{
							MyWindow::_mode_df1 = true;
							MyWindow::_f_1d.load(filename);
							MyWindow::display();
						}
					}
				}
				break;
			case IDM_2D_GRAPH  :MyWindow::_2d_graph(); break;
			case IDM_2D_CONTOUR :MyWindow::_2d_contour(); break;
			case IDM_3D_ISOSURFACE_ZERO :MyWindow::_3d_isosurface_zero(); break;
			case IDM_3D_ISOSURFACE_DX_P :MyWindow::_3d_isosurface_plus_dx(); break;
			case IDM_3D_ISOSURFACE_DX_M :MyWindow::_3d_isosurface_minus_dx(); break;
			case IDM_3D_GRAPH         :MyWindow::_3d_graph        (); break;	
			case IDM_3D_GRAPH_XCENTER :MyWindow::_3d_graph_xcenter(); break;
			case IDM_3D_GRAPH_YCENTER :MyWindow::_3d_graph_ycenter(); break;
			case IDM_3D_GRAPH_ZCENTER :MyWindow::_3d_graph_zcenter(); break;
			}
			return 0;
		}

	case WM_PAINT:
		{
			hdc = BeginPaint(hwnd,&ps);
			MyWindow::render();
			//
			int w = MyWindow::m_display.m_imax+1;
			int h = MyWindow::m_display.m_jmax+1;
			
			//
			BITMAPINFO bminfo;
			BITMAPINFOHEADER& bmInfoHeader = bminfo.bmiHeader;

			bmInfoHeader.biSize    = 40;
			bmInfoHeader.biWidth    = MyWindow::m_display.m_imax+1;
			bmInfoHeader.biHeight   = MyWindow::m_display.m_jmax+1;
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
			unsigned char* buf = MyWindow::m_rgb_buf;
				
			for (int j=0; j < h; j++)
			for (int i=0; i < w; i++)
			{
				int rgb = MyWindow::m_display.m_rgb[i][j];
				unsigned char r,g,b; _RGB_SPLIT(rgb,r,g,b);
				*buf = b; buf++;
				*buf = g; buf++;
				*buf = r; buf++;
			}

			//
			SetDIBitsToDevice(hdc,0,0,w,h,0,0,0,h,MyWindow::m_rgb_buf,&bminfo,DIB_RGB_COLORS);
			
			EndPaint(hwnd,&ps);
			
			return 0;
		}

	case WM_SIZE:
		{
			RECT rect_client; GetClientRect(hwnd,&rect_client);
			RECT rect_window; GetWindowRect(hwnd,&rect_window);

			int width = rect_window.right - rect_window.left 
					  - rect_client.right + rect_client.left;

			int height = rect_window.bottom - rect_window.top 
					   - rect_client.bottom + rect_client.top;

			MoveWindow(hwnd,5,5,width +MyWindow::m_display.m_imax,
								height+MyWindow::m_display.m_jmax,TRUE);
			return 0;
		}

	case WM_KEYDOWN:
		switch( wParam )
		{
		case VK_PRIOR :
			MyWindow::m_display.m_scaling.x *=1.1;
			MyWindow::m_display.m_scaling.y *=1.1;
			MyWindow::m_display.m_scaling.z *=1.1;
			break;

		case VK_NEXT :
			MyWindow::m_display.m_scaling.x /=1.1;
			MyWindow::m_display.m_scaling.y /=1.1;
			MyWindow::m_display.m_scaling.z /=1.1;
			break;

		case VK_LEFT : MyWindow::m_display.rotate_to_y_minus(); break;
		case VK_RIGHT: MyWindow::m_display.rotate_to_y_plus (); break;
		case VK_UP   : MyWindow::m_display.rotate_to_x_minus(); break;
		case VK_DOWN : MyWindow::m_display.rotate_to_x_plus (); break;
		};
		MyWindow::display();

	case WM_CHAR:
		switch( wParam)
		{
		case 'z' : MyWindow::m_display.m_rotation.set_rotation(      0,1,0,0); break;
		case 'x' : MyWindow::m_display.m_rotation.set_rotation(-2*PI/3,1,1,1); break;
		case 'y' : MyWindow::m_display.m_rotation.set_rotation( 2*PI/3,1,1,1); break;
		case 'o' : MyWindow::m_display.m_rotation.set_rotation(  -PI/2,1,0,0); break;

		case 'l': MyWindow::m_display.m_translation.x -= .1; break;
		case 'r': MyWindow::m_display.m_translation.x += .1; break;
		case 'b': MyWindow::m_display.m_translation.y -= .1; break;
		case 't': MyWindow::m_display.m_translation.y += .1; break;

		case 'R' :
			printf("the currect rotation : %f + (%f,%f,%f)\n", MyWindow::m_display.m_rotation.a,
				MyWindow::m_display.m_rotation.u1,
				MyWindow::m_display.m_rotation.u2,
				MyWindow::m_display.m_rotation.u3);
			break;

		case 's':
			MyWindow::m_display.save_bmp(); break;
		};

		MyWindow::display();

	}

	return DefWindowProc(hwnd,message,wParam,lParam);
}