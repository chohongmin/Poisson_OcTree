#pragma once
#include "./library/STACK.h"
#include "./library/ARRAY.h"
#include "library/visualization/Display.h"


struct OcFace{
	int imin,jmin,kmin,imax,jmax,kmax;
	int index_prnt;
	int index_child[2][2];
	int index_cell[2];
	bool leaf;

	OcFace();
	void set_child(int a, int b, OcFace& child) const;
	int size() const{return (imax-imin+jmax-jmin+kmax-kmin)/2;};
	int type_ijk() const;
};

struct OcCell{
	int  imin,jmin,kmin,imax,jmax,kmax;
	int  index_prnt;
	int  index_child[2][2][2];
	int  index_face [3][2];
	bool leaf;
	
	OcCell();
	void set_child(int i, int j, int k, OcCell& child) const;
	void set_face(int type, int l, OcFace& F) const;
	int size() const{return imax-imin;};
};



struct OcTree{
	STACK<OcCell> cells;
	STACK<OcFace> faces;
	int N;
	// tree generation
	int num_cells() const{return cells.size();};
	int num_faces() const{return faces.size();};
	OcTree(int N);
	void split_c(int index_c);
	void split_f(int index_f);
	int search_leafcell(double i, double j, double k) const;
	void set_face_cell_relation(int index_f);
	void split_all_cells();
	void generate_ot(double(*phi)(double, double,double), double Lip, int Csize_min = 1);
	void visualize_domain( double (*phi)(double,double,double), Display& disp);
	void visualize_xsection(Display& disp);
	void visualize_ysection(Display& disp);
	void visualize_zsection(Display& disp);
	void visualize_ysection(double (*phi)(double, double, double), const ARRAY<double>& U, double scale);
	void visualize_zsection(double (*phi)(double, double, double), const ARRAY<double>& U, double scale);
	// physical coordinates
	double xmin, xmax, ymin, ymax, zmin, zmax;
	double x_fr_i(double i) const { return (xmax - xmin) / N * i + xmin; }
	double y_fr_j(double j) const { return (ymax - ymin) / N * j + ymin; }
	double z_fr_k(double k) const { return (zmax - zmin) / N * k + zmin; }

	double i_fr_x(double x) const { return N * (x - xmin) / (xmax - xmin); }
	double j_fr_y(double y) const { return N * (y - ymin) / (ymax - ymin); }
	double k_fr_z(double z) const { return N * (z - zmin) / (zmax - zmin); }

	double xsize(int index_c)const { return (xmax - xmin) / N * (cells[index_c].imax - cells[index_c].imin); }
	double ysize(int index_c)const { return (ymax - ymin) / N * (cells[index_c].jmax - cells[index_c].jmin); }
	double zsize(int index_c)const { return (zmax - zmin) / N * (cells[index_c].kmax - cells[index_c].kmin); }

	double xcen(int index_c)const { return x_fr_i((cells[index_c].imin + cells[index_c].imax) / 2.); }
	double ycen(int index_c)const { return y_fr_j((cells[index_c].jmin + cells[index_c].jmax) / 2.); }
	double zcen(int index_c)const { return z_fr_k((cells[index_c].kmin + cells[index_c].kmax) / 2.); }

	// Laplacian on Rectangle
	double calculate_G(const ARRAY<double>& u, int index_f) const;
	void   calculate_G(const ARRAY<double>& u, ARRAY<double>& Gu) const;
	double calculate_flux(const ARRAY<double>& f, int index_f) const;
	void calculate_D(const ARRAY<double>& f, ARRAY<double>& D) const;
	void calculate_W(const ARRAY<double>& f, int index_F, double& Wf, double& delta) const;
	void calculate_W(const ARRAY<double>& f, ARRAY<double>& Wf) const;
	void set_W(int index_F, double Wf_value, ARRAY<double>& Wf) const;
	double inner_product(const ARRAY<double>& x, const ARRAY<double>& y) const;
	void projection_sumzero(ARRAY<double>& x) const;
	void solve_Laplacian(const ARRAY<double>& f, ARRAY<double>& u)const;

	// Irregular domains
	void calculate_H( double (*phi)(double,double,double), ARRAY<double>& H)const;
	void calculate_H( const ARRAY<double>& f, const ARRAY<double>& H, ARRAY<double>& Hf) const;
	void solve_Laplacian(const ARRAY<double>& H, const ARRAY<double>& f, ARRAY<double>& u)const;

	// 
	double face_Area (int index_F)const;
	double face_Delta(int index_F)const;
	// FSI
	void calculate_GH_J(const ARRAY<double>& WH, ARRAY<double>& GH_x, 
		                                         ARRAY<double>& GH_y, 
		                                         ARRAY<double>& GH_z, ARRAY<double>& Jx,
		                                                              ARRAY<double>& Jy,
		                                                              ARRAY<double>& Jz, const double cx, 
		                                                                                 const double cy,
		                                                                                 const double cz) const;
	void solve_FSI    ( const ARRAY<double>& H, const ARRAY<double>& RHS, 
		           double minv, double Iinv[3][3], double cx, double cy, double cz, ARRAY<double>& pressure) const;
	/*
	// Matrix A=DWHG generation
	void calculate_Gmatrix(int index_f, double scale, SparseSymmMatrix& A, int i);
	double calculate_WDelta(int index_f) const;
	void calculate_WHGmatrix(int index_f, double scale, const ARRAY<double>& H, SparseSymmMatrix& A, int i);
	void calculate_DWHGmatrix(const ARRAY<double>& H, SparseSymmMatrix& A);
	*/


};

