#include "OcTree.h"
#include <math.h>

#define PI 3.14159265358979323846264338328

double     u(double x, double y, double z) { return  (-2*y-2*z)* exp(-x-y*y-z*z); }
double     v(double x, double y, double z) { return  exp(-x-y*y-z*z); }
double     w(double x, double y, double z) { return  exp(-x-y*y-z*z); }
double ustar(double x, double y, double z) { return u(x,y,z)-sin(x)*cos(y)*cos(z); }
double vstar(double x, double y, double z) { return v(x,y,z)-cos(x)*sin(y)*cos(z); }
double wstar(double x, double y, double z) { return w(x,y,z)-cos(x)*cos(y)*sin(z); }

int main() {
	//----------------------------------------------
	// generation of quadtree
	//----------------------------------------------
	OcTree ot(1024);
	ot.split_c(0);
	ot.split_c(8);
	ot.split_c(4);
	ot.split_all_cells();
	ot.split_all_cells();
	ot.split_all_cells();
	ot.split_all_cells();

	//ot.visualize();
	ot.xmin = 0, ot.xmax = PI;
	ot.ymin = 0, ot.ymax = PI;
	ot.zmin = 0, ot.zmax = PI;

	ARRAY<double>   Ustar(ot.faces.size());
	ARRAY<double>  WUstar(ot.faces.size());
	ARRAY<double> DWUstar(ot.cells.size());
	for(int index_f=0; index_f<ot.num_faces(); index_f++){
		const OcFace& F = ot.faces[index_f];
		double x = ot.x_fr_i((F.imin+F.imax)/2.);
		double y = ot.y_fr_j((F.jmin+F.jmax)/2.);
		double z = ot.z_fr_k((F.kmin+F.kmax)/2.);
		     if(F.type_ijk()==0) Ustar[index_f] = ustar(x,y,z);
		else if(F.type_ijk()==1) Ustar[index_f] = vstar(x,y,z);
		else                     Ustar[index_f] = wstar(x,y,z);
	}
	ot.calculate_W(Ustar,WUstar);
	ot.calculate_D(WUstar,DWUstar);

	//----------------------------------------------
	// adjoint test
	//----------------------------------------------
	ARRAY<double> f1(ot.cells.size());
	ARRAY<double> f2(ot.cells.size());
	for (int i = 0; i < ot.cells.size(); i++) {
		f1[i] = ((double)rand()) / RAND_MAX;
		f2[i] = ((double)rand()) / RAND_MAX;
	}
	ARRAY<double> Gf1(ot.faces.size()); ot.calculate_G(f1, Gf1);
	ARRAY<double> Gf2(ot.faces.size()); ot.calculate_G(f2, Gf2);
	ARRAY<double> DGf1(ot.cells.size()); ot.calculate_D(Gf1, DGf1);
	ARRAY<double> DGf2(ot.cells.size()); ot.calculate_D(Gf2, DGf2);

	double f1_DGf2 = ot.inner_product(f1, DGf2);
	double DGf1_f2 = ot.inner_product(DGf1, f2);

	
	//----------------------------------------------
	// calculation of p through linear system DbetaGp=DU*
	//----------------------------------------------
	ARRAY<double> p(ot.cells.size());
	ot.solve_Laplacian(DWUstar, p);

	//----------------------------------------------
	// calculation of U = U*-WGp
	//----------------------------------------------
	ARRAY<double>  Gp(ot.faces.size());
	ARRAY<double> WGp(ot.faces.size());
	ARRAY<double>   U(ot.faces.size());
	ot.calculate_G( p,  Gp);
	ot.calculate_W(Gp, WGp);

	for (int index_f = 0; index_f < ot.faces.size(); index_f++) {
		const OcFace& F = ot.faces[index_f];
		if (F.leaf) U[index_f] = Ustar[index_f] - WGp[index_f];
		else		U[index_f] = 0;
	}
	//----------------------------------------------
	// error = U-Uexct
	//----------------------------------------------
	ARRAY<double>  Error(ot.faces.size()); int cnt = 0;
	double sum = 0;
	for (int index_f = 0; index_f < ot.faces.size(); index_f++) {
		const OcFace& F = ot.faces[index_f];
		if (F.leaf) {
			double x = ot.x_fr_i((F.imin + F.imax) / 2.);
			double y = ot.y_fr_j((F.jmin + F.jmax) / 2.);
			double z = ot.z_fr_k((F.kmin + F.kmax) / 2.);
			     if(F.type_ijk()==0) Error[index_f] = u(x,y,z)-U[index_f];
			else if(F.type_ijk()==1) Error[index_f] = v(x,y,z)-U[index_f];
			else                     Error[index_f] = w(x,y,z)-U[index_f];
			sum += Error[index_f] * Error[index_f];
			cnt++;
		}
		else
			Error[index_f] = 0;
	}
	printf("error in L2=%f\n", sqrt(sum / cnt));
	return 0;
}
