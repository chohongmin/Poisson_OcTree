#include "OcTree.h"
#include <math.h>

double phi (double x,double y,double z){return 1-sqrt(x*x+y*y+z*z);}
double u   (double x,double y,double z){return x*x*z+3*y*y*z;}
double v   (double x,double y,double z){return  -2*x*y*z;}
double w   (double x,double y,double z){return  -x*x*x-x*y*y;}
double ustar(double x,double y,double z){return u(x,y,z)+exp(x-y+z);}
double vstar(double x,double y,double z){return v(x,y,z)-exp(x-y+z);}
double wstar(double x,double y,double z){return w(x,y,z)+exp(x-y+z);}

int main() {
	//----------------------------------------------
	// 기본 octree 생성
	//----------------------------------------------
	OcTree ot(1024);
	ot.xmin = -1.6, ot.xmax = 1.6;
	ot.ymin = -1.6, ot.ymax = 1.6;
	ot.zmin = -1.6, ot.zmax = 1.6;
	ot.generate_ot(phi,1,32);
	//ot.split_all_cells();
	//ot.split_all_cells();

	//----------------------------------------------
	// Ustar
	//----------------------------------------------
	ARRAY<double> Ustar(ot.num_faces());
	for(int index_F=0;index_F<ot.num_faces();index_F++){
		const OcFace& F=ot.faces[index_F];
		double x = ot.x_fr_i((F.imin+F.imax)/2.);
		double y = ot.y_fr_j((F.jmin+F.jmax)/2.);
		double z = ot.z_fr_k((F.kmin+F.kmax)/2.);
		     if(F.type_ijk()==0) Ustar[index_F]=ustar(x,y,z);
		else if(F.type_ijk()==1) Ustar[index_F]=vstar(x,y,z);
		else                     Ustar[index_F]=wstar(x,y,z);
	}
	//----------------------------------------------
	// visualization
	//----------------------------------------------
	ot.visualize_zsection(phi,Ustar,0.1);

	//----------------------------------------------
	// Heaviside function
	//----------------------------------------------
	ARRAY<double> H(ot.num_faces());
	ot.calculate_H(phi,H);

	//----------------------------------------------
	// adjoint test
	//----------------------------------------------
	ARRAY<double> f1(ot.cells.size());
	ARRAY<double> f2(ot.cells.size());
	for (int i = 0; i < ot.cells.size(); i++) {
		f1[i] = ((double)rand()) / RAND_MAX;
		f2[i] = ((double)rand()) / RAND_MAX;
	}
	ARRAY<double>    Gf1(ot.faces.size()); ot.calculate_G(f1,Gf1);
	ARRAY<double>    Gf2(ot.faces.size()); ot.calculate_G(f2,Gf2);
	ARRAY<double>   HGf1(ot.faces.size()); ot.calculate_H(Gf1,H,HGf1);
	ARRAY<double>   HGf2(ot.faces.size()); ot.calculate_H(Gf2,H,HGf2);
	ARRAY<double>  WHGf1(ot.faces.size()); ot.calculate_W(HGf1, WHGf1);
	ARRAY<double>  WHGf2(ot.faces.size()); ot.calculate_W(HGf2, WHGf2);
	ARRAY<double> DWHGf1(ot.faces.size()); ot.calculate_D(WHGf1,DWHGf1);
	ARRAY<double> DWHGf2(ot.faces.size()); ot.calculate_D(WHGf2,DWHGf2);

	double f1_DWHGf2 = ot.inner_product(f1,DWHGf2);
	double DWHGf1_f2 = ot.inner_product(DWHGf1,f2);


	//----------------------------------------------
	// calculation of U* and DWHU*
	//----------------------------------------------
	ARRAY<double>   Hustar(ot.faces.size());
	ARRAY<double>  WHustar(ot.faces.size());
	ARRAY<double> DWHustar(ot.cells.size());

	ot.calculate_H(Ustar,H,Hustar);
	ot.calculate_W(Hustar,WHustar);
	ot.calculate_D(WHustar,DWHustar);

	//----------------------------------------------
	// calculation of p through linear system DbetaGp=DU*
	//----------------------------------------------
	ARRAY<double> p(ot.cells.size());
	ot.solve_Laplacian(H,DWHustar,p);

	//----------------------------------------------
	// calculation of U = U*-WGp
	//----------------------------------------------
	ARRAY<double>  Gp(ot.faces.size());
	ARRAY<double> WGp(ot.faces.size());
	ARRAY<double>  U (ot.faces.size());
	ot.calculate_G( p,  Gp);
	ot.calculate_W(Gp,WGp);

	for (int index_f = 0; index_f < ot.faces.size(); index_f++) {
		const OcFace& F = ot.faces[index_f];
		if(F.leaf && H[index_f]>0){
			U[index_f] = Ustar[index_f] - WGp[index_f];
		}
		else
			U[index_f]=0;
	}
	//----------------------------------------------
	// visualization
	//----------------------------------------------
	ARRAY<double>  WU (ot.faces.size());
	ot.calculate_W(U,WU);
	ot.visualize_zsection(phi,WU,1);

	//----------------------------------------------
	// error = U-Uexct
	//----------------------------------------------
	ARRAY<double>  Error(ot.faces.size()); 
	double sum=0;
	for (int index_f = 0; index_f < ot.faces.size(); index_f++) {
		const OcFace& F = ot.faces[index_f];
		if(F.leaf && H[index_f]>0){
			double x = ot.x_fr_i((F.imin + F.imax) / 2.);
			double y = ot.y_fr_j((F.jmin + F.jmax) / 2.);
			double z = ot.z_fr_k((F.kmin + F.kmax) / 2.);

			     if(F.type_ijk()==0) Error[index_f] = u(x,y,z) - U[index_f];
			else if(F.type_ijk()==1) Error[index_f] = v(x,y,z) - U[index_f];
			else                     Error[index_f] = w(x,y,z) - U[index_f];  
			
			sum += Error[index_f]*Error[index_f]*H[index_f]
				   *ot.face_Area(index_f)*ot.face_Delta(index_f);
		}
		else
			Error[index_f]=0;
	}
	printf("error in L2=%f\n", sqrt(sum));
	return 0;
}
