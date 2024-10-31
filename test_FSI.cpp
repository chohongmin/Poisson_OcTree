#include "OcTree.h"
#include <math.h>

#define PI 3.14159265358979323846264338328

double phi(double x, double y, double z){ return sqrt(x*x+y*y+(z-0.1)*(z-0.1))-1;}

int main(){
	//----------------------------------------------
	// generation of octree
	//----------------------------------------------
	OcTree ot(2048);
	ot.xmin=-2; ot.xmax=2;
	ot.ymin=-2; ot.ymax=2;
	ot.zmin=-2; ot.zmax=2;
	ot.generate_ot(phi,1,128);
	ot.split_all_cells();
	ot.split_all_cells();
	ot.split_all_cells();

	//----------------------------------------------
	// initilization of U*,v*,w*
	//----------------------------------------------
	ARRAY<double> Ustar(ot.num_faces());
	for(int index_F=0;index_F<ot.num_faces();index_F++){
		const OcFace& F = ot.faces[index_F];
		if(F.type_ijk()==2)
			Ustar[index_F]=(F.kmin==0 || F.kmin==ot.N)?0:-1;
		else
			Ustar[index_F]=0;
	}
	double vstarx= 0, wstarx=0, cx=0;
	double vstary= 0, wstary=0, cy=0;
	double vstarz=-1, wstarz=0, cz=0;

	ot.visualize_ysection(phi,Ustar,0.1);

	//----------------------------------------------
	// calculation of H, GH and J 
	//----------------------------------------------
	ARRAY<double>  H(ot.num_faces()); ot.calculate_H(phi,H);
	ARRAY<double> WH(ot.num_faces()); ot.calculate_W(H,WH);
	ARRAY<double>  GHx(ot.num_cells()), GHy(ot.num_cells()), GHz(ot.num_cells());
	ARRAY<double>   Jx(ot.num_cells()),  Jy(ot.num_cells()),  Jz(ot.num_cells());
	ot.calculate_GH_J(WH,GHx,GHy,GHz,Jx,Jy,Jz,cx,cy,cz);

	//----------------------------------------------
	// calculation of DWHU - vstar*GH - wstar*J
	//----------------------------------------------
	ARRAY<double>    HUstar(ot.num_faces());
	ARRAY<double>   WHUstar(ot.num_faces());
	ARRAY<double>  DWHUstar(ot.num_cells());
	ot.calculate_H(  Ustar,H,HUstar);
	ot.calculate_W( HUstar, WHUstar);
	ot.calculate_D(WHUstar,DWHUstar);

	ARRAY<double> RHS(ot.num_cells());
	for(int index_C=0;index_C<ot.num_cells();index_C++){
		const OcCell& C = ot.cells[index_C];
		if(C.leaf)
			RHS[index_C] = DWHUstar[index_C] - vstarx*GHx[index_C] - wstarx*Jx[index_C]
			- vstary*GHy[index_C] - wstary*Jy[index_C]
			- vstarz*GHz[index_C] - wstarz*Jz[index_C];
		else
			RHS[index_C]=0;
	}
	//----------------------------------------------
	// solve DWHGp - GH*minv*GHp - J*Iinv*Jp = RHS
	//----------------------------------------------
	ARRAY<double> p(ot.num_cells());
	double minv = 3./(8*PI);
	double Iinv[3][3]={{15./(16*PI),0,0},
		{0,15./(16*PI),0},
		{0,0,15./(16*PI)}};
	ot.solve_FSI(H,RHS,minv,Iinv,cx,cy,cz,p);

	//----------------------------------------------
	// calculation of U = U*-WGp
	//----------------------------------------------
	ARRAY<double>  Gp(ot.num_faces());
	ARRAY<double> WGp(ot.num_faces());
	ARRAY<double>  U (ot.num_faces());
	ARRAY<double> WU (ot.num_faces());
	ot.calculate_G( p,  Gp);
	ot.calculate_W(Gp, WGp);

	for (int index_F = 0; index_F < ot.num_faces(); index_F++){
		const OcFace& F = ot.faces[index_F];
		if(F.leaf && H[index_F]>0){
			U[index_F] = Ustar[index_F] - WGp[index_F];
		}
		else U[index_F]=0;
	}
	ot.calculate_W(U,WU);
	ot.visualize_ysection(phi,WU,0.1);
	//----------------------------------------------
	// calculation of v = <p,GH>, w=<p,J>
	//----------------------------------------------
	// v = minv * <GH,p>
	double vx = minv * ot.inner_product(p,GHx);
	double vy = minv * ot.inner_product(p,GHy);
	double vz = minv * ot.inner_product(p,GHz);
	// w = Iinv * <J,p>
	double tempx = ot.inner_product(p, Jx);
	double tempy = ot.inner_product(p, Jy);
	double tempz = ot.inner_product(p, Jz);
	double wx = Iinv[0][0]*tempx + Iinv[0][1]*tempy + Iinv[0][2]*tempz;
	double wy = Iinv[1][0]*tempx + Iinv[1][1]*tempy + Iinv[1][2]*tempz;
	double wz = Iinv[2][0]*tempx + Iinv[2][1]*tempy + Iinv[2][2]*tempz;

	//----------------------------------------------
	// save v,w,WHU
	//----------------------------------------------
	FILE* fp = fopen("error.txt","wt");
	fprintf(fp,"%f,%f,%f\n",vx,vy,vz);
	fprintf(fp,"%f,%f,%f\n",wx,wy,wz);
	for (int index_f = 0; index_f < ot.num_faces(); index_f++) {
		if(index_f==6749)
			bool stop=true;
		const OcFace& F = ot.faces[index_f];
		double face_area  = ot.face_Area (index_f);
		double face_delta = ot.face_Delta(index_f);

		if(F.leaf && H[index_f]>0 && face_delta>0 ){
			fprintf(fp,"l,%d,%f,%f\n",index_f, WU[index_f], face_area*face_delta*H[index_f]);
		}
		if(F.leaf==false && ot.faces[F.index_child[0][0]].leaf==true
			&& ot.faces[F.index_child[0][1]].leaf==true
			&& ot.faces[F.index_child[1][0]].leaf==true
			&& ot.faces[F.index_child[1][1]].leaf==true  && face_delta>0 ){
			double sum_wu=0; double sum_h=0;
			for(int i=0;i<2;i++)
				for(int j=0;j<2;j++){
					double wu = WU[F.index_child[i][j]];
					double  h =  H[F.index_child[i][j]];
					if(h>0){ sum_wu+=wu*h; sum_h+=h;}
				}
			double wu = sum_wu/sum_h;
			if (sum_h >0) fprintf(fp,"p,%d,%f\n",index_f, wu);
		}
	}
	return 0;
}
