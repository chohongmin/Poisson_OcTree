#include "./library/geometry/MarchingCubes.h"
#include "OcTree.h"

void OcCell::set_child(int i, int j, int k, OcCell& child) const{
	child.imin = (i == 0) ? imin : (imin + imax) / 2;
	child.jmin = (j == 0) ? jmin : (jmin + jmax) / 2;
	child.kmin = (k == 0) ? kmin : (kmin + kmax) / 2;
	child.imax = (i == 1) ? imax : (imin + imax) / 2;
	child.jmax = (j == 1) ? jmax : (jmin + jmax) / 2;
	child.kmax = (k == 1) ? kmax : (kmin + kmax) / 2;
}

OcCell::OcCell() {
	index_prnt = -1;
	for(int i=0;i<2;i++)
	for(int j=0;j<2;j++)
	for(int k=0;k<2;k++){
		index_child[i][j][k]=-1;
	}
	for(int i=0;i<3;i++)
	for(int j=0;j<2;j++){
		index_face[i][j] = -1;
	}
	leaf = true;
}

OcFace::OcFace() {
	index_prnt = -1;
	for (int i = 0; i < 2; i++){
		for (int j=0; j<2; j++){
		index_child[i][j] = -1;
		}
		index_cell[i]  = -1;
	}
	leaf = true;
}

int OcFace::type_ijk() const{
	if(imin==imax) return 0;
	if(jmin==jmax) return 1;
	if(kmin==kmax) return 2;
	assert (false);
	return -1;
}


void OcFace::set_child(int a, int b, OcFace& child) const{
	int type = type_ijk();
	
	if (type==0){
		child.imin = imin; child.imax = imax;
		child.jmin = (a == 0) ? jmin : (jmin + jmax) / 2; 
		child.jmax = (a == 1) ? jmax : (jmin + jmax) / 2; 
		child.kmin = (b == 0) ? kmin : (kmin + kmax) / 2;
		child.kmax = (b == 1) ? kmax : (kmin + kmax) / 2;
	}
	if (type==1) {
		child.jmin = jmin; child.jmax = jmax;
		child.imin = (a == 0) ? imin : (imin + imax) / 2;
		child.imax = (a == 1) ? imax : (imin + imax) / 2;
		child.kmin = (b == 0) ? kmin : (kmin + kmax) / 2;
		child.kmax = (b == 1) ? kmax : (kmin + kmax) / 2;
	}
	if (type==2) {
		child.kmin = kmin; child.kmax = kmax;
		child.imin = (a == 0) ? imin : (imin + imax) / 2;
		child.imax = (a == 1) ? imax : (imin + imax) / 2;
		child.jmin = (b == 0) ? jmin : (jmin + jmax) / 2;
		child.jmax = (b == 1) ? jmax : (jmin + jmax) / 2;
	}
}

double OcTree::face_Area(int index_F)const{
	const OcFace& F = faces[index_F];
	double xmin = x_fr_i(F.imin), xmax=x_fr_i(F.imax);
	double ymin = y_fr_j(F.jmin), ymax=y_fr_j(F.jmax);
	double zmin = z_fr_k(F.kmin), zmax=z_fr_k(F.kmax);
	if(F.type_ijk()==0) return (ymax-ymin)*(zmax-zmin);
	if(F.type_ijk()==1) return (xmax-xmin)*(zmax-zmin);
	if(F.type_ijk()==2) return (ymax-ymin)*(xmax-xmin);
	assert(false); return 0;
}

double OcTree::face_Delta(int index_F)const{
	const OcFace& F = faces[index_F];
	int c0 = F.index_cell[0];
	int c1 = F.index_cell[1];
	if( c0==-1 || c1==-1 ) return 0;
	if(F.type_ijk()==0) return .5*(xsize(c0)+xsize(c1));
	if(F.type_ijk()==1) return .5*(ysize(c0)+ysize(c1));
	if(F.type_ijk()==2) return .5*(zsize(c0)+zsize(c1));
	assert(false); return 0;
}

OcTree::OcTree( int N ){
	this->N = N;
	OcCell C;
	C.imin = 0; C.imax = N;
	C.jmin = 0; C.jmax = N;
	C.kmin = 0; C.kmax = N;

	for(int type=0; type<3; type++)
	for(int l=0; l<2; l++){
		OcFace F; C.set_face(type,l,F);
		F.index_cell[1-l]=0;
		C.index_face[type][l] = num_faces();
		faces.push(F);
	}
	cells.push(C);
}

void OcCell::set_face(int type, int l, OcFace &F) const{
	if(type==0){
		F.jmin = jmin; F.jmax = jmax;
		F.kmin = kmin; F.kmax = kmax;
		F.imin = (l==0)? imin : imax;
		F.imax = (l==0)? imin : imax;
	}
	if(type==1){
		F.imin = imin; F.imax = imax;
		F.kmin = kmin; F.kmax = kmax;
		F.jmin = (l==0)? jmin : jmax;
		F.jmax = (l==0)? jmin : jmax;
	}
	if(type==2){
		F.imin = imin; F.imax = imax;
		F.jmin = jmin; F.jmax = jmax;
		F.kmin = (l==0)? kmin : kmax;
		F.kmax = (l==0)? kmin : kmax;
	}
}

void OcTree::split_f(int index_f){
	OcFace F = faces[index_f]; 
	if(F.leaf==false) return;
	F.leaf = false;
	for(int a=0; a<2; a++)
	for(int b=0; b<2; b++){
		OcFace child;
		F.set_child(a,b,child); child.index_prnt = index_f;
		F.index_child[a][b] = num_faces();
		faces.push(child);
	}
	faces[index_f] = F;
}

int OcTree::search_leafcell(double i, double j, double k) const{
	int index_c = 0;
	if(i<0 || j<0 || k<0|| i>N || j>N || k>N) return -1;
	while(cells[index_c].leaf == false){
		const OcCell& C = cells[index_c];
		int icen = (C.imin + C.imax) /2;
		int jcen = (C.jmin + C.jmax) /2;
		int kcen = (C.kmin + C.kmax) /2;
		int ioff = (i<icen)? 0:1;
		int joff = (j<jcen)? 0:1;
		int koff = (k<kcen)? 0:1;
		index_c = C.index_child[ioff][joff][koff];
	}
	return index_c;
}

void OcTree::set_face_cell_relation(int index_f){
	OcFace& F = faces[index_f];
	
	int icen = (F.imin + F.imax) /2;
	int jcen = (F.jmin + F.jmax) /2;
	int kcen = (F.kmin + F.kmax) /2;

	double ioff = (F.imin==F.imax)?1e-5:0;
	double joff = (F.jmin==F.jmax)?1e-5:0;
	double koff = (F.kmin==F.kmax)?1e-5:0;

	for(int l=0; l<2; l++){
		int index_c=search_leafcell(icen+ioff*(l-0.5),jcen+joff*(l-0.5),kcen+koff*(l-0.5));
		if(index_c!=-1){
			OcCell&C = cells[index_c];
			F.index_cell[l]=index_c;
		
			if(F.size()==C.size()){
				int type_ijk=F.type_ijk();
				C.index_face[type_ijk][1-l]=index_f;
			}
		}
	}
}

void OcTree::split_c(int index_c){
	// split cells
	OcCell C=cells[index_c]; C.leaf = false;
	for(int i=0; i<2; i++)
	for(int j=0; j<2; j++)
	for(int k=0; k<2; k++){
		OcCell child;
		C.set_child(i,j,k,child); child.index_prnt = index_c;
		C.index_child[i][j][k] = num_cells();
		cells.push(child);
	}
	cells[index_c]=C;

	// split boundary faces
	for(int type=0; type<3; type++)
	for(int l=0; l<2; l++){
		int index_f = C.index_face[type][l];
		split_f(index_f);
		for(int a=0; a<2; a++)
		for(int b=0; b<2; b++){
			OcFace F = faces[index_f];
			set_face_cell_relation(F.index_child[a][b]);
		}
	}

	// split inside faces
	for(int type=0; type<3; type++){
		OcFace F;
		F.imin=C.imin; F.imax=C.imax; 
		F.jmin=C.jmin; F.jmax=C.jmax; 
		F.kmin=C.kmin; F.kmax=C.kmax;

		if(type==0){F.imin = F.imax = (C.imin + C.imax) /2;	}
		if(type==1){F.jmin = F.jmax = (C.jmin + C.jmax) /2;	}
		if(type==2){F.kmin = F.kmax = (C.kmin + C.kmax) /2;	}

		for(int a=0; a<2; a++)
		for(int b=0; b<2; b++){
			int index_f = num_faces();
			OcFace child;
			F.set_child(a,b,child);
			faces.push(child);
			set_face_cell_relation(index_f);
		}
	}
}

double OcTree:: calculate_G (const ARRAY <double> & u, int index_f) const{
	const OcFace &F = faces[index_f];
	if (F.leaf){ int type_ijk = F.type_ijk();
		if (F.index_cell[0]==-1||F.index_cell[1]==-1)
			return 0;
		else if (type_ijk==0){
			int c0 = F.index_cell[0];
			int c1 = F.index_cell[1];
				return(u[c1]-u[c0])/(xsize(c1)*.5+xsize(c0)*.5);
		}
		else if (type_ijk == 1) {
			int c0 = F.index_cell[0];
			int c1 = F.index_cell[1];
			return(u[c1] - u[c0]) / (ysize(c1) * .5 + ysize(c0) * .5);
		}
		else if (type_ijk == 2) {
			int c0 = F.index_cell[0];
			int c1 = F.index_cell[1];
			return(u[c1] - u[c0]) / (zsize(c1) * .5 + zsize(c0) * .5);
		}
		else{
		assert(false);
		return 0;
		}
	}
}
void OcTree::split_all_cells() {
	int size_c = cells.size();
	for (int index_c = 0; index_c < size_c; index_c++) {
		if (cells[index_c].leaf)
			split_c(index_c);
	}
}

double OcTree::calculate_flux(const ARRAY<double>& f, int index_f) const {
	const OcFace& F = faces[index_f];
	int type_ijk = F.type_ijk();
	if(F.leaf){
		      if(type_ijk==0) return f[index_f]*(y_fr_j(F.jmax)-y_fr_j(F.jmin))*(z_fr_k(F.kmax)-z_fr_k(F.kmin));
		 else if(type_ijk==1) return f[index_f]*(x_fr_i(F.imax)-x_fr_i(F.imin))*(z_fr_k(F.kmax)-z_fr_k(F.kmin));
		 else                 return f[index_f]*(x_fr_i(F.imax)-x_fr_i(F.imin))*(y_fr_j(F.jmax)-y_fr_j(F.jmin));
	}
	else{
		return calculate_flux(f,F.index_child[0][0])
			+calculate_flux(f, F.index_child[0][1])
			+calculate_flux(f, F.index_child[1][0])
			+calculate_flux(f, F.index_child[1][1]);
	}
}

void OcTree::calculate_D(const ARRAY<double>& f, ARRAY<double>& D) const {
	for (int index_d=0; index_d<cells.size(); index_d++){
		const OcCell& C=cells[index_d];
		if(C.leaf==true){
			D[index_d]=calculate_flux(f, C.index_face[0][1])
					  -calculate_flux(f, C.index_face[0][0])
					  +calculate_flux(f, C.index_face[1][1])
					  -calculate_flux(f, C.index_face[1][0])
					  +calculate_flux(f, C.index_face[2][1])
					  -calculate_flux(f, C.index_face[2][0]);
		}
		else D[index_d]=0;
	}
}


void OcTree::calculate_W(const ARRAY<double>& f, int index_F, double& Wf, double& delta) const {
	const OcFace & F=faces[index_F];
	if(F.leaf) {
		Wf=f[index_F];
		int c0=F.index_cell[0];
		int c1=F.index_cell[1];
		if(c0==-1||c1==-1) {delta=1; return;}
			 if(F.type_ijk()==0) delta=.5*(xsize(c0)+xsize(c1));
		else if(F.type_ijk()==1) delta=.5*(ysize(c0)+ysize(c1));
		   else delta=.5*(zsize(c0)+zsize(c1));
		return;
	}
	double Wf0, delta0, Wf1, delta1, Wf2, delta2, Wf3, delta3;
	calculate_W(f, F.index_child[0][0],Wf0, delta0);
	calculate_W(f, F.index_child[0][1], Wf1, delta1);
	calculate_W(f, F.index_child[1][0], Wf2, delta2);
	calculate_W(f, F.index_child[1][1], Wf3, delta3);

	Wf=(Wf0*delta0+Wf1*delta1+Wf2*delta2+Wf3*delta3)/(delta0+delta1+delta2+delta3);
	delta=(delta0 + delta1 + delta2 + delta3)/4;
}


void OcTree::set_W(int index_F, double Wf_value, ARRAY<double>& Wf) const{
	const OcFace &F=faces[index_F];
	Wf[index_F]=Wf_value;
	if(F.leaf==false) {
		set_W(F.index_child[0][0], Wf_value, Wf);
		set_W(F.index_child[0][1], Wf_value, Wf);
		set_W(F.index_child[1][0], Wf_value, Wf);
		set_W(F.index_child[1][1], Wf_value, Wf);
	}
}


void OcTree::calculate_G(const ARRAY<double>& u, ARRAY<double>& Gu) const {
	for (int i = 0; i < faces.size(); i++) {
		const OcFace& F = faces[i];
		if (F.leaf)
			Gu[i] = calculate_G(u, i);
		else
			Gu[i] = 0;
	}
}


void OcTree::calculate_W(const ARRAY<double>& f, ARRAY<double>& Wf) const {
	for (int index_F = 0; index_F < faces.size(); index_F++) {
		const OcFace& F = faces[index_F];
		int c0 = F.index_cell[0];
		int c1 = F.index_cell[1];
		if (c0 == -1 || c1 == -1) Wf[index_F] = f[index_F];
		else {
			const OcCell& C0 = cells[c0]; int S0 = C0.imax - C0.imin;
			const OcCell& C1 = cells[c1]; int S1 = C1.imax - C1.imin;
			if (C0.leaf || C1.leaf) {
				if (C0.leaf && C1.leaf) {
					if (S0 == S1) Wf[index_F] = f[index_F];
					else;
				}
				else {
					if (S0 == S1) {
						double delta;
						calculate_W(f, index_F, Wf[index_F], delta);
						set_W(index_F, Wf[index_F], Wf);
					}
					else;
				}
			}
			else Wf[index_F] = f[index_F];
		}
	}
}

double OcTree::inner_product(const ARRAY<double>& x, const ARRAY<double>& y) const {
	double sum = 0;
	for (int index_c = 0; index_c < cells.size(); index_c++)
		if (cells[index_c].leaf)
			sum += x[index_c] * y[index_c];
	return sum;
}

void OcTree::projection_sumzero(ARRAY<double>& x) const {
	double sum = 0; int cnt = 0;
	for (int index_c = 0; index_c < cells.size(); index_c++)
		if (cells[index_c].leaf) {
			sum += x[index_c];
			cnt++;
		}
	double avg = sum / cnt;
	for (int index_c = 0; index_c < cells.size(); index_c++)
		if (cells[index_c].leaf)
			x[index_c] -= avg;
}

void OcTree::solve_Laplacian(const ARRAY<double>& f, ARRAY<double>& u) const {
	ARRAY<double>& x = u;
	ARRAY<double>  r(cells.size());
	ARRAY<double>  p(cells.size());
	ARRAY<double> Ap(cells.size());
	for (int i = 0; i < cells.size(); i++) {
		x[i] = 0; r[i] = f[i]; p[i] = f[i];
	}

	projection_sumzero(r);
	projection_sumzero(p);
	double rr = inner_product(r, r);
	double rr0 = rr; int cnt = 0;
	ARRAY<double>  Gp(faces.size());
	ARRAY<double> WGp(faces.size());

	while (rr > rr0 * 1e-15) {
		// Ap=DWGp
		calculate_G(p, Gp);
		calculate_W(Gp, WGp);
		calculate_D(WGp, Ap);

		projection_sumzero(Ap);
		double pAp = inner_product(p, Ap);
		double alpha = rr / pAp;

		for (int i = 0; i < cells.size(); i++) {
			x[i] += alpha * p[i];
			r[i] -= alpha * Ap[i];
		}

		projection_sumzero(r);
		double rr_next = inner_product(r, r);
		double beta_ = rr_next / rr;
		rr = rr_next;

		for (int i = 0; i < cells.size(); i++)
			p[i] = r[i] + beta_ * p[i];

		if (cnt % 10 == 0)
			printf("cnt=%d, rr=%f\n", cnt, rr);
		cnt++;
	}
}

void OcTree::visualize_domain( double (*phi)(double,double,double), Display& disp){
	disp.set_physical_coordinates(xmin,xmax,ymin,ymax,zmin,zmax);
	disp.make_it_white();
	for(int index_c=0;index_c<num_cells();index_c++){
		const OcCell& C = cells[index_c];
		double xmin=x_fr_i(C.imin), xmax=x_fr_i(C.imax);
		double ymin=y_fr_j(C.jmin), ymax=y_fr_j(C.jmax);
		double zmin=z_fr_k(C.kmin), zmax=z_fr_k(C.kmax);
		double phiC[2][2][2];
		for(int i=0;i<2;i++){ double x=(i==0)?xmin:xmax;
		for(int j=0;j<2;j++){ double y=(j==0)?ymin:ymax;
		for(int k=0;k<2;k++){ double z=(k==0)?zmin:zmax;
			phiC[i][j][k] = phi(x,y,z);
		}}}
		STACK<POINT3> pts;
		MarchingCubes::isosurface(phiC,xmin,xmax,ymin,ymax,zmin,zmax,pts);
		for(int n=0;n<pts.size();n+=3){
			POINT3 N; POINT3::cross(pts[n+1]-pts[n],
			                        pts[n+2]-pts[n], N);
			N /= N.abs();
			int rgb_shading = disp.Phong_Shading(N, 0xff0000);
			disp.draw_triangle(pts[n],pts[n+1],pts[n+2], rgb_shading);
		}
	}
}

void OcTree::visualize_xsection( Display& disp ){
	for(int index_f=0; index_f<num_faces();index_f++){
		const OcFace& F = faces[index_f];
		if(F.type_ijk()==0 && F.imin == N/2){
			double x = x_fr_i(F.imin);
			double ymin = y_fr_j(F.jmin), ymax=y_fr_j(F.jmax);
			double zmin = z_fr_k(F.kmin), zmax=z_fr_k(F.kmax);
			disp.draw_line(POINT3(x,ymin,zmin),POINT3(x,ymin,zmax),0x000000);
			disp.draw_line(POINT3(x,ymin,zmin),POINT3(x,ymax,zmin),0x000000);
			disp.draw_line(POINT3(x,ymax,zmax),POINT3(x,ymin,zmax),0x000000);
			disp.draw_line(POINT3(x,ymax,zmax),POINT3(x,ymax,zmin),0x000000);
		}
	}
}

void OcTree::visualize_ysection( Display& disp ){
	for(int index_f=0; index_f<num_faces();index_f++){
		const OcFace& F = faces[index_f];
		if(F.type_ijk()==1 && F.jmin == N/2){
			double y = y_fr_j(F.jmin);
			double xmin = x_fr_i(F.imin), xmax=x_fr_i(F.imax);
			double zmin = z_fr_k(F.kmin), zmax=z_fr_k(F.kmax);
			disp.draw_line(POINT3(xmin,y,zmin),POINT3(xmin,y,zmax),0x000000);
			disp.draw_line(POINT3(xmin,y,zmin),POINT3(xmax,y,zmin),0x000000);
			disp.draw_line(POINT3(xmax,y,zmax),POINT3(xmin,y,zmax),0x000000);
			disp.draw_line(POINT3(xmax,y,zmax),POINT3(xmax,y,zmin),0x000000);
		}
	}
}

void OcTree::visualize_zsection( Display& disp ){
	for(int index_f=0; index_f<num_faces();index_f++){
		const OcFace& F = faces[index_f];
		if(F.type_ijk()==2 && F.kmin == N/2){
			double z = z_fr_k(F.kmin);
			double ymin = y_fr_j(F.jmin), ymax=y_fr_j(F.jmax);
			double xmin = x_fr_i(F.imin), xmax=x_fr_i(F.imax);
			disp.draw_line(POINT3(xmin,ymin,z),POINT3(xmax,ymin,z),0x000000);
			disp.draw_line(POINT3(xmin,ymin,z),POINT3(xmin,ymax,z),0x000000);
			disp.draw_line(POINT3(xmax,ymax,z),POINT3(xmax,ymin,z),0x000000);
			disp.draw_line(POINT3(xmax,ymax,z),POINT3(xmin,ymax,z),0x000000);
		}
	}
}

void OcTree::generate_ot(double(*phi)(double, double,double), double Lip, int Csize_min) {
	bool change_on_ot = true;
	while (change_on_ot) {
		change_on_ot = false;
		for (int i = cells.size() - 1; i >= 0; i--) {
			const OcCell& C = cells[i];
			double diag = sqrt(xsize(i)*xsize(i)+ysize(i)*ysize(i)+zsize(i)*zsize(i));
			if (C.leaf && (C.imax - C.imin) > Csize_min) {
				double phiv_abs_min = 1e8;
				for(int i=0;i<2;i++){ double x=x_fr_i((i==0)?C.imin:C.imax);
				for(int j=0;j<2;j++){ double y=y_fr_j((j==0)?C.jmin:C.jmax);
				for(int k=0;k<2;k++){ double z=z_fr_k((k==0)?C.kmin:C.kmax);
					double phiv=phi(x,y,z);
					if (fabs(phiv) < phiv_abs_min)
						phiv_abs_min = fabs(phiv);
				}}}
				if (phiv_abs_min <= Lip * diag / 2) {
					split_c(i); change_on_ot = true;
				}
			}
		}
	}
}

POINT3 interpolate_phizero( double phi0, double phi1, POINT3 P0, POINT3 P1 ){
	assert(phi0*phi1<0);
	return (P0*phi1-P1*phi0)/(phi1-phi0);
}


void OcTree::visualize_ysection(double (*phi)(double, double, double), const ARRAY<double>& U, double scale) {
	Display disp;
	disp.set_physical_coordinates(xmin, xmax, ymin, ymax, zmin, zmax);
	disp.make_it_white();
	disp.m_rotation.set_rotation(-PI/2,1,0,0);
	for(int index_F=0;index_F<num_faces();index_F++){
		const OcFace& F = faces[index_F];
		if(F.jmin==N/2 && F.jmax==N/2 && F.leaf){
			//
			double y    = y_fr_j(N / 2);
			double xmin = x_fr_i(F.imin), xmax = x_fr_i(F.imax);
			double zmin = z_fr_k(F.kmin), zmax = z_fr_k(F.kmax);
			double phiF[2][2]; POINT3 vertices[2][2];
			for(int i=0;i<2;i++){double x=(i==0)?xmin:xmax;
			for(int k=0;k<2;k++){double z=(k==0)?zmin:zmax;
			phiF    [i][k] = phi   (x,y,z);
			vertices[i][k] = POINT3(x,y,z);
			}}
			//
			int index_C=-1;
			     if(F.index_cell[0]==-1) index_C=F.index_cell[1];
			else if(F.index_cell[1]==-1) index_C=F.index_cell[0];
			else {
				const OcCell& C0=cells[F.index_cell[0]];
				const OcCell& C1=cells[F.index_cell[1]];
				if(C0.size()==F.size() && C0.leaf) index_C = F.index_cell[0];
				if(C1.size()==F.size() && C1.leaf) index_C = F.index_cell[1];
				assert(index_C != -1);
			}
			//
			const OcCell& C = cells[index_C];
			double uiph = U[C.index_face[0][1]];
			double uimh = U[C.index_face[0][0]];
			double wkph = U[C.index_face[2][1]];
			double wkmh = U[C.index_face[2][0]];
			//
			disp.draw_line(POINT3(xmin,y,zmin),POINT3(xmax,y,zmin),0x000000);
			disp.draw_line(POINT3(xmin,y,zmin),POINT3(xmin,y,zmax),0x000000);
			disp.draw_line(POINT3(xmax,y,zmax),POINT3(xmax,y,zmin),0x000000);
			disp.draw_line(POINT3(xmax,y,zmax),POINT3(xmin,y,zmax),0x000000);
			//
			double x = (xmin+xmax)/2, u=(uiph+uimh)/2.;
			double z = (zmin+zmax)/2, w=(wkph+wkmh)/2.;
			disp.draw_line(POINT3(x,y,z),POINT3(x+scale*u,y,z+scale*w),0xff0000);
			//
			STACK<POINT3> pts;
			if(phiF[0][0]*phiF[1][0]<0) pts.push(interpolate_phizero(phiF[0][0],phiF[1][0],vertices[0][0],vertices[1][0]));
			if(phiF[0][0]*phiF[0][1]<0) pts.push(interpolate_phizero(phiF[0][0],phiF[0][1],vertices[0][0],vertices[0][1]));
			if(phiF[1][1]*phiF[1][0]<0) pts.push(interpolate_phizero(phiF[1][1],phiF[1][0],vertices[1][1],vertices[1][0]));
			if(phiF[1][1]*phiF[0][1]<0) pts.push(interpolate_phizero(phiF[1][1],phiF[0][1],vertices[1][1],vertices[0][1]));
			//
			if(pts.size()==2)
				disp.draw_line(pts[0],pts[1],0x0000ff);
		}
	}
	disp.save_bmp();
}

void OcTree::visualize_zsection(double (*phi)(double, double, double), const ARRAY<double>& U, double scale) {
	Display disp;
	disp.set_physical_coordinates(xmin, xmax, ymin, ymax, zmin, zmax);
	disp.make_it_white();
	for(int index_F=0;index_F<num_faces();index_F++){
		const OcFace& F = faces[index_F];
		if(F.kmin==N/2 && F.kmax==N/2 && F.leaf){
			//
			double z    = z_fr_k(N / 2);
			double xmin = x_fr_i(F.imin), xmax = x_fr_i(F.imax);
			double ymin = y_fr_j(F.jmin), ymax = y_fr_j(F.jmax);
			double phiF[2][2]; POINT3 vertices[2][2];
			for(int i=0;i<2;i++){double x=(i==0)?xmin:xmax;
			for(int j=0;j<2;j++){double y=(j==0)?ymin:ymax;
				phiF    [i][j] = phi(x,y,z);
				vertices[i][j] = POINT3(x,y,z);
			}}
			//
			int index_C=-1;
			     if(F.index_cell[0]==-1) index_C=F.index_cell[1];
			else if(F.index_cell[1]==-1) index_C=F.index_cell[0];
			else {
				const OcCell& C0=cells[F.index_cell[0]];
				const OcCell& C1=cells[F.index_cell[1]];
				if(C0.size()==F.size() && C0.leaf) index_C = F.index_cell[0];
				if(C1.size()==F.size() && C1.leaf) index_C = F.index_cell[1];
				assert(index_C != -1);
			}
			//
			const OcCell& C = cells[index_C];
			double uiph = U[C.index_face[0][1]];
			double uimh = U[C.index_face[0][0]];
			double vjph = U[C.index_face[1][1]];
			double vjmh = U[C.index_face[1][0]];
			//
			disp.draw_line(POINT3(xmin,ymin,z),POINT3(xmax,ymin,z),0x000000);
			disp.draw_line(POINT3(xmin,ymin,z),POINT3(xmin,ymax,z),0x000000);
			disp.draw_line(POINT3(xmax,ymax,z),POINT3(xmax,ymin,z),0x000000);
			disp.draw_line(POINT3(xmax,ymax,z),POINT3(xmin,ymax,z),0x000000);
			//
			double x = (xmin+xmax)/2, u=(uiph+uimh)/2.;
			double y = (ymin+ymax)/2, v=(vjph+vjmh)/2.;
			disp.draw_line(POINT3(x,y,z),POINT3(x+scale*u,y+scale*v,z),0xff0000);
			//
			STACK<POINT3> pts;
			if(phiF[0][0]*phiF[1][0]<0) pts.push(interpolate_phizero(phiF[0][0],phiF[1][0],vertices[0][0],vertices[1][0]));
			if(phiF[0][0]*phiF[0][1]<0) pts.push(interpolate_phizero(phiF[0][0],phiF[0][1],vertices[0][0],vertices[0][1]));
			if(phiF[1][1]*phiF[1][0]<0) pts.push(interpolate_phizero(phiF[1][1],phiF[1][0],vertices[1][1],vertices[1][0]));
			if(phiF[1][1]*phiF[0][1]<0) pts.push(interpolate_phizero(phiF[1][1],phiF[0][1],vertices[1][1],vertices[0][1]));
			//
			if(pts.size()==2)
				disp.draw_line(pts[0],pts[1],0x0000ff);
		}
	}
	disp.save_bmp();
}

// Heaviside value on triangle
double calculate_H( double phi0, double phi1, double phi2 ){
	int count_pos=0;
	if(phi0>0) count_pos++;
	if(phi1>0) count_pos++;
	if(phi2>0) count_pos++;
	if(count_pos==0) return 0;
	if(count_pos==3) return 1;
	bool sign_flip=false;
	if(count_pos==1){ sign_flip=true;
		phi0=-phi0; phi1=-phi1; phi2=-phi2;
	}
	double H=0;
	if(phi0<=0) H=1-phi0/(phi0-phi1)*phi0/(phi0-phi2);
	if(phi1<=0) H=1-phi1/(phi1-phi2)*phi1/(phi1-phi0);
	if(phi2<=0) H=1-phi2/(phi2-phi0)*phi2/(phi2-phi1);
	return sign_flip ? (1-H) : H;
}

// Heaviside value on rectangle
double calculate_H( double phi00, double phi01, double phi10, double phi11 ){
	double phi_middle = (phi00+phi01+phi10+phi11)*.25;
	return (calculate_H(phi_middle,phi00,phi10)
		   +calculate_H(phi_middle,phi10,phi11)
		   +calculate_H(phi_middle,phi11,phi01)
		   +calculate_H(phi_middle,phi01,phi00))*.25;
}

// Heavise values on all the faces
void OcTree::calculate_H( double(*phi)(double, double, double), ARRAY<double>& H) const{
	for(int index_F=0;index_F<num_faces();index_F++){
		const OcFace& F = faces[index_F];
		if(F.leaf){
			double xmin=x_fr_i(F.imin), xmax=x_fr_i(F.imax);
			double ymin=y_fr_j(F.jmin), ymax=y_fr_j(F.jmax);
			double zmin=z_fr_k(F.kmin), zmax=z_fr_k(F.kmax);
			if(F.type_ijk()==0) 
				H[index_F]=::calculate_H(phi(xmin,ymin,zmin),
					                     phi(xmin,ymin,zmax),
					                     phi(xmin,ymax,zmin),
					                     phi(xmin,ymax,zmax));
			else if(F.type_ijk()==1) 
				H[index_F]=::calculate_H(phi(xmin,ymin,zmin),
									     phi(xmin,ymin,zmax),
									     phi(xmax,ymin,zmin),
									     phi(xmax,ymin,zmax));
			else if(F.type_ijk()==2) 
				H[index_F]=::calculate_H(phi(xmin,ymin,zmin),
									     phi(xmin,ymax,zmin),
									     phi(xmax,ymin,zmin),
									     phi(xmax,ymax,zmin));
			else assert(false);
		}
		else
			H[index_F]=0;
	}
}

// irregular domain
void OcTree::calculate_H(const ARRAY<double>& f, const ARRAY<double>& H, ARRAY<double>& Hf )const{
	for(int index_F=0;index_F<faces.size();index_F++){
		const OcFace& F = faces[index_F];
		if(F.leaf)
			Hf[index_F] = f[index_F]*H[index_F];
		else
			Hf[index_F]=0;
	}
}

void OcTree::solve_Laplacian(const ARRAY<double>& H, const ARRAY<double>& f, ARRAY<double>& u) const {
	ARRAY<double>& x = u;
	ARRAY<double>  r(cells.size());
	ARRAY<double>  p(cells.size());
	ARRAY<double> Ap(cells.size());
	for (int i = 0; i < cells.size(); i++) {
		x[i] = 0; r[i] = f[i]; p[i] = f[i];
	}

	projection_sumzero(r);
	projection_sumzero(p);
	double rr = inner_product(r, r);
	double rr0 = rr; int cnt = 0;
	ARRAY<double>   Gp(faces.size());
	ARRAY<double>  HGp(faces.size());
	ARRAY<double> WHGp(faces.size());

	while (rr > rr0 * 1e-15) {
		// Ap=DWGp
		calculate_G( p,Gp);
		calculate_H(H,Gp,HGp);
		calculate_W(HGp, WHGp);
		calculate_D(WHGp, Ap);

		projection_sumzero(Ap);
		double pAp = inner_product(p, Ap);
		double alpha = rr / pAp;

		for (int i = 0; i < cells.size(); i++) {
			x[i] += alpha * p[i];
			r[i] -= alpha * Ap[i];
		}

		projection_sumzero(r);
		double rr_next = inner_product(r, r);
		double beta_ = rr_next / rr;
		rr = rr_next;

		for (int i = 0; i < cells.size(); i++)
			p[i] = r[i] + beta_ * p[i];

		if (cnt % 10 == 0)
			printf("cnt=%d, rr=%e\n", cnt, rr);
		cnt++;
	}
}

	// FSI
void OcTree::calculate_GH_J(const ARRAY<double>& WH, ARRAY<double>& GHx, ARRAY<double>& GHy, ARRAY<double>& GHz, 
	                                                 ARRAY<double>&  Jx, ARRAY<double>&  Jy, ARRAY<double>&  Jz, 
	                       const double cx, const double cy, const double cz) const{
	for(int index_C=0;index_C<num_cells();index_C++){
		const OcCell& C = cells[index_C];
		if(C.leaf){
			double HF[3][2];
			for(int type_ijk=0;type_ijk<3;type_ijk++)
			for(int        n=0;       n<2;       n++)
				HF[type_ijk][n]=WH[C.index_face[type_ijk][n]];
			double volume = xsize(index_C)*ysize(index_C)*zsize(index_C);
			GHx[index_C]=(HF[0][1]-HF[0][0])/xsize(index_C)*volume;
			GHy[index_C]=(HF[1][1]-HF[1][0])/ysize(index_C)*volume;
			GHz[index_C]=(HF[2][1]-HF[2][0])/zsize(index_C)*volume;
			POINT3 J = POINT3::cross(POINT3(x_fr_i((C.imin+C.imax)/2.)-cx,
				                            y_fr_j((C.jmin+C.jmax)/2.)-cy,
				                            z_fr_k((C.kmin+C.kmax)/2.)-cz),
				                     POINT3(GHx[index_C],GHy[index_C],GHz[index_C]));
			Jx[index_C]=J.x;
			Jy[index_C]=J.y;
			Jz[index_C]=J.z;
		}
		else{
			GHx[index_C]=0;
			GHy[index_C]=0;
			GHz[index_C]=0;
			 Jx[index_C]=0;
			 Jy[index_C]=0;
			 Jz[index_C]=0;
		}
	}
}
void OcTree::solve_FSI    ( const ARRAY<double>& H, const ARRAY<double>& RHS, 
		           double minv, double Iinv[3][3], double cx, double cy, double cz, ARRAY<double>& pressure) const{
	ARRAY<double> GHx(num_cells());
	ARRAY<double> GHy(num_cells());
	ARRAY<double> GHz(num_cells());
	ARRAY<double>  Jx(num_cells());
	ARRAY<double>  Jy(num_cells());
	ARRAY<double>  Jz(num_cells()); 

	// calculate WH, GH and J
	ARRAY<double> WH(num_faces());
	calculate_W(H,WH);	
	calculate_GH_J(WH,GHx,GHy,GHz,Jx,Jy,Jz,cx,cy,cz);
	// 기존의  CG routine
	ARRAY<double>& x = pressure;
	ARRAY<double>  r(num_cells());
	ARRAY<double>  p(num_cells());
	ARRAY<double> Ap(num_cells());
	for (int i = 0; i < cells.size(); i++) {
		x[i] = 0; r[i] = RHS[i]; p[i] = RHS[i];
	}
	projection_sumzero(r);
	projection_sumzero(p);
	double rr = inner_product(r, r);
	double rr0 = rr; int cnt = 0;
	ARRAY<double>   Gp(num_faces());
	ARRAY<double>  HGp(num_faces());
	ARRAY<double> WHGp(num_faces());

	while (rr > rr0 * 1e-15) {
		// Ap=DWHGp
		calculate_G(p, Gp);
		calculate_H(H, Gp, HGp);
		calculate_W(HGp, WHGp);
		calculate_D(WHGp, Ap);
		// v = minv * <GH,p>
		double vx = minv * inner_product(p,GHx);
		double vy = minv * inner_product(p,GHy);
		double vz = minv * inner_product(p,GHz);
		// w = Iinv * <J,p>
		double tempx = inner_product(p, Jx);
		double tempy = inner_product(p, Jy);
		double tempz = inner_product(p, Jz);
		double wx = Iinv[0][0]*tempx + Iinv[0][1]*tempy + Iinv[0][2]*tempz;
		double wy = Iinv[1][0]*tempx + Iinv[1][1]*tempy + Iinv[1][2]*tempz;
		double wz = Iinv[2][0]*tempx + Iinv[2][1]*tempy + Iinv[2][2]*tempz;
		// A = DWHG - GH*GH - J*J
		//
		for (int i=0;i<num_cells();i++)
			Ap[i]-=GHx[i]*vx+GHy[i]*vy+GHz[i]*vz
			       +Jx[i]*wx+ Jy[i]*wy+ Jz[i]*wz;

		projection_sumzero(Ap);
		double pAp = inner_product(p, Ap);
		double alpha = rr / pAp;

		for (int i = 0; i < cells.size(); i++) {
			x[i] += alpha * p[i];
			r[i] -= alpha * Ap[i];
		}

		projection_sumzero(r);
		double rr_next = inner_product(r, r);
		double beta_ = rr_next / rr;
		rr = rr_next;

		for (int i = 0; i < cells.size(); i++)
			p[i] = r[i] + beta_ * p[i];

		if (cnt % 10 == 0)
			printf("cnt=%d, rr=%e\n", cnt, rr);
		cnt++;
	}
}
/*
void OcTree::visualize() {
	Display disp;
	disp.set_physical_coordinates(0, N, 0, N);
	disp.make_it_white();
	for (int index = 0; index < cells.size(); index++) {
		OcCell cell = cells[index];
		if (cell.leaf) {
			disp.draw_line(cell.imin, cell.jmin, cell.imax, cell.jmin, 0x000000);
			disp.draw_line(cell.imin, cell.jmin, cell.imin, cell.jmax, 0x000000);
			disp.draw_line(cell.imax, cell.jmax, cell.imax, cell.jmin, 0x000000);
			disp.draw_line(cell.imax, cell.jmax, cell.imin, cell.jmax, 0x000000);
		}
	}
	disp.save_bitmap_file();
}

void OcTree::split_all_cells() {
	int size_c = cells.size();
	for (int index_c = 0; index_c < size_c; index_c++) {
		if (cells[index_c].leaf)
			split_c(index_c);
	}
}






void OcTree::visualize( const ARRAY<double>& U, double scale, double (*phi)(double,double) ) const {
	Display disp;
	disp.set_physical_coordinates(xmin, xmax, ymin, ymax);
	disp.make_it_white();
	for (int index = 0; index < cells.size(); index++) {
		OcCell C = cells[index];
		// drawing cell
		if (C.leaf) {
			double xmin = x_fr_i(C.imin), xmax=x_fr_i(C.imax);
			double ymin = y_fr_j(C.jmin), ymax=y_fr_j(C.jmax);
			disp.draw_line(xmin, ymin, xmax, ymin, 0x000000);
			disp.draw_line(xmin, ymin, xmin, ymax, 0x000000);
			disp.draw_line(xmax, ymax, xmax, ymin, 0x000000);
			disp.draw_line(xmax, ymax, xmin, ymax, 0x000000);

			double u = (U[C.index_face[0][0]]+U[C.index_face[0][1]])/2;
			double v = (U[C.index_face[1][0]]+U[C.index_face[1][1]])/2;
			disp.draw_line((xmin+xmax)/2        ,(ymin+ymax)/2        ,
			               (xmin+xmax)/2+u*scale,(ymin+ymax)/2+v*scale, 0xff0000); 
			
			double phi00=phi(xmin,ymin), phi10=phi(xmax,ymin);
			double phi01=phi(xmin,ymax), phi11=phi(xmax,ymax);
			double x[4],y[4]; int count=0;
			if(phi00*phi10<0){ x[count]=(xmin*phi10-xmax*phi00)/(phi10-phi00); y[count]=ymin; count++;}
			if(phi00*phi01<0){ y[count]=(ymin*phi01-ymax*phi00)/(phi01-phi00); x[count]=xmin; count++;}
			if(phi10*phi11<0){ y[count]=(ymin*phi11-ymax*phi10)/(phi11-phi10); x[count]=xmax; count++;}
			if(phi01*phi11<0){ x[count]=(xmin*phi11-xmax*phi01)/(phi11-phi01); y[count]=ymax; count++;}
			if(count==2)
				disp.draw_line(x[0],y[0],x[1],y[1],0xff0000);
		}
	}
	disp.save_bitmap_file();
}

double OcTree::calculate_G(const ARRAY<double>& u, int index_f) const {

	const OcFace& F = faces[index_f];
	if (F.leaf) {
		if (F.index_cell[0] == -1 || F.index_cell[1] == -1)
			return 0;
		else if (F.imax == F.imin) {
			int c0 = F.index_cell[0];
			int c1 = F.index_cell[1];
			return (u[c1] - u[c0]) / (xsize(c1) * .5 + xsize(c0) * .5);
		}
		else if (F.jmax == F.jmin) {
			int c0 = F.index_cell[0];
			int c1 = F.index_cell[1];
			return(u[c1] - u[c0]) / (ysize(c1) * .5 + ysize(c0) * .5);
		}
		else {
			assert(false);
			return 0;
		}
	}
	else
		assert(false);
}

void OcTree::calculate_G(const ARRAY<double>& u, ARRAY<double>& Gu) const{
	for(int i=0;i<faces.size();i++){
		const OcFace& F=faces[i];
		if(F.leaf)
			Gu[i] = calculate_G(u,i);
		else 
			Gu[i]=0;
	}
}

double OcTree::calculate_flux(const ARRAY<double>& f, int index_f) const{
	const OcFace& F = faces[index_f];
	if(F.leaf){
		if(F.imax==F.imin)
			return f[index_f]*(y_fr_j(F.jmax)-y_fr_j(F.jmin));
		else 
			return f[index_f]*(x_fr_i(F.imax)-x_fr_i(F.imin));
	}
	else{
		return calculate_flux(f,F.index_child[0])
		     + calculate_flux(f,F.index_child[1]);
	}
}

void OcTree::calculate_D(const ARRAY<double>&f, ARRAY<double>& D) const{
		for(int index_d = 0; index_d<cells.size(); index_d++){
		const OcCell&C = cells[index_d];
		if (C.leaf==true){
			D[index_d]=calculate_flux(f,C.index_face[0][1])
					- calculate_flux(f, C.index_face[0][0])
					+ calculate_flux(f, C.index_face[1][1])
					- calculate_flux(f, C.index_face[1][0]);
		}
		else D[index_d]=0;
	}
}






void OcTree::solve_Laplacian( const ARRAY<double>& beta, const ARRAY<double>& f, ARRAY<double>& u ) const{
	ARRAY<double>& x = u;
	ARRAY<double>  r(cells.size());
	ARRAY<double>  p(cells.size());
	ARRAY<double> Ap(cells.size());
	for(int i=0;i<cells.size();i++){
		x[i] = 0; r[i]=f[i]; p[i]=f[i]; 
	}

	projection_sumzero(r);
	projection_sumzero(p);
	double rr = inner_product(r,r);
	double rr0 = rr; int cnt=0;
	ARRAY<double> Gp(faces.size());

	while(rr>rr0*1e-15){
		for(int index_f=0;index_f<faces.size();index_f++){
			if(faces[index_f].leaf){
				Gp[index_f] = beta[index_f]*calculate_G(p,index_f);
			}
		}
		calculate_D(Gp,Ap);

		projection_sumzero(Ap);
		double pAp  = inner_product(p,Ap);
		double alpha = rr/pAp;

		for(int i=0;i<cells.size();i++){
			x[i]+= alpha* p[i];
			r[i]-= alpha*Ap[i];
		}

		projection_sumzero(r);
		double rr_next = inner_product(r,r);
		double beta_ = rr_next/rr;
		rr = rr_next;

		for(int i=0;i<cells.size();i++)
			p[i] = r[i] +beta_*p[i];

		if(cnt%10==0)
			printf("cnt=%d, rr=%f\n", cnt, rr);
		cnt++;
	}
}
// Laplacian : 2nd order method
void OcTree::calculate_W(const ARRAY<double>& f, int index_F, double& Wf, double& delta) const{
	const OcFace& F= faces[index_F];
	if(F.leaf) {
		Wf=f[index_F];
		int c0=F.index_cell[0];
		int c1=F.index_cell[1];
		if(c0==-1 || c1==-1) {	delta=1; return;	}
		if(F.imin==F.imax) delta=.5*(xsize(c0)+xsize(c1));
		else			   delta=.5*(ysize(c0)+ysize(c1));
		return;
	}
	double Wf0, delta0;
	calculate_W(f, F.index_child[0], Wf0, delta0);
	double Wf1, delta1;
	calculate_W(f, F.index_child[1], Wf1, delta1);
	Wf=(Wf0*delta0+Wf1*delta1)/(delta0+delta1);
	delta=(delta0+delta1)/2;
}

void OcTree::calculate_W(const ARRAY<double>& f, ARRAY<double>& Wf) const {
	for(int index_F=0; index_F<faces.size(); index_F++) {
		const OcFace& F = faces[index_F];
		int c0 = F.index_cell[0];
		int c1 = F.index_cell[1];
		if(c0==-1||c1==-1) Wf[index_F]=f[index_F];
		else{
			const OcCell& C0= cells[c0]; int S0= C0.imax-C0.imin;
			const OcCell& C1 = cells[c1]; int S1=C1.imax-C1.imin;
			if(C0.leaf||C1.leaf) {
				if(C0.leaf&&C1.leaf) {
					if(S0==S1) Wf[index_F]=f[index_F];
					else ;
				}
				else {
					if(S0==S1) {
						double delta;
						calculate_W(f, index_F, Wf[index_F], delta);
						set_W(index_F, Wf[index_F], Wf);
					}
					else ;
				}
			}
			else Wf[index_F]=f[index_F];
		}
	}
}

void OcTree::set_W(int index_F, double Wf_value, ARRAY<double>& Wf) const {
	const OcFace& F= faces[index_F];
	Wf[index_F]=Wf_value;
	if(F.leaf==false) {
		set_W(F.index_child[0], Wf_value, Wf);
		set_W(F.index_child[1], Wf_value, Wf);
	}
}

void OcTree::solve_Laplacian_1st(const ARRAY<double>& H, const ARRAY<double>& f, ARRAY<double>& u) const {
	ARRAY<double>& x = u;
	ARRAY<double>  r(cells.size());
	ARRAY<double>  p(cells.size());
	ARRAY<double> Ap(cells.size());
	for (int i = 0; i < cells.size(); i++) {
		x[i] = 0; r[i] = f[i]; p[i] = f[i];
	}

	projection_sumzero(r);
	projection_sumzero(p);
	double rr = inner_product(r, r);
	double rr0 = rr; int cnt = 0;
	ARRAY<double>  Gp(faces.size());
	ARRAY<double> HGp(faces.size());

	while (rr > rr0 * 1e-15) {
		// Ap=DWGp
		calculate_G( p,Gp);
		calculate_H(H,Gp,HGp);
		calculate_D(HGp, Ap);

		projection_sumzero(Ap);
		double pAp = inner_product(p, Ap);
		double alpha = rr / pAp;

		for (int i = 0; i < cells.size(); i++) {
			x[i] += alpha * p[i];
			r[i] -= alpha * Ap[i];
		}

		projection_sumzero(r);
		double rr_next = inner_product(r, r);
		double beta_ = rr_next / rr;
		rr = rr_next;

		for (int i = 0; i < cells.size(); i++)
			p[i] = r[i] + beta_ * p[i];

		if (cnt % 10 == 0)
			printf("cnt=%d, rr=%f\n", cnt, rr);
		cnt++;
	}
}



// irregular domain
void OcTree::calculate_H(const ARRAY<double>& f, const ARRAY<double>& H, ARRAY<double>& Hf )const{
	for(int index_F=0;index_F<faces.size();index_F++){
		const OcFace& F = faces[index_F];
		if(F.leaf)
			Hf[index_F] = f[index_F]*H[index_F];
		else
			Hf[index_F]=0;
	}
}

// FSI

void OcTree::calculate_GH_J(const ARRAY<double>& H, ARRAY<double>& GH_x, ARRAY<double>& GH_y, ARRAY<double>& J, const double cx, const double cy) const{
	for(int i=0; i<cells.size(); i++){
		const OcCell& C = cells[i];
		const OcFace& FR = faces[C.index_face[0][1]];
		const OcFace& FL = faces[C.index_face[0][0]];
		const OcFace& FU = faces[C.index_face[1][1]];
		const OcFace& FD = faces[C.index_face[1][0]];

		if(FR.leaf && FL.leaf && FU.leaf && FD.leaf){
		GH_x[i] = calculate_flux(H, C.index_face[0][1]) - calculate_flux(H, C.index_face[0][0]);
		GH_y[i] = calculate_flux(H, C.index_face[1][1]) - calculate_flux(H, C.index_face[1][0]);

		double x = (x_fr_i(C.imax)+x_fr_i(C.imin))/2;
		double y = (y_fr_j(C.jmax)+y_fr_j(C.jmin))/2;

		J[i] = (x-cx)*GH_y[i] - (y-cy)*GH_x[i];
		}
		else{
			GH_x[i]=0; GH_y[i]=0; J[i]=0;
		}
	}
}

void OcTree::solve_FSI(const ARRAY<double>& H, const ARRAY<double>& RHS, 
                   double minv, double Iinv, double cx, double cy, ARRAY<double>& pressure) const{
	ARRAY<double> GH_x(cells.size());
	ARRAY<double> GH_y(cells.size());
	ARRAY<double> J   (cells.size()); calculate_GH_J(H,GH_x,GH_y,J,cx,cy);
	// 기존의  CG routine
	ARRAY<double>& x = pressure;
	ARRAY<double>  r(cells.size());
	ARRAY<double>  p(cells.size());
	ARRAY<double> Ap(cells.size());
	for (int i = 0; i < cells.size(); i++) {
		x[i] = 0; r[i] = RHS[i]; p[i] = RHS[i];
	}
	projection_sumzero(r);
	projection_sumzero(p);
	double rr = inner_product(r, r);
	double rr0 = rr; int cnt = 0;
	ARRAY<double>  Gp(faces.size());
	ARRAY<double> HGp(faces.size());

	while (rr > rr0 * 1e-15) {
		// Ap=DWGp
		calculate_G(p, Gp);
		calculate_H(H, Gp, HGp);
		calculate_D(HGp, Ap);
		//
		double vx,vy,w;
		vx = minv*inner_product(p,GH_x);
		vy = minv*inner_product(p,GH_y);
		w  = Iinv*inner_product(p,J);
		//
		for(int i=0;i<cells.size();i++)
			Ap[i]=-Ap[i]+GH_x[i]*vx+GH_y[i]*vy+J[i]*w;

		projection_sumzero(Ap);
		double pAp = inner_product(p, Ap);
		double alpha = rr / pAp;

		for (int i = 0; i < cells.size(); i++) {
			x[i] += alpha * p[i];
			r[i] -= alpha * Ap[i];
		}

		projection_sumzero(r);
		double rr_next = inner_product(r, r);
		double beta_ = rr_next / rr;
		rr = rr_next;

		for (int i = 0; i < cells.size(); i++)
			p[i] = r[i] + beta_ * p[i];

		if (cnt % 10 == 0)
			printf("cnt=%d, rr=%f\n", cnt, rr);
		cnt++;
	}
}

void OcTree::solve_FSI_2nd(const ARRAY<double>& H, const ARRAY<double>& RHS,
	double minv, double Iinv, double cx, double cy, ARRAY<double>& pressure) const {
	ARRAY<double> GH_x(cells.size());
	ARRAY<double> GH_y(cells.size());
	ARRAY<double> J(cells.size()); calculate_GH_J(H, GH_x, GH_y, J, cx, cy);
	// 기존의  CG routine
	ARRAY<double>& x = pressure;
	ARRAY<double>  r(cells.size());
	ARRAY<double>  p(cells.size());
	ARRAY<double> Ap(cells.size());
	for (int i = 0; i < cells.size(); i++) {
		x[i] = 0; r[i] = RHS[i]; p[i] = RHS[i];
	}
	projection_sumzero(r);
	projection_sumzero(p);
	double rr = inner_product(r, r);
	double rr0 = rr; int cnt = 0;
	ARRAY<double>   Gp(faces.size());
	ARRAY<double>  HGp(faces.size());
	ARRAY<double> WHGp(faces.size());

	while (rr > rr0 * 1e-15) {
		// Ap=DWGp
		calculate_G(p, Gp);
		calculate_H(H, Gp, HGp);
		calculate_W(HGp, WHGp);
		calculate_D(WHGp, Ap);
		//
		double vx, vy, w;
		vx = minv * inner_product(p, GH_x);
		vy = minv * inner_product(p, GH_y);
		w  = Iinv * inner_product(p, J   );
		//
		for (int i = 0; i < cells.size(); i++)
			Ap[i] = -Ap[i] + GH_x[i] * vx + GH_y[i] * vy + J[i] * w;

		projection_sumzero(Ap);
		double pAp = inner_product(p, Ap);
		double alpha = rr / pAp;

		for (int i = 0; i < cells.size(); i++) {
			x[i] += alpha * p[i];
			r[i] -= alpha * Ap[i];
		}

		projection_sumzero(r);
		double rr_next = inner_product(r, r);
		double beta_ = rr_next / rr;
		rr = rr_next;

		for (int i = 0; i < cells.size(); i++)
			p[i] = r[i] + beta_ * p[i];

		if (cnt % 10 == 0)
			printf("cnt=%d, rr=%f\n", cnt, rr);
		cnt++;
	}
}

void OcTree::calculate_Gmatrix(int index_f, double scale, SparseSymmMatrix& A, int i){
	const OcFace& F = faces[index_f];
	if(F.leaf){
		int c0 = F.index_cell[0];
		int c1 = F.index_cell[1];
		if(c0 != -1 && c1 != -1){
			double Delta = (F.imin == F.imax)? 0.5*(xsize(c0) + xsize(c1)) : 0.5 * (ysize(c0) + ysize(c1));
			A.add_element(i,c0,-scale/Delta);
			A.add_element(i,c1,+scale/Delta);
		}
	}
}

double OcTree::calculate_WDelta(int index_F) const{
	const OcFace& F = faces[index_F];
	if(F.leaf){
		int c0 = F.index_cell[0];
		int c1 = F.index_cell[1];
		if(c0 == -1 || c1 == -1) 
			return 1;
		if(F.imin==F.imax)
			return 0.5*(xsize(c0) + xsize(c1));
		else
			return 0.5*(ysize(c0) + ysize(c1));
	}
	else
		return 0.5*(calculate_WDelta(F.index_child[0]) + calculate_WDelta(F.index_child[1]));
}

void OcTree::calculate_WHGmatrix(int index_f, double scale, const ARRAY<double>& H, SparseSymmMatrix& A, int i){
	const OcFace& F = faces[index_f];
	if(F.leaf){
		calculate_Gmatrix(index_f,scale*H[index_f],A,i);
	}
	else {
		double Delta0 = calculate_WDelta(F.index_child[0]);
		double Delta1 = calculate_WDelta(F.index_child[1]);
		calculate_WHGmatrix(F.index_child[0], scale*Delta0/(Delta0+Delta1), H, A, i);
		calculate_WHGmatrix(F.index_child[1], scale*Delta1/(Delta0+Delta1), H, A, i);
	}
}

void OcTree::calculate_DWHGmatrix(const ARRAY<double>& H, SparseSymmMatrix & A){
	A.set_dimension(num_cells());
	for(int i=0; i<num_cells(); i++){
		const OcCell& C = cells[i];
		if(C.leaf){
			for(int type_ij=0; type_ij<2; type_ij++)
			for(int k=0; k<2; k++){
				int index_F = C.index_face[type_ij][k];
				OcFace F = faces[index_F];
				if(F.index_cell[0] != -1 && F.index_cell[1] != -1){
					int iprime = F.index_cell[k];
					OcCell Cprime = cells[iprime];
					int Cprimesize = (type_ij==0)? (Cprime.jmax - Cprime.jmin):(Cprime.imax - Cprime.imin);
					while(F.size() < Cprimesize){
						index_F = F.index_prnt;
						F = faces[index_F];
					}
					double Delta = (type_ij==0)? ysize(i) : xsize(i);
					calculate_WHGmatrix(index_F, Delta * ((k==1)? 1:(-1)), H, A, i);
				}
			}
		}
		else
			A.add_element(i,i,1);
	}
}
*/
