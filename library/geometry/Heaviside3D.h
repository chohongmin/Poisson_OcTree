#ifndef HEAVISIDE3D_H
#define HEAVISIDE3D_H

#include "../ARRAY3D.h"
#include "../uniform_grid/GRID3D.h"
/*----------------------------------------------------------------------------*\
 * calculate the Heaviside function 
\*----------------------------------------------------------------------------*/
static void calculate_H_Towers(const ARRAY3DwithBUFFER<double>& Phi, 
							         ARRAY3DwithBUFFER<double>& Hiphjphkph, double dx, 
																			double dy, 
																			double dz)
{
	int isize = Phi.isize; 
	int jsize = Phi.jsize; 
	int ksize = Phi.ksize;
	Hiphjphkph.resize(isize-1,jsize-1,ksize-1,2); 
	// 
	for(int i=0;i<isize-1;i++)
	for(int j=0;j<jsize-1;j++)
	for(int k=0;k<ksize-1;k++)
	{
		int count_positive=0;
		double phi_imhjmhkmh = Phi(i  ,j  ,k  ); if(phi_imhjmhkmh>0) count_positive++;
		double phi_imhjmhkph = Phi(i  ,j  ,k+1); if(phi_imhjmhkph>0) count_positive++;
		double phi_imhjphkmh = Phi(i  ,j+1,k  ); if(phi_imhjphkmh>0) count_positive++;
		double phi_imhjphkph = Phi(i  ,j+1,k+1); if(phi_imhjphkph>0) count_positive++;
		double phi_iphjmhkmh = Phi(i+1,j  ,k  ); if(phi_iphjmhkmh>0) count_positive++;
		double phi_iphjmhkph = Phi(i+1,j  ,k+1); if(phi_iphjmhkph>0) count_positive++;
		double phi_iphjphkmh = Phi(i+1,j+1,k  ); if(phi_iphjphkmh>0) count_positive++;
		double phi_iphjphkph = Phi(i+1,j+1,k+1); if(phi_iphjphkph>0) count_positive++;

			 if( count_positive==0 ) Hiphjphkph(i,j,k)=0;
		else if( count_positive==8 ) Hiphjphkph(i,j,k)=1;  
		else
		{
			double phi_iph = .25*(phi_iphjmhkmh+phi_iphjmhkph+phi_iphjphkmh+phi_iphjphkph);
			double phi_imh = .25*(phi_imhjmhkmh+phi_imhjmhkph+phi_imhjphkmh+phi_imhjphkph);
			double phi_jph = .25*(phi_imhjphkmh+phi_imhjphkph+phi_iphjphkmh+phi_iphjphkph);
			double phi_jmh = .25*(phi_imhjmhkmh+phi_imhjmhkph+phi_iphjmhkmh+phi_iphjmhkph);
			double phi_kph = .25*(phi_imhjmhkph+phi_imhjphkph+phi_iphjmhkph+phi_iphjphkph);
			double phi_kmh = .25*(phi_imhjmhkmh+phi_imhjphkmh+phi_iphjmhkmh+phi_iphjphkmh);

			double r_iph = (phi_iph>0) ? phi_iph : 0;
			double r_imh = (phi_imh>0) ? phi_imh : 0;
			double r_jph = (phi_jph>0) ? phi_jph : 0;
			double r_jmh = (phi_jmh>0) ? phi_jmh : 0;
			double r_kph = (phi_kph>0) ? phi_kph : 0;
			double r_kmh = (phi_kmh>0) ? phi_kmh : 0;

			Hiphjphkph(i,j,k) =((r_iph-r_imh)/dx*(phi_iph-phi_imh)/dx
				              + (r_jph-r_jmh)/dy*(phi_jph-phi_jmh)/dy
							  + (r_kph-r_kmh)/dz*(phi_kph-phi_kmh)/dz)/
						       ((phi_iph-phi_imh)/dx*(phi_iph-phi_imh)/dx
				              +(phi_jph-phi_jmh)/dy*(phi_jph-phi_jmh)/dy
							  +(phi_kph-phi_kmh)/dz*(phi_kph-phi_kmh)/dz+1E-10);
		}
	}
	Hiphjphkph.extrapolate_constant();
}

#endif
