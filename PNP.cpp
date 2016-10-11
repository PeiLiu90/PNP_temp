#include "PNP.h"
#include "ChaseMethod.h"
#include "math.h"
#include <iostream>
using namespace std;

ostream & operator<<(ostream & os, const PNP & a)
{
	for(int i=0;i<=a._n;i++)
	{
		os<<a._x[i];
		for(int j=0;j<a._s;j++)
		{
			os<<"  "<<a._density[j][i]<<"  "<<a._theta[j][i];
		}
		os<<"  "<<a._phi[i]<<endl;
	}
	return os;
}

void PNP::SetDomain(const int & n, const double & l1, const double & l2)
{
	_n=n;
	_x = new double [n+1];
	_phi = new double [n+1];
	_phi_new = new double [n+1];
	_rhs = new double [n+1];
	_a = new double [n+1];
	_b = new double [n+1];
	_c = new double [n+1];
	_L_min = l1;
	_L_max = l2;
	_h = (l2-l1)/n;
	for(int i=0;i<=n;i++)
	{
		_x[i]= l1 + double(i)*_h;
	}
}

void PNP::SetInitial(const int & s, double * z, double ** c, double ** theta)
{
	_s = s;
	_z = new double [s];
	_density = new double * [s];
	_density_new = new double * [s];
	_theta = new double * [s];
	_theta_new = new double * [s];
	for(int i=0;i<_s;i++)
	{
		_density[i]=new double [_n+1];
		_density_new[i]=new double [_n+1];
		_theta[i] = new double [_n+1];
		_theta_new[i] = new double [_n+1];
		_z[i]=z[i];
		for(int j=0;j<=_n;j++)
		{
			_density[i][j]=c[i][j];
			_theta[i][j]=theta[i][j];
		}
	}
}

void PNP::SetLaplacian()
{
	_a[0]=_c[0]=0.;
	_b[0]=1.;
	for(int i=1;i<_n;i++)
	{
		_a[i]=-1.;
		_b[i]=2.;
		_c[i]=-1.;
	}
	_a[_n]=_c[_n]=0.;
	_b[_n]=1.;
}

void PNP::ComputePotential()
{
	_phi[0]=-1.;
	for(int i=1;i<_n;i++)
	{
		_phi[i]=0.;
		for(int j=0;j<_s;j++)
		{
			_phi[i]+=_density[j][i]*_z[j]*_h*_h;
		}
	}
	_phi[_n]=1.;
	ChaseMethod(_a,_b,_c,_n+1,_phi);
}

void PNP::ComputePhinew()
{
	_phi_new[0]=-1.;
	for(int i=1;i<_n;i++)
	{
		_phi_new[i]=0.;
		for(int j=0;j<_s;j++)
		{
			_phi_new[i]+=_density_new[j][i]*_z[j];
		}
	}
	_phi_new[_n]=1.;
	ChaseMethod(_a,_b,_c,_n+1,_phi_new);
}

void PNP::SetTime(const int & m, const double & tau)
{
	_m=m;
	_tau=tau;
	_lambda=_tau/_h/_h;
	_T=_tau*double(_m);
}

void PNP::UpdateDensity()
{
	for(int i=0;i<_s;i++)
	{
		_density_new[i][0]=_density[i][0]+_lambda/2. * (_theta[i][1]+_theta[i][0])*( (_density[i][1]-_density[i][0])+ _z[i]*(_density[i][1]+_density[i][0])*(_phi[1]-_phi[0])/2. );
		_theta_new[i][0]=_theta[i][0] + _lambda * 4. * sqrt(_theta[i][0]) /_density[i][0] * ( _theta[i][1] - _theta[i][0]);
		for(int j=1;j<_n;j++)
		{
			_density_new[i][j]=_density[i][j]+_lambda/2. * ( (_theta[i][j+1]+_theta[i][j])*( (_density[i][j+1]-_density[i][j])+ _z[i]* (_density[i][j+1]+_density[i][j])*(_phi[j+1]-_phi[j])/2. )
				                         - (_theta[i][j]+_theta[i][j-1])* ( (_density[i][j]-_density[i][j-1]) + _z[i]* (_density[i][j]+_density[i][j-1])*(_phi[j]-_phi[j-1])/2. ) );
			_theta_new[i][j]=_theta[i][j] + _lambda * 4. * sqrt(_theta[i][j]) /_density[i][j] * ( _theta[i][j+1]+_theta[i][j-1] - 2. *_theta[i][j])
										  + _lambda * _theta[i][j]/4. * ( (_density[i][j+1]-_density[i][j-1])/_density[i][j] + _z[i]* (_phi[j+1]-_phi[j-1]) )* (_theta[i][j+1]-_theta[i][j-1])
										  + _lambda * pow(_theta[i][j],1.5)/_density[i][j] * ( (_density[i][j+1]-_density[i][j-1])/_density[i][j] + _z[i]* (_phi[j+1]-_phi[j-1]) ) *( (_density[i][j+1]-_density[i][j-1])/_density[i][j] + _z[i]* (_phi[j+1]-_phi[j-1]) )  ;
			
		}
		_density_new[i][_n]=_density[i][_n]-_lambda/2. * (_theta[i][_n]+_theta[i][_n-1])*( (_density[i][_n]-_density[i][_n-1]) +_z[i]*(_density[i][_n]+_density[i][_n-1])*(_phi[_n]-_phi[_n-1])/2. );\
		_theta_new[i][_n]=_theta[i][_n] + _lambda * 4. * sqrt(_theta[i][_n]) /_density[i][_n] * ( _theta[i][_n-1] - _theta[i][_n]);
	}
	for(int i=0;i<_s;i++)
	{
		for(int j=0;j<=_n;j++)
		{
			_density[i][j]=_density_new[i][j];
			_theta[i][j]=_theta_new[i][j];
		}
	}
}

void PNP::Update()
{
	for(int i=0;i<_s;i++)
	{
		for(int j=0;j<=_n;j++)
		{
			_density_new[i][j]=_density[i][j];
			_theta_new[i][j]=_theta[i][j];
		}
	}
	ComputePhinew();
	for(int j=0;j<=_n;j++)
	{
		_phi[j]=_phi_new[j];
	}
	for(int k=0;k<2;k++)
	{
		for(int i=0;i<_s;i++)
		{
			_rhs[0]=_tau/_h/(4.*_h) * ( (_theta[i][1]+_theta[i][0])*(_density[i][1]-_density[i][0])+ (_theta_new[i][1]+_theta_new[i][0])*(_density_new[i][1]-_density_new[i][0])
										+_z[i]*(_density[i][1]+_density[i][0])*(_phi[1]-_phi[0]) + _z[i]*(_density_new[i][1]+_density_new[i][0])*(_phi_new[1]-_phi_new[0]) );
			for(int j=1;j<_n;j++)
			{
				_density_new[i][j]=_density[i][j]+_tau/_h/(2.*_h) * ( (_theta[i][j+1]+_theta[i][j])*(_density[i][j+1]-_density[i][j])+ _z[i]* (_density[i][j+1]+_density[i][j])*(_phi[j+1]-_phi[j])
					                         - (_theta[i][j]+_theta[i][j-1])*(_density[i][j]-_density[i][j-1]) - _z[i]* (_density[i][j]+_density[i][j-1])*(_phi[j]-_phi[j-1]) );
			}
			_density_new[i][_n]=_density[i][_n]-_tau/_h/(2.*_h) * ( (_theta[i][_n]+_theta[i][_n-1])*(_density[i][_n]-_density[i][_n-1]) +_z[i]*(_density[i][_n]+_density[i][_n-1])*(_phi[_n]-_phi[_n-1]) );
		}
	}
}
	
	
	
void PNP::Solve()
{
	for(int i=0;i<_m;i++)
	{
		cout<<i<<endl;
		
		ComputePotential();
		UpdateDensity();
	}
}



















