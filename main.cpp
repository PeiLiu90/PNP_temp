#include "PNP.h"
#include "math.h"
#include <iostream>
using namespace std;

const int s=2;
const int n=100;
//const int m=n*n*10;
const double Ls=-3;
const double Lg=3;
const double h = (Lg-Ls)/double(n);
double * z = new double [s];
double ** c = new double * [s];
double ** theta = new double * [s];
const double T = 5.;
const double tau = h*h*0.01;
const int m=floor(T/tau);

int main()
{
	for(int i=0;i<s;i++)
	{
		c[i]=new double [n+1];
		theta[i]=new double [n+1];
		z[i]=(i*2-1);//+1 and -1
		for(int j=0;j<=n;j++)
		{
			c[i][j]=1.;
			theta[i][j]=0.5+ double(j)*1./double(n);
		}
	}
	
	PNP test;
	test.SetDomain(n,Ls,Lg);
	test.SetInitial(s,z,c,theta);
	test.SetLaplacian();
	test.SetTime(m,tau);
	
	test.Solve();
	
	
	ofstream os("density.txt");
	os<<test;
	os.close();

	return 0;
}