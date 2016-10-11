#include <fstream>
#include <iostream>
using namespace std;

class PNP
{
	friend ostream & operator<<(ostream &, const PNP &);
public:
	PNP(){};
	void SetDomain(const int &, const double &, const double &);
	void SetInitial(const int &, double *, double **, double **);
	void SetLaplacian();
	void ComputePotential();
	void ComputePhinew();
	void SetTime(const int &, const double &);
	void UpdateDensity();
	void Update();
	
	void Solve();
private:
	int _s; //species
	// space discretization
	int _n;
	double _L_min;
	double _L_max;
	double _h;
	double * _x; // x coordinate
	// time discretization
	int _m;
	double _tau;
	double _lambda;
	double _T;
	double * _sigma_l;//surface charge density
	double * _sigma_r;
	
	double * _z;//valence
	
	double ** _density;// density of each species at different position
	double ** _density_new;
	double ** _theta;   //temperature
	double ** _theta_new;
	double * _phi;//potential
	double * _phi_new;
	//tridiagnal matrix to solve _phi
	double * _a;
	double * _b;
	double * _c;
	//temp
	double * _rhs;
};