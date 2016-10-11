/*
This file is used to solve the equation Ax=f,
with A is tridiagonale matrix:
b a 0 0 0 0
c b a 0 0 0
0 c b a 0 0
0 0 c b a 0
0 0 0 c b a
0 0 0 0 c b

or

b a 0 0 0 c
c b a 0 0 0
0 c b a 0 0
0 0 c b a 0
0 0 0 c b a
a 0 0 0 c b
*/


#ifndef CHASEMETHOD_H
#define CHASEMETHOD_H


//for the first type

void ChaseMethod(const double * a, const double * b, const double * c, const int N, double * f)
{
    double *beta=new double [N];
    beta[0]=a[0]/b[0];
    f[0]=f[0]/b[0];
    for(int i=1;i<N;i++)
    {
        beta[i]=a[i]/(b[i]-c[i]*beta[i-1]);
        f[i]=(f[i]-f[i-1]*c[i])/(b[i]-c[i]*beta[i-1]);
    }
    for(int i=N-2;i>=0;i--)
    {
        f[i]=f[i]-beta[i]*f[i+1];
    }
    delete [] beta;
}

void ChaseMethod(const double & a, const double & b, const double & c, const int N, double * f)
{
    double *beta=new double [N];
    beta[0]=a/b;
    f[0]=f[0]/b;
    for(int i=1;i<N;i++)
    {
        beta[i]=a/(b-c*beta[i-1]);
        f[i]=(f[i]-f[i-1]*c)/(b-c*beta[i-1]);
    }
    for(int i=N-2;i>=0;i--)
    {
        f[i]=f[i]-beta[i]*f[i+1];
    }
    delete [] beta;
}

//for the second type

void ChaseMethod2(const double * a, const double * b, const double * c, const int N, double * f)
{
    double * beta=new double [N-1];
    double * t=new double [N-2];//the last column is nonzero
    double s=a[N-1];//the last row has one nonzero element except c[N-1], b[N-1]
    double temp=b[N-1];//b[N-1] has changed value
    beta[0]=a[0]/b[0];
    f[0]=f[0]/b[0];
    t[0]=c[0]/b[0];
    for(int i=1;i<N-2;i++)
    {
        temp-=s*t[i-1];
        f[N-1]-=f[i-1]*s;
        s=-beta[i-1]*s;
        beta[i]=a[i]/(b[i]-c[i]*beta[i-1]);
        t[i]=-c[i]*t[i-1]/(b[i]-c[i]*beta[i-1]);
        f[i]=(f[i]-f[i-1]*c[i])/(b[i]-c[i]*beta[i-1]);
    }
    temp-=s*t[N-3];
    f[N-1]-=f[N-3]*s;
    s=c[N-1]-beta[N-3]*s;
    beta[N-2]=(a[N-2]-t[N-3]*c[N-2])/(b[N-2]-c[N-2]*beta[N-3]);
    f[N-2]=(f[N-2]-f[N-3]*c[N-2])/(b[N-2]-c[N-2]*beta[N-3]);
    f[N-1]=(f[N-1]-f[N-2]*s)/(temp-s*beta[N-2]);
    f[N-2]=f[N-2]-beta[N-2]*f[N-1];
    for(int i=N-3;i>=0;i--)
    {
        f[i]=f[i]-beta[i]*f[i+1]-t[i]*f[N-1];
    }
    delete [] beta;
    delete [] t;
}

void ChaseMethod2(const double a, const double b, const double c, const int N, double * f)
{
    double * beta=new double [N-1];
    double * t=new double [N-2];//the last column is nonzero
    double s=a;//the last row has one nonzero element except c[N-1], b[N-1]
    double temp=b;//b[N-1] has changed value
    beta[0]=a/b;
    f[0]=f[0]/b;
    t[0]=c/b;
    for(int i=1;i<N-2;i++)
    {
        temp-=s*t[i-1];
        f[N-1]-=f[i-1]*s;
        s=-beta[i-1]*s;
        beta[i]=a/(b-c*beta[i-1]);
        t[i]=-c*t[i-1]/(b-c*beta[i-1]);
        f[i]=(f[i]-f[i-1]*c)/(b-c*beta[i-1]);
    }
    temp-=s*t[N-3];
    f[N-1]-=f[N-3]*s;
    s=c-beta[N-3]*s;
    beta[N-2]=(a-t[N-3]*c)/(b-c*beta[N-3]);
    f[N-2]=(f[N-2]-f[N-3]*c)/(b-c*beta[N-3]);
    f[N-1]=(f[N-1]-f[N-2]*s)/(temp-s*beta[N-2]);
    f[N-2]=f[N-2]-beta[N-2]*f[N-1];
    for(int i=N-3;i>=0;i--)
    {
        f[i]=f[i]-beta[i]*f[i+1]-t[i]*f[N-1];
    }
    delete [] beta;
    delete [] t;
}
#endif
