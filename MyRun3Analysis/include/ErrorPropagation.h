#ifndef ERROTPROGATION
#define ERROTPROGATION
#include "math.h"

double Error_Ratio(double x, double ex, double y, double ey, double rho){
    //z=x/y
    //rho is correlation term, which is sqrt(N_y/N_x)
    //Cov(x,y) = rho*ex*ey
    double Contain = pow(ex/y,2)
        +pow(x*ey/(y*y),2)
        +2*(1./y)*(-1.*x/(y*y))*rho*ex*ey;
    if(Contain<0){
        rho=1;
        Contain = pow(ex/y,2)
        +pow(x*ey/(y*y),2)
        +2*(1./y)*(-1.*x/(y*y))*rho*ex*ey;
        if(Contain<0)Contain=4;
    }
    return sqrt(
        Contain
    );
}

double Error_vNL(double N, double eN, double D1, double eD1){
    //vNL = Numerator/sqrt(Denominator1)
    double err = sqrt(
        pow(eN/sqrt(D1),2)
        +pow(N*eD1/(2*pow(D1,3./2.)),2)
    );
    if(err<=0||err>5)printf("Warning: err for vNL is %f\n",err);
    return err;
}

double Error_Rho(double N, double eN, double D1, double eD1, double D2, double eD2){
    //Rho = Numerator/sqrt(Denominator1*Denominator2)
    double err = sqrt(
        pow(eN/sqrt(D1*D2),2)
        +pow(N*D2*eD1/(2*pow(D1*D2,3./2.)),2)
        +pow(N*D1*eD2/(2*pow(D1*D2,3./2.)),2)
    );
    if(err<=0||err>5)printf("Warning: err for Rho is %f\n",err);
    return err;
}

double Error_Chi(double N, double eN, double D1, double eD1){
    //Chi = Numerator/Denominator1
    double err = Error_Ratio(N,eN,D1,eD1,0.);
    if(err<=0||err>5)printf("Warning: err for Chi is %f\n",err);
    return err;
}

#endif
