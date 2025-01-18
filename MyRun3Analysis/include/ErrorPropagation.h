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

double Error_SCnm(double N, double eN, double N1, double eN1, double N2, double eN2){
    //SC = (N-N1*N2)

    double err = sqrt(
        pow(eN,2) + pow(eN1*N2,2) + pow(eN2*N1,2)
    );
    if(err<=0||err>5)printf("Warning: err for SC is %f\n",err);
    return err;
}

double Error_SCklm(double N_klm, double eN_klm, double N_kl, double eN_kl, 
double N_km, double eN_km, double N_lm, double eN_lm,
double N_k, double eN_k, double N_l, double eN_l, double N_m, double eN_m
){
    //SC(k,l,m) = N_klm - N_kl*N_m - N_km*N_l - N_lm*N_k + 2*N_k*N_l*N_m

    double err = sqrt(
        pow(eN_klm,2)+pow(eN_kl*N_m,2)+pow(eN_km*N_l,2)+pow(eN_lm*N_k,2)
        +pow((2*N_l*N_m-N_lm)*eN_k,2)
        +pow((2*N_k*N_m-N_km)*eN_l,2)
        +pow((2*N_k*N_l-N_kl)*eN_m,2)
    );
    if(err<=0||err>5)printf("Warning: err for SC is %f\n",err);
    return err;
}

double Error_NSCklm(double N_klm, double eN_klm, double N_kl, double eN_kl, 
double N_km, double eN_km, double N_lm, double eN_lm,
double N_k, double eN_k, double N_l, double eN_l, double N_m, double eN_m,
double G_k, double eG_k, double G_l, double eG_l, double G_m, double eG_m
){
    //SC(k,l,m) = N_klm - N_kl*N_m - N_km*N_l - N_lm*N_k + 2*N_k*N_l*N_m
    double scklm = N_klm - N_kl*N_m - N_km*N_l - N_lm*N_k + 2*N_k*N_l*N_m;
    double normal = G_k*G_l*G_m;
    double err = sqrt(
        pow(eN_klm/normal,2)+pow(eN_kl*N_m/normal,2)+pow(eN_km*N_l/normal,2)+pow(eN_lm*N_k/normal,2)
        +pow((2*N_l*N_m-N_lm)*eN_k/normal,2)
        +pow((2*N_k*N_m-N_km)*eN_l/normal,2)
        +pow((2*N_k*N_l-N_kl)*eN_m/normal,2)
        +pow(scklm*eG_k/(normal*G_k),2)
        +pow(scklm*eG_l/(normal*G_l),2)
        +pow(scklm*eG_m/(normal*G_m),2)
    );
    if(err<=0||err>5)printf("Warning: err for NSC is %f\n",err);
    return err;
}

double Error_NSC(double N, double eN, double N1, double eN1, double N2, double eN2, double D1, double eD1, double D2, double eD2){
    //NSC = (N-N1*N2)/(D1*D2)

    double err = sqrt(
        pow(eN/(D1*D2),2)
        +pow(eN1*N2/(D1*D2),2)
        +pow(eN2*N1/(D1*D2),2)
        +pow(eD1*(N-N1*N2)/(D1*D1*D2),2)
        +pow(eD2*(N-N1*N2)/(D1*D2*D2),2)
    );
    if(err<=0||err>5)printf("Warning: err for NSC is %f\n",err);
    return err;
}

double Error_CN10(double cor10, double cor10e, double cor8, double cor8e, double cor6, double cor6e, double cor4, double cor4e, double cor2, double cor2e){
    double parts[5];
    parts[0] = cor10e;
    parts[1] = -25.*cor2*cor8e;
    parts[2] = -100.*cor4*cor6e + 400.*cor2*cor2*cor6e;
    parts[3] = -100.*cor6*cor4e + 900.*2.*cor4*cor2*cor4e - 3600.*cor2*cor2*cor2*cor4e;
    parts[4] = -25*cor8*cor2e + 400*cor6*2.*cor2*cor2e - 900.*cor4*cor4*cor2e - 3600.*cor4*3.*cor2*cor2*cor2e + 2880.*5.*cor2*cor2*cor2*cor2*cor2e;
    double retval = 0;
    for (int i = 0; i < 5; i++)
        retval += TMath::Power(parts[i], 2);
    return TMath::Sqrt(retval);
}

#endif
