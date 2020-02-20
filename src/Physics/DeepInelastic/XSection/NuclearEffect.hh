#ifndef _spectralFunT_H_
#define _spectralFunT_H_
#include "LHAPDF/LHAPDF.h"
#include <iostream>
#include <string>
#include <cmath>
#include <complex>
#include <algorithm>
#include "param.hh"
#include "Nucleon.hh"
using namespace LHAPDF;
using namespace Param;

class SpectralFunction {
      public:
      double gamman(double x, double q2);
      double sigma(double x, double q2);
        void DSh1int(double& DDDK0,double& DDDK,double& drho,double& DMU,double& spect);
       void DSh1int_proton(double& DDDK0,double& DDDK,double& drho,double& DMU,double& spect);
void DSh1int_neutron(double& DDDK0,double& DDDK,double& drho,double& DMU,double& spect);

        double findzero(double DK, double drho, double DMU);
        double Deriv(double DK0,double DK,double drho);
        double ff(double DK0,double DK,double drho);
        double DRG20R(double A, double B, int N,double CF[2000]);
      double reGS48(double X1,double X2,double CF[48]);
      void  STGS48(double& X1,double& X2,double X[48]);
      void DSG20R(double& A, double& B, int& N,  double X[2000], int& NP);
double startNucleus();
//void  startNucleus(double& rho0A,double& rho0p,double& rho0n);
double  startProton();
double  startNeutron();

double  Density(double );
double  Densityp(double r);
double  Densityn(double r);


//       private:
        void  ULIND(double& QZR,double& Q,double& XKF,std::complex<double>& CUFUN);
        double UIM(double W00,double q,double K1,double K2,double DM);
        double  DIMASM(double DK0,double DK);
        double  DFINM(double DQ,double DK0,double DK);
        void  SELA(double& E,double& SEL);
        double  reasi_int(double a,double b);
        void DIMAm1(double& DDK0,double& DDK,double& DMU,double& DIMA1);
};
class PionCloud{
      public:
// For Pions: helper function
 void forvlvt(double& p_zero, double& p_mom, double& Dzero, double& FormFactor, double&  DTILDEzero, double& FormFactorTILDE);

void forvlvtrho(double& p_zero, double& p_mom, double& Dzerorho, double& FormFactorRho, double&  DTILDEzerorho, double& FormFactorTILDErho);

void VlVtprima(double& p_zero, double& p_mom, double& Vlprima, double& Vtprima);
void deltaimD2(double& p_zero,double& p_mom,double& drho,double& deltaIm);
double F2pion(double x,double Q2);
double xseapion(double x,double Q2);
double xvalencepion(double x,double Q2);
private:

void ULIND2(double& QZR,double& q,double& drho,std::complex<double>& CUFUN);

void lindhard(double& q_zero,double& q_mod,double& drho,double& k_fermi,std::complex<double>& lind);
void delta_lind(double q_zero, double q_mod, double drho, double k_fermi,std::complex<double>& dlt_lind);
void deltaimDrho(double& p_zero,double& p_mom,double& drho,double& deltaImrho);
void deltaImDrho2(double& p_zero,double& p_mom,double& drho,double& deltaImrho);
void deltaimD(double& p_zero,double& p_mom,double& drho,double& deltaIm);
 void RCH1rho(int& m,double& h,double& drho,double& p_zero,double& p_mom,double D[100][100]);
 void RCH1(int& m,double& h,double& drho,double& p_zero,double& p_mom,double D[100][100]);
 void Difirstandrho(double& h,double& drho,double& p_zero,double& p_mom,double& D0,double& D0rho);
   std::complex<double> CDminusD0(double p_zero,double p_mom,double drho);
std::complex<double> CDminusD0rho(double p_zero,double p_mom,double drho);
};

class ShadowingEffect{
public:
void c2A_c3A_ct1A(double& x,double& q2, std::complex<double>& C2A,std::complex<double> & C3A, std::complex<double>& CT1A);
double deltaR2(double x,double aq2,std::string setname, int imem,int iset,const PDF* pdf);
double deltaR2_F2(double x,double aq2,std::string setname, int imem,int iset,const PDF* pdf);
double deltaR_F3(double x,double q2,std::string setname, int imem,int iset,const PDF* pdf);

};

#endif
