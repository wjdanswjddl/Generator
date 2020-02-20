#ifndef _Nucleon_H_
#define _Nucleon_H_
#include "LHAPDF/LHAPDF.h"
#include <iostream>
#include <string>
#include <cmath>
#include <complex>
#include <algorithm>
#include "param.hh"

using namespace LHAPDF;
using namespace Param;
class NucleonNLO{
      public:
     void f3(double x, double& aq2,const std::string setname, int& imem,int& iset,const PDF* pdf,double& f3Nucleon,double& f3Proton,double& f3Neutron);


void f2(double x, double& q,const std::string& setname, int& imem,int& iset,const PDF* pdf,double& f2proton, double& f2neutron, double& f2nucleon);

void f1(double x,double& aq2,const std::string setname, int& imem,int& iset,const PDF* pdf,double& f1Nucleon,double& f1Proton,double& f1Neutron);

void Fl(double& x, double& aq2,const std::string setname, int& imem,int& iset,const PDF* pdf,double& flNucleon,double& flProton,double& flNeutron);

void f3pdf(double& x,double& aq2,double& as,double& f3pdfx,double& f3pdfxP,double& f3pdfxN,const std::string& setname, int& imem,int& iset,const PDF* pdf);
void xf2nuc(double& x,double& q,double& as,double& xsinglet,double& xnonsinglet,double& xgluon,
	    const std::string& setname, int& imem,int& iset,const PDF* pdf);
void xf2proton(double& x,double& q,double& as,double& xsingproton,double& xnonsingproton,double& xgluon_P,
	    const std::string& setname, int& imem,int& iset,const PDF* pdf);
void xf2neutron(double& x,double& q,double& as,double& xsingneutron,double& xnonsingneutron,double& xgluon_N,
	    const std::string& setname, int& imem,int& iset,const PDF* pdf);
};
class TargetMassEffect{
      public:

double gamman(double x, double q2);
double sigma(double x, double q2);

void f1tmc(double x, double& q2gev2,const std::string setname, int& imem,int& iset,const PDF* pdf,double& f1tmc_Nuc,double& f1tmc_P,double& f1tmc_N);
void f2tmc(double& x,double& q2gev2,const std::string setname, int& imem,int& iset,const PDF* pdf,double& f2P_TMC,double& f2N_TMC,double& f2NucTMC);
void f3tmc(double x, double& q2gev2,const std::string setname, int& imem,int& iset,const PDF* pdf,double& f3tmc_Nuc,double& f3tmc_P,double& f3tmc_N);

};

class HigherTwist4{
public:
void pdfs_nucleon(double& x,double& q,double& as,double& xsinglet,double& xnonsinglet,double& xgluon, const std::string& setname, int& imem,int& iset,const PDF* pdf);

void pdfs_proton(double& x,double& aq2,double& as,double& xsingproton,double& xnonsingproton,double& xgluon_P, const std::string& setname, int& imem,int& iset,const PDF* pdf);

void pdfs_neutron(double& x,double& aq2,double& as,double& xsingneutron,double& xnonsingneutron,double& xgluon_N, const std::string& setname, int& imem,int& iset,const PDF* pdf);

void htns_nucleon(double& x, double& aq2,const std::string setname, int& imem,int& iset,const PDF* pdf,double& htNucleon,double& htProton,double& htNeutron);
void f2ht(double& x, double& aq2,const std::string setname, int& imem,int& iset,const PDF* pdf,double& f2htNucleon, double& f2htProton, double& f2htNeutron);

// for F1 HT
void htns_f1(double x, double aq2,std::string setname, int imem,int iset,const PDF* pdf,double& h1tNucleon,double& h1tProton,double& h1tNeutron);
void f1ht(double& x, double& aq2,const std::string setname, int& imem,int& iset,const PDF* pdf,double& f1htNucleon, double& f1htProton, double& f1htNeutron);
////// For F3 HT
void f3ht(double& x, double& aq2,const std::string setname, int& imem,int& iset,const PDF* pdf,double& f3htNucleon, double& f3htProton, double& f3htNeutron);
};
class ForShadowinfF1{
public:
double  Flongitudinal(double x, double qgev,std::string setname, int imem,int iset,const PDF* pdf);
double f1nucshad(double x,double q,std::string setname, int imem,int iset,const PDF* pdf);
};
class Integration{
public:
  double DRG20R(double A, double B, int N,double CF[2000]);
      double reGS48(double X1,double X2,double CF[48]);
      void  STGS48(double& X1,double& X2,double X[48]);
      void DSG20R(double& A, double& B, int& N,  double X[2000], int& NP);
};
#endif
