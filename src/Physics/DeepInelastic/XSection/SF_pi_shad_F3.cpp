#include "LHAPDF/LHAPDF.h"
#include <complex>
#include <algorithm>
#include "param.hh"
#include "NuclearEffect.hh"
#include "Nucleon.hh"
using namespace LHAPDF;
using namespace std;
using namespace Param;

// Declare helper functions defined in SF_pi_shad_F1.cpp
bool IswcutforPionCloud(double x,double q2, double p0,double p3,double p,double q0);
bool IsWacceptforShad(double& x,double& q2gev2);
void  D2(double& x, double& q2,double& p,double& r, double p0,const std::string setname, int& imem,int& iset,const PDF* pdf,double& F2facNuc,double& F2facProton,double& F2facNeutron);

// Nuclear F3
// for Spectral Function which include Fermi motion, Pauli Blocking and Nucleon Binding

//This is the subroutine where Nucleon structure functions are called


void D2(double& x, double& q2,double& p,double& r, double p0,std::string setname, int& imem,int& iset,const PDF* pdf,double& F3FacNuc,double& F3FacP,double& F3FacN);

// For  Isoscalar target

// calculates structue function convulated with spectral function. q2 in fm.

void spectralF3_Iso(double& x,double& q2,const std::string& setname, int& imem,int& iset,const PDF* pdf,double& resaf);

// For Non Isoscalar target
// calculates structue function convulated with spectral function. q2 in fm.


void spectralF3_proton(double& x,double& q2,const std::string& setname, int& imem,int& iset,const PDF* pdf,double& resafp);

void spectralF3_neutron(double& x,double& q2,const std::string& setname, int& imem,int& iset,const PDF* pdf,double& resafn);


//Shadowing Functions q2 in fm
double dshadowingf3(double x, double q2,std::string setname, int imem,int iset,const PDF* pdf);

bool IsWaccept(double x,double& q2, double& p0,double& p3,double& p,double& q0)
{
TargetMassEffect TMC;
double xw2fm=p0*p0-p*p-q2+2.0*q0*(p0-TMC.gamman(x,q2)*p3);
//cout<<xw2fm<<endl;

double wlimit_fm=(wlimit/GeVtfm)*(wlimit/GeVtfm); //in fm

          return (xw2fm > wlimit_fm);
}

void spectralF3_Iso(double& x,double& q2,const std::string& setname, int& imem,int& iset,const PDF* pdf,double& resaf)
{
 SpectralFunction spectral;
  double specd[201][201][201][2];
  double spec1d[200][200][200];

  double  rd[2000],pd[2000],p0d[2000],rp0ad[2000],rpad[2000],rrad[2000];
     //int itime=0.0;

  double xm=DM;
  // bool isNLO,isTMC_NLO,isTMC_NLO_HT;
  double rlim=12.0;
  double plimf=5.0;
  double p0lim1=-2.0;

  int nr=1;
  int np=1;
  int np0=1;
  int np0a=8;

  int nrp,npp,np0pa,np0p;
  double r,drho,fpf,ep0,p,p0lim2,DMU,silly,ep,p0,fker,dk0r;
  double resa1,resa2,resa3,resa,spect,res0a,F3FacNuc,F3FacP,F3FacN;
  double rlow=0.0;
  double pflow=0.0;
  double pfup,der,dif,p0up;

  ///////// r loop
  spectral.DSG20R(rlow,rlim,nr,rd,nrp);
  for(int ir=1;ir<=nrp;++ir){
    r=rd[ir];
    drho=spectral.Density(r);
    double pf=pow((3.0/2.0*pi*pi*drho),(1.0/3.0));

    if(drho > 0.10){fpf=0.8;}
    else if(drho > 0.05){fpf=0.6;}
    else if(drho > 0.020){fpf=0.4;}
    else {fpf=0.0;}
    ep0=0.0;
    p=0.0;
    spectral.DSh1int(ep0,p,drho,DMU,silly);
    p0lim2=DMU;
    resa1=0.0;

    pfup=fpf*pf;

    if(fpf < 0.10){goto g123;}

    ///////////////////p loop
    spectral.DSG20R(pflow,pfup,np,pd,npp);
    for(int ip=1;ip<=npp;++ip){
      p=pd[ip];
      ep=sqrt(p*p+xm*xm);
      fker=xm/ep;
      /////////// p0 loop
      spectral.DSG20R(p0lim1,p0lim2,np0a,p0d,np0pa);
      for(int ip0=1;ip0 <= np0pa;++ip0){
	ep0=p0d[ip0];
	p0=ep0+xm;
	spectral.DSh1int(ep0,p,drho,DMU,spect);
	spec1d[ip0][ip][ir]=spect;
      D2(x, q2, p, r, p0,setname,imem,iset,pdf,F3FacNuc,F3FacP,F3FacN);

      rp0ad[ip0]=spect*F3FacNuc;

      }
      res0a=spectral.DRG20R(p0lim1,p0lim2,np0a,rp0ad);
      rpad[ip]=res0a*p*p*fker;

    }
    resa1= spectral.DRG20R(pflow,pfup,np,rpad);
    // second region


  g123:spectral.DSG20R(pfup,pf,np,pd,npp);
    for(int ip=1;ip<=npp;++ip){
      p=pd[ip];
      ep=sqrt(p*p+xm*xm);
      fker=xm/ep;
      dk0r=spectral.findzero(p,drho,DMU);
      der=spectral.Deriv(dk0r,p,drho);
      dif=DMU-dk0r;
      p0up=dk0r-dif;
      spectral.DSG20R(p0lim1,p0up,np0,p0d,np0p);
      for(int ip0=1;ip0<=np0p;++ip0){
	ep0=p0d[ip0];
	p0=ep0+xm;
	spectral.DSh1int(ep0,p,drho,DMU,spect);
	specd[ip0][ip][ir][0]=spect;
       D2(x, q2, p, r, p0,setname,imem,iset,pdf,F3FacNuc,F3FacP,F3FacN);
	rp0ad[ip0]=spect*F3FacNuc;
      }
      res0a= spectral.DRG20R(p0lim1,p0up,np0,rp0ad);
      D2(x, q2, p, r, dk0r+xm,setname,imem,iset,pdf,F3FacNuc,F3FacP,F3FacN);
      double d2fac=F3FacNuc;
      rpad[ip]=res0a*p*p*fker;
      rpad[ip]=rpad[ip]+d2fac*p*p*fker/der;

    }
    resa2 = spectral.DRG20R(pfup,pf,np,rpad);


    //c 3rd region
    spectral.DSG20R(pf,plimf,np,pd,npp);
    for(int ip=1;ip<=npp;++ip){
      p=pd[ip];
      ep=sqrt(p*p+xm*xm);
      fker=xm/ep;
      spectral.DSG20R(p0lim1,p0lim2,np0,p0d,np0p);
      for(int ip0=1;ip0 <= np0p;++ip0){
	ep0=p0d[ip0];
	p0=ep0+xm;
	spectral.DSh1int(ep0,p,drho,DMU,spect);
	specd[ip0][ip][ir][1]=spect;

      D2(x, q2, p, r, p0,setname,imem,iset,pdf,F3FacNuc,F3FacP,F3FacN);
	rp0ad[ip0]=spect*F3FacNuc;
      }



      res0a= spectral.DRG20R(p0lim1,p0lim2,np0,rp0ad);
      rpad[ip]=res0a*p*p*fker;
    }
    resa3= spectral.DRG20R(pf,plimf,np,rpad);

    resa=resa1+resa2+resa3;

    rrad[ir]=resa*r*r;
  }
  double resaf1= spectral.DRG20R(rlow,rlim,nr,rrad);


  double fac=(4.0*pi)*(2.0*pi)/((2.0*pi)*(2.0*pi)*(2.0*pi))*4.0/(double)A;


  resaf=resaf1*fac*(double)A;
}

void spectralF3_proton(double& x,double& q2,const std::string& setname, int& imem,int& iset,const PDF* pdf,double& resafp)
{
 SpectralFunction spectral;
  double specd[201][201][201][2];
  double spec1d[200][200][200];

  double  rd[2000],pd[2000],p0d[2000],rp0ad[2000],rpad[2000],rrad[2000];
     //int itime=0.0;

  double xm=DM;
  // bool isNLO,isTMC_NLO,isTMC_NLO_HT;
  double rlim=2.5*aproton;
  double plimf=5.0;
  double p0lim1=-2.0;

  int nr=1;
  int np=1;
  int np0=1;
  int np0a=8;

  int nrp,npp,np0pa,np0p;
  double r,drho,fpf,ep0,p,p0lim2,DMU,silly,ep,p0,fker,dk0r;
  double resa1,resa2,resa3,resa,spect,res0a,F3FacNuc,F3FacP,F3FacN;
  double rlow=0.0;
  double pflow=0.0;
  double pfup,der,dif,p0up;

  ///////// r loop
  spectral.DSG20R(rlow,rlim,nr,rd,nrp);
  for(int ir=1;ir<=nrp;++ir){
    r=rd[ir];
    drho=spectral.Densityp(r);
    double pf=pow((3.0*pi*pi*drho),(1.0/3.0));

    if(drho > 0.10){fpf=0.8;}
    else if(drho > 0.05){fpf=0.6;}
    else if(drho > 0.020){fpf=0.4;}
    else {fpf=0.0;}
    ep0=0.0;
    p=0.0;
    spectral.DSh1int_proton(ep0,p,drho,DMU,silly);
    p0lim2=DMU;
    resa1=0.0;

    pfup=fpf*pf;

    if(fpf < 0.10){goto g123;}

    ///////////////////p loop
    spectral.DSG20R(pflow,pfup,np,pd,npp);
    for(int ip=1;ip<=npp;++ip){
      p=pd[ip];
      ep=sqrt(p*p+xm*xm);
      fker=xm/ep;
      /////////// p0 loop
      spectral.DSG20R(p0lim1,p0lim2,np0a,p0d,np0pa);
      for(int ip0=1;ip0 <= np0pa;++ip0){
	ep0=p0d[ip0];
	p0=ep0+xm;
	spectral.DSh1int_proton(ep0,p,drho,DMU,spect);
	spec1d[ip0][ip][ir]=spect;
      D2(x, q2, p, r, p0,setname,imem,iset,pdf,F3FacNuc,F3FacP,F3FacN);

      rp0ad[ip0]=spect*F3FacP;

      }
      res0a=spectral.DRG20R(p0lim1,p0lim2,np0a,rp0ad);
      rpad[ip]=res0a*p*p*fker;

    }
    resa1= spectral.DRG20R(pflow,pfup,np,rpad);
    // second region


  g123:spectral.DSG20R(pfup,pf,np,pd,npp);
    for(int ip=1;ip<=npp;++ip){
      p=pd[ip];
      ep=sqrt(p*p+xm*xm);
      fker=xm/ep;
      dk0r=spectral.findzero(p,drho,DMU);
      der=spectral.Deriv(dk0r,p,drho);
      dif=DMU-dk0r;
      p0up=dk0r-dif;
      spectral.DSG20R(p0lim1,p0up,np0,p0d,np0p);
      for(int ip0=1;ip0<=np0p;++ip0){
	ep0=p0d[ip0];
	p0=ep0+xm;
	spectral.DSh1int_proton(ep0,p,drho,DMU,spect);
	specd[ip0][ip][ir][0]=spect;
       D2(x, q2, p, r, p0,setname,imem,iset,pdf,F3FacNuc,F3FacP,F3FacN);
	rp0ad[ip0]=spect*F3FacP;
      }
      res0a= spectral.DRG20R(p0lim1,p0up,np0,rp0ad);
      D2(x, q2, p, r, dk0r+xm,setname,imem,iset,pdf,F3FacNuc,F3FacP,F3FacN);
      double d2fac=F3FacP;
      rpad[ip]=res0a*p*p*fker;
      rpad[ip]=rpad[ip]+d2fac*p*p*fker/der;

    }
    resa2 = spectral.DRG20R(pfup,pf,np,rpad);


    //c 3rd region
    spectral.DSG20R(pf,plimf,np,pd,npp);
    for(int ip=1;ip<=npp;++ip){
      p=pd[ip];
      ep=sqrt(p*p+xm*xm);
      fker=xm/ep;
      spectral.DSG20R(p0lim1,p0lim2,np0,p0d,np0p);
      for(int ip0=1;ip0 <= np0p;++ip0){
	ep0=p0d[ip0];
	p0=ep0+xm;
	spectral.DSh1int_proton(ep0,p,drho,DMU,spect);
	specd[ip0][ip][ir][1]=spect;

      D2(x, q2, p, r, p0,setname,imem,iset,pdf,F3FacNuc,F3FacP,F3FacN);
	rp0ad[ip0]=spect*F3FacP;
      }



      res0a= spectral.DRG20R(p0lim1,p0lim2,np0,rp0ad);
      rpad[ip]=res0a*p*p*fker;
    }
    resa3= spectral.DRG20R(pf,plimf,np,rpad);

    resa=resa1+resa2+resa3;

    rrad[ir]=resa*r*r;
  }
  double resaf1= spectral.DRG20R(rlow,rlim,nr,rrad);


  double fac=(4.0*pi)*(2.0*pi)/((2.0*pi)*(2.0*pi)*(2.0*pi))*2.0/(double)Z;


  resafp=resaf1*fac*(double)Z;
}

void spectralF3_neutron(double& x,double& q2,const std::string& setname, int& imem,int& iset,const PDF* pdf,double& resafn)
{
 SpectralFunction spectral;
  double specd[201][201][201][2];
  double spec1d[200][200][200];

  double  rd[2000],pd[2000],p0d[2000],rp0ad[2000],rpad[2000],rrad[2000];
     //int itime=0.0;

  double xm=DM;
  // bool isNLO,isTMC_NLO,isTMC_NLO_HT;
  double rlim=2.5*aneutron;
  double plimf=5.0;
  double p0lim1=-2.0;

  int nr=1;
  int np=1;
  int np0=1;
  int np0a=8;

  int nrp,npp,np0pa,np0p;
  double r,drho,fpf,ep0,p,p0lim2,DMU,silly,ep,p0,fker,dk0r;
  double resa1,resa2,resa3,resa,spect,res0a,F3FacNuc,F3FacP,F3FacN;
  double rlow=0.0;
  double pflow=0.0;
  double pfup,der,dif,p0up;

  ///////// r loop
  spectral.DSG20R(rlow,rlim,nr,rd,nrp);
  for(int ir=1;ir<=nrp;++ir){
    r=rd[ir];
    drho=spectral.Densityn(r);
    double pf=pow((3.0*pi*pi*drho),(1.0/3.0));

    if(drho > 0.10){fpf=0.8;}
    else if(drho > 0.05){fpf=0.6;}
    else if(drho > 0.020){fpf=0.4;}
    else {fpf=0.0;}
    ep0=0.0;
    p=0.0;
    spectral.DSh1int_neutron(ep0,p,drho,DMU,silly);
    p0lim2=DMU;
    resa1=0.0;

    pfup=fpf*pf;

    if(fpf < 0.10){goto g123;}

    ///////////////////p loop
    spectral.DSG20R(pflow,pfup,np,pd,npp);
    for(int ip=1;ip<=npp;++ip){
      p=pd[ip];
      ep=sqrt(p*p+xm*xm);
      fker=xm/ep;
      /////////// p0 loop
      spectral.DSG20R(p0lim1,p0lim2,np0a,p0d,np0pa);
      for(int ip0=1;ip0 <= np0pa;++ip0){
	ep0=p0d[ip0];
	p0=ep0+xm;
	spectral.DSh1int_neutron(ep0,p,drho,DMU,spect);
	spec1d[ip0][ip][ir]=spect;
      D2(x, q2, p, r, p0,setname,imem,iset,pdf,F3FacNuc,F3FacP,F3FacN);

      rp0ad[ip0]=spect*F3FacN;

      }
      res0a=spectral.DRG20R(p0lim1,p0lim2,np0a,rp0ad);
      rpad[ip]=res0a*p*p*fker;

    }
    resa1= spectral.DRG20R(pflow,pfup,np,rpad);
    // second region


  g123:spectral.DSG20R(pfup,pf,np,pd,npp);
    for(int ip=1;ip<=npp;++ip){
      p=pd[ip];
      ep=sqrt(p*p+xm*xm);
      fker=xm/ep;
      dk0r=spectral.findzero(p,drho,DMU);
      der=spectral.Deriv(dk0r,p,drho);
      dif=DMU-dk0r;
      p0up=dk0r-dif;
      spectral.DSG20R(p0lim1,p0up,np0,p0d,np0p);
      for(int ip0=1;ip0<=np0p;++ip0){
	ep0=p0d[ip0];
	p0=ep0+xm;
	spectral.DSh1int_neutron(ep0,p,drho,DMU,spect);
	specd[ip0][ip][ir][0]=spect;
       D2(x, q2, p, r, p0,setname,imem,iset,pdf,F3FacNuc,F3FacP,F3FacN);
	rp0ad[ip0]=spect*F3FacN;
      }
      res0a= spectral.DRG20R(p0lim1,p0up,np0,rp0ad);
      D2(x, q2, p, r, dk0r+xm,setname,imem,iset,pdf,F3FacNuc,F3FacP,F3FacN);
      double d2fac=F3FacN;
      rpad[ip]=res0a*p*p*fker;
      rpad[ip]=rpad[ip]+d2fac*p*p*fker/der;

    }
    resa2 = spectral.DRG20R(pfup,pf,np,rpad);


    //c 3rd region
    spectral.DSG20R(pf,plimf,np,pd,npp);
    for(int ip=1;ip<=npp;++ip){
      p=pd[ip];
      ep=sqrt(p*p+xm*xm);
      fker=xm/ep;
      spectral.DSG20R(p0lim1,p0lim2,np0,p0d,np0p);
      for(int ip0=1;ip0 <= np0p;++ip0){
	ep0=p0d[ip0];
	p0=ep0+xm;
	spectral.DSh1int_neutron(ep0,p,drho,DMU,spect);
	specd[ip0][ip][ir][1]=spect;

      D2(x, q2, p, r, p0,setname,imem,iset,pdf,F3FacNuc,F3FacP,F3FacN);
	rp0ad[ip0]=spect*F3FacN;
      }



      res0a= spectral.DRG20R(p0lim1,p0lim2,np0,rp0ad);
      rpad[ip]=res0a*p*p*fker;
    }
    resa3= spectral.DRG20R(pf,plimf,np,rpad);

    resa=resa1+resa2+resa3;

    rrad[ir]=resa*r*r;
  }
  double resaf1= spectral.DRG20R(rlow,rlim,nr,rrad);


  double fac=(4.0*pi)*(2.0*pi)/((2.0*pi)*(2.0*pi)*(2.0*pi))*2.0/(double)N;


  resafn=resaf1*fac*(double)N;
}





//Shadowing Effect
 //    Q2 is given in fm**(-2)



//////////////// For Shadowing F3

//c Computing of shadowing effect for F3 following the model by
//c Kulagin and Petti, Nuclear Physics A 765, 126 (2006)
//c and Phys.Rev.D76:094023,2007

	double dshadowingf3(double x, double q2,std::string setname, int imem,int iset,const PDF* pdf)
{
      double q2gev2=q2*GeVtfm*GeVtfm;
	//double qgev=sqrt(q2gev2);
      double f3Nucleon,f3Proton,f3Neutron,shadow,deltR3=0.0;
      NucleonNLO Nucleon;
      ShadowingEffect shadF3;
        Nucleon.f3(x,q2gev2,setname, imem,iset,pdf,f3Nucleon,f3Proton,f3Neutron);
  if(noWCut) {goto g12;}

          else{if(IsWacceptforShad(x,q2gev2)){

          deltR3=shadF3.deltaR_F3(x,q2gev2,setname,imem,iset,pdf);
              g12:shadow= f3Nucleon*deltR3; }
else
   {shadow=0.0;}
 }
	return shadow;
}




