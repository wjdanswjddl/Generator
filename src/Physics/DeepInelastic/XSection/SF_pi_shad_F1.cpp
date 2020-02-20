#include "LHAPDF/LHAPDF.h"
#include <complex>
#include <algorithm>
#include "param.hh"
#include "NuclearEffect.hh"
#include "Nucleon.hh"

using namespace LHAPDF;
using namespace std;
using namespace Param;



// for Spectral Function which include Fermi motion, Pauli Blocking and Nucleon Binding

//This is the subroutine where Nucleon structure functions are called
void  D2(double& x, double& q2,double& p,double& r, double p0,const std::string setname, int& imem,int& iset,const PDF* pdf,double& F1facNuc,double& F1facProton,double& F1facNeutron);

// For  Isoscalar target

// calculates structue function convulated with spectral function. q2 in fm.
void spectralF1_Iso(double& x,double& q2,const std::string& setname, int& imem,int& iset,const PDF* pdf,double& resaf);

// For Non Isoscalar target
// calculates structue function convulated with spectral function. q2 in fm.

void spectralF1_proton(double& x,double& q2,const std::string& setname, int& imem,int& iset,const PDF* pdf,double& resafp);


void spectralF1_neutron(double& x,double& q2,const std::string& setname, int& imem,int& iset,const PDF* pdf,double& resafn);

// for pion q2 in fm
double F1ApionEM(double x,double q2);

//Shadowing Functions q2 in fm
double dshadowingF1(double x, double q2,std::string setname, int imem,int iset,const PDF* pdf);

void start_spectral_function() {
  SpectralFunction spectral;
  spectral.startNucleus();
  spectral.startProton();
  spectral.startNeutron();
}

// calculate W and check for given W limit
bool IsWaccept(double& x,double& q2, double& p0,double& p3,double& p,double& q0)
{
TargetMassEffect TMC;
double xw2fm=p0*p0-p*p-q2+2.0*q0*(p0-TMC.gamman(x,q2)*p3);
//cout<<xw2fm<<endl;

double wlimit_fm=(wlimit/GeVtfm)*(wlimit/GeVtfm); //in fm

          return (xw2fm > wlimit_fm);   }

// W for pion cloud
bool IswcutforPionCloud(double x,double q2, double p0,double p3,double p,double q0)
{
double gamma=1.0;
double  w2fm=p0*p0-p*p-q2-2.0*q0*(p0-gamma*p3);
//cout<<xw2fm<<endl;
double wlimit_fm=(wlimit/GeVtfm)*(wlimit/GeVtfm); //in fm
         return (w2fm > wlimit_fm);
}
// W for Shadowing
bool IsWacceptforShad(double& x,double& q2gev2)
{
   double w2ngev=Mn*Mn+q2gev2*(1.0/x-1.0);
       double wlimit2=wlimit*wlimit;
          return (w2ngev > wlimit2);
}


void D2(double& x, double& q2,double& p,double& r, double p0,const std::string setname, int& imem,int& iset,const PDF* pdf,double& F1facNuc,double& F1facProton,double& F1facNeutron)
{
SpectralFunction spectral;
NucleonNLO Nucleon;
TargetMassEffect TMC;
HigherTwist4 HTwist;
   double cosid[2500],Nucleonfcosid[2500],Protonfcosid[2500],Neutronfcosid[2500];
  double q2gev2=q2*GeVtfm*GeVtfm;

  int nt=1,iker;
  double cosmin=-1.0;
  double cosmax=1.0;
   double f2proton,f2neutron,f2Nucleon,f1Nucleon,f1Proton,f1Neutron;

  double cosi,p3,xx,q0,fact,f2P_TMC,f2N_TMC,f2NucTMC,f2htNucleon,f2htProton,f2htNeutron;
  double f1tmc_Nuc,f1tmc_P,f1tmc_N,f1htNucleon,f1htProton,f1htNeutron;
  int ntp;
  spectral.DSG20R(cosmin,cosmax,nt,cosid,ntp);
  for (int it=1;it<=ntp;++it)
    {
      cosi=cosid[it];
      p3=p*cosi;
      // xx is x_nucleon
      q0=q2/(2.0*DM*x);
      iker=1.0;
      xx=DM*x/(p0-TMC.gamman(x,q2)*p3);
      fact=(p*p-p3*p3)/(q0*(p0-TMC.gamman(x,q2)*p3)*2.00);
         if(noWCut) {goto g12;}

          else {if(IsWaccept(x,q2,p0,p3,p,q0)){
       g12:if(isNLO)
      {
       Nucleon.f2( xx,q2gev2,setname, imem,iset,pdf,f2proton,f2neutron,f2Nucleon);
        Nucleon.f1( xx,q2gev2,setname, imem,iset,pdf,f1Nucleon,f1Proton,f1Neutron);
               Nucleonfcosid[it]=f1Nucleon+fact*f2Nucleon;
      Protonfcosid[it]=f1Proton+fact*f2proton;
      Neutronfcosid[it]=f1Neutron+fact*f2neutron;

}
      else if(isTMC_NLO)
      {
      TMC.f1tmc(xx,q2gev2,setname,imem,iset,pdf,f1tmc_Nuc,f1tmc_P,f1tmc_N);
       TMC.f2tmc(xx,q2gev2,setname,imem,iset,pdf,f2P_TMC,f2N_TMC,f2NucTMC);
      Nucleonfcosid[it]=f1tmc_Nuc+fact*f2NucTMC;
      Protonfcosid[it]=f1tmc_P+fact*f2P_TMC;
      Neutronfcosid[it]=f1tmc_N+fact*f2N_TMC;
    }
      else if(isTMC_NLO_HT)
     {
       HTwist.f1ht(xx,q2gev2,setname,imem,iset,pdf,f1htNucleon,f1htProton,f1htNeutron);
       HTwist.f2ht(xx,q2gev2,setname,imem,iset,pdf,f2htNucleon,f2htProton,f2htNeutron);
       Nucleonfcosid[it]=f1htNucleon+fact*f2htNucleon;
      Protonfcosid[it]=f1htProton+fact*f2htProton;
      Neutronfcosid[it]=f1htNeutron+fact*f2htNeutron;

     }

}
else{
Nucleonfcosid[it]=0.0;
 Protonfcosid[it]=0.0;
Neutronfcosid[it]=0.0;
}
}
}
  F1facNuc=spectral.DRG20R(cosmin,cosmax,nt,Nucleonfcosid);
  F1facProton=spectral.DRG20R(cosmin,cosmax,nt,Protonfcosid);
  F1facNeutron=spectral.DRG20R(cosmin,cosmax,nt,Neutronfcosid);

//cout<<F1facNuc<<endl;

}





void spectralF1_Iso(double& x, double& q2, const std::string& setname, int& imem,int& iset,const PDF* pdf,double& resaf)
{
 SpectralFunction spectral;
  double specd[201][201][201][2];
  double spec1d[200][200][200];

  double  rd[2000],pd[2000],p0d[2000],rp0ad[2000],rpad[2000],rrad[2000];
     //int itime=0.0;

  double xm=DM;
   double rlim=12.0;
  double plimf=5.0;
  double p0lim1=-2.0;

  int nr=1;
  int np=1;
  int np0=1;
  int np0a=8;

  int nrp,npp,np0pa,np0p;
  double r,drho,fpf,ep0,p,p0lim2,DMU,silly,ep,p0,fker,dk0r;
  double resa1,resa2,resa3,resa,spect,res0a,F1facNuc, F1facProton, F1facNeutron;
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
      D2(x, q2, p, r, p0,setname,imem,iset,pdf,F1facNuc, F1facProton, F1facNeutron);
//      cout<<F1facNuc<<endl;
      rp0ad[ip0]=spect*F1facNuc;

      }
      res0a=spectral.DRG20R(p0lim1,p0lim2,np0a,rp0ad);
      rpad[ip]=res0a*p*p*fker;

    }
    resa1= spectral.DRG20R(pflow,pfup,np,rpad);
//     cout<<resa1<<endl;
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
       D2(x, q2, p, r, p0,setname,imem,iset,pdf,F1facNuc,F1facProton,F1facNeutron);
	rp0ad[ip0]=spect*F1facNuc;
      }
      res0a= spectral.DRG20R(p0lim1,p0up,np0,rp0ad);
      D2(x, q2, p, r, dk0r+xm,setname,imem,iset,pdf,F1facNuc,F1facProton,F1facNeutron);
      double d2fac=F1facNuc;
      rpad[ip]=res0a*p*p*fker;
      rpad[ip]=rpad[ip]+d2fac*p*p*fker/der;

    }
    resa2 = spectral.DRG20R(pfup,pf,np,rpad);
      // cout<<resa2<<endl;


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

      D2(x, q2, p, r, p0,setname,imem,iset,pdf,F1facNuc, F1facProton, F1facNeutron);
	rp0ad[ip0]=spect*F1facNuc;
      }

 //     cout<<spect<<" "<<F1facNuc<<endl;


      res0a= spectral.DRG20R(p0lim1,p0lim2,np0,rp0ad);
      rpad[ip]=res0a*p*p*fker;
    }
    resa3= spectral.DRG20R(pf,plimf,np,rpad);
    //    cout<<spect<<" "<<F1facNuc<<endl;

    resa=resa1+resa2+resa3;

    rrad[ir]=resa*r*r;
  }
   resaf= spectral.DRG20R(rlow,rlim,nr,rrad);


  double fac=(4.0*pi)*(2.0*pi)/((2.0*pi)*(2.0*pi)*(2.0*pi))*4.0/(double)A;


  resaf=resaf*fac*(double)A;
//cout<<resaf<<endl;
}

//proton

void spectralF1_proton(double& x,double& q2,const std::string& setname, int& imem,int& iset,const PDF* pdf,double& resafp)
{
SpectralFunction spectral;
  double specd[201][201][201][2];
  double spec1d[200][200][200];

  double  rd[2000],pd[2000],p0d[2000],rp0ad[2000],rpad[2000],rrad[2000];
     //int itime=0.0;

  double xm=DM;

  double rlim=12.0;
  double plimf=5.0;
  double p0lim1=-2.0;

  int nr=1;
  int np=1;
  int np0=1;
  int np0a=8;

  //       double x=xw;
  //       double q2=q2w;
  //       gamma=gammaw
  int nrp,npp,np0pa,np0p;
  double r,drho,fpf,ep0,p,p0lim2,DMU,silly,ep,p0,fker,dk0r;
  double resa1,resa2,resa3,resa,spect,res0a,F1facNuc,F1facProton,F1facNeutron;
;
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
      //        fker=1.0;
      fker=xm/ep;
      /////////// p0 loop
      spectral.DSG20R(p0lim1,p0lim2,np0a,p0d,np0pa);
      for(int ip0=1;ip0 <= np0pa;++ip0){
	ep0=p0d[ip0];
	p0=ep0+xm;
//	 if(itime != 137){
	spectral.DSh1int_proton(ep0,p,drho,DMU,spect);
	spec1d[ip0][ip][ir]=spect;//}
      D2(x, q2, p, r, p0,setname,imem,iset,pdf,F1facNuc,F1facProton,F1facNeutron);
	rp0ad[ip0]=spect*F1facProton;
}
      res0a=spectral.DRG20R(p0lim1,p0lim2,np0a,rp0ad);
      rpad[ip]=res0a*p*p*fker;

    }
    resa1= spectral.DRG20R(pflow,pfup,np,rpad);
    // second region


  g123:
     spectral.DSG20R(pfup,pf,np,pd,npp);
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
//	 if(itime != 137){
	spectral.DSh1int_proton(ep0,p,drho,DMU,spect);
	specd[ip0][ip][ir][0]=spect;//}
     D2(x, q2, p, r, p0,setname,imem,iset,pdf,F1facNuc,F1facProton,F1facNeutron);
	rp0ad[ip0]=spect*F1facProton;
      }
      res0a= spectral.DRG20R(p0lim1,p0up,np0,rp0ad);
      D2(x, q2, p, r, dk0r+xm,setname,imem,iset,pdf,F1facNuc,F1facProton,F1facNeutron);

      rpad[ip]=res0a*p*p*fker;
      rpad[ip]=rpad[ip]+F1facProton*p*p*fker/der;

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
	//	 if(itime != 137){
	spectral.DSh1int_proton(ep0,p,drho,DMU,spect);
	specd[ip0][ip][ir][1]=spect;//}

D2(x, q2, p, r, p0,setname,imem,iset,pdf,F1facNuc,F1facProton,F1facNeutron);

	rp0ad[ip0]=spect*F1facProton;
      }



      res0a= spectral.DRG20R(p0lim1,p0lim2,np0,rp0ad);
      rpad[ip]=res0a*p*p*fker;
      //         cout<<rpad[ip]<<endl;
    }
    resa3= spectral.DRG20R(pf,plimf,np,rpad);

    resa=resa1+resa2+resa3;

    rrad[ir]=resa*r*r;
  }
  double resaf1= spectral.DRG20R(rlow,rlim,nr,rrad);


  double fac=(4.0*pi)*(2.0*pi)/((2.0*pi)*(2.0*pi)*(2.0*pi))*2.0/(double)Z;


  resafp=resaf1*fac*(double)Z;
//         cout<<resaf1<<"  "<<fac<<"  "<<A<<endl;
  //c to distinguish first time
 //        itime=137;
}

/// Neutron

void spectralF1_neutron(double& x,double& q2,const std::string& setname, int& imem,int& iset,const PDF* pdf,double& resafn)
{
SpectralFunction spectral;
  double specd[201][201][201][2];
  double spec1d[200][200][200];

  double  rd[2000],pd[2000],p0d[2000],rp0ad[2000],rpad[2000],rrad[2000];
     //int itime=0.0;

  double xm=DM;

  double rlim=12.0;
  double plimf=5.0;
  double p0lim1=-2.0;

  int nr=1;
  int np=1;
  int np0=1;
  int np0a=8;

  //       double x=xw;
  //       double q2=q2w;
  //       gamma=gammaw
  int nrp,npp,np0pa,np0p;
  double r,drho,fpf,ep0,p,p0lim2,DMU,silly,ep,p0,fker,dk0r;
  double resa1,resa2,resa3,resa,spect,res0a,F1facNuc,F1facProton,F1facNeutron;

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
      //        fker=1.0;
      fker=xm/ep;
      /////////// p0 loop
      spectral.DSG20R(p0lim1,p0lim2,np0a,p0d,np0pa);
      for(int ip0=1;ip0 <= np0pa;++ip0){
	ep0=p0d[ip0];
	p0=ep0+xm;
//	 if(itime != 137){
	spectral.DSh1int_neutron(ep0,p,drho,DMU,spect);
	spec1d[ip0][ip][ir]=spect;//}
      D2(x, q2, p, r, p0,setname,imem,iset,pdf,F1facNuc,F1facProton,F1facNeutron);
	rp0ad[ip0]=spect*F1facNeutron;
}
      res0a=spectral.DRG20R(p0lim1,p0lim2,np0a,rp0ad);
      rpad[ip]=res0a*p*p*fker;

    }
    resa1= spectral.DRG20R(pflow,pfup,np,rpad);
    // second region


  g123:
     spectral.DSG20R(pfup,pf,np,pd,npp);
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
//	 if(itime != 137){
	spectral.DSh1int_neutron(ep0,p,drho,DMU,spect);
	specd[ip0][ip][ir][0]=spect;//}
     D2(x, q2, p, r, p0,setname,imem,iset,pdf,F1facNuc,F1facProton,F1facNeutron);
	rp0ad[ip0]=spect*F1facNeutron;
      }
      res0a= spectral.DRG20R(p0lim1,p0up,np0,rp0ad);
      D2(x, q2, p, r, dk0r+xm,setname,imem,iset,pdf,F1facNuc,F1facProton,F1facNeutron);

      rpad[ip]=res0a*p*p*fker;
      rpad[ip]=rpad[ip]+F1facNeutron*p*p*fker/der;

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
	//	 if(itime != 137){
	spectral.DSh1int_neutron(ep0,p,drho,DMU,spect);
	specd[ip0][ip][ir][1]=spect;//}

D2(x, q2, p, r, p0,setname,imem,iset,pdf,F1facNuc,F1facProton,F1facNeutron);

	rp0ad[ip0]=spect*F1facNeutron;
      }



      res0a= spectral.DRG20R(p0lim1,p0lim2,np0,rp0ad);
      rpad[ip]=res0a*p*p*fker;
      //         cout<<rpad[ip]<<endl;
    }
    resa3= spectral.DRG20R(pf,plimf,np,rpad);

    resa=resa1+resa2+resa3;

    rrad[ir]=resa*r*r;
  }
  double resaf1= spectral.DRG20R(rlow,rlim,nr,rrad);


  double fac=(4.0*pi)*(2.0*pi)/((2.0*pi)*(2.0*pi)*(2.0*pi))*2.0/(double)N;


  resafn=resaf1*fac*(double)N;
//         cout<<resaf1<<"  "<<fac<<"  "<<A<<endl;
  //c to distinguish first time
 //        itime=137;
}





//*******************************************************************
      double  F1ApionEM(double x,double aq2)
{
      SpectralFunction spectral;
      PionCloud pionEffect;


      double  rd[2000],pd[2000],dcosid[2000],p0d[2000];
      double rrd[2000],rpd[2000],rdcosid[2000],rp0d[2000];
//c We have an spectral for r from 0 to infinity, we must put a
//c cutoff in rlim
         double rlim=12.0;
// we have another spectral for p (momentum) from 0 to infinity,
//c and we must put another cut, let's take for instance plim=20.d0 and
//c let's see what happen.
         double plim=20.0;
         double dcosmin=0.0;
         double dcosmax=1.0;
//c         gamma=dsqrt(1.d0+4.d0*0.939d0**2*x**2/Q2gev2)
//c bjorken limit
         double gamma=1.0;
         double xxpi3=pow((2.0*pi),3);
         int nr=1;
         int np=4,npp,nrr,np00,ncosii;
         int ncosi=1;
         int np0=2;
         double p0,lim,p3,res,dcostheta,fact,aux;
         double Dzero,FormFactor, DTILDEzero,FormFactorTILDE,drho,deltaIm;
         double f1,fact1,xpion;

	   double q2=aq2/(GeVtfm*GeVtfm); // in fm
	   double q0=q2/(2.0*DM*x);
	   double rlow=0.0;
         double plow=0.0;
         spectral.DSG20R(rlow,rlim,nr,rd,nrr);
         spectral.DSG20R(plow,plim,np,pd,npp);

         int kk=0;
         for(int j=1;j<=nrr;++j){
            double r=rd[j];
            drho=spectral.Density(r);
            for(int k=1;k<=npp;k++){
               double p=pd[k];
               double p0limsup=(p*gamma)-(DM*x);
	       if(p0limsup > 0.0){
                  spectral.DSG20R(plow,p0limsup,np0,p0d,np00);
                  for(int m=1;m<=np00;++m){
                      p0=p0d[m];
                     pionEffect.deltaimD2(p0,p,drho,deltaIm);
                   lim=pow(10,-25);
		     if(abs(deltaIm)>lim){
                        kk=kk+1;
                        dcosmin=(DM*x+p0)/(gamma*p);
                        dcosmax=std::min(1.0,(DM+p0)/(gamma*p));
                        if(dcosmin < dcosmax){
                           spectral.DSG20R(dcosmin,dcosmax,ncosi,dcosid,ncosii);
                           for(int l=1;l<=ncosii;++l){
                              dcostheta=dcosid[l];
                              p3=p*dcostheta;
                              xpion=DM*x/(gamma*p3-p0);

                    fact=(p*p-p3*p3)/(q0*(gamma*p3-p0)*2.0);
         		f1=pionEffect.F2pion(xpion,aq2)/2.00/xpion;
		         fact1=f1+fact*pionEffect.F2pion(xpion,aq2);

                             if(noWCut){goto g12;}

          else {if(IswcutforPionCloud(xpion,q2,p0,p3,p,q0)){


                   g12:rdcosid[l]=-(6.0*4.0*pi*2.0*DM*fact1)/xxpi3;
} else{
           rdcosid[l]=0.0;
}

}
                           }
               res=spectral.DRG20R(dcosmin,dcosmax,ncosi,rdcosid);
}
                        else{ res=0.0; }

                        rp0d[m]=res*deltaIm;}
		     else{ rp0d[m]=0.0; }

                  }
                  res=spectral.DRG20R(plow,p0limsup,np0,rp0d);
//c new
//c substraction nucleon Lindhard
                  double resxxx=0.0;
                  double p0x=p*p/(2.0*DM);
                  if(p0x < p0limsup){
                     dcosmin=(DM*x+p0x)/(gamma*p);
                     dcosmax=min(1.0,(DM+p0x)/(gamma*p));

                     if(dcosmin < dcosmax){
                        spectral.DSG20R(dcosmin,dcosmax,ncosi,dcosid,ncosii);
                        for(int l=1;l<=ncosii;++l){
                           dcostheta=dcosid[l];
                           p3=p*dcostheta;
                             xpion=DM*x/(gamma*p3-p0x);

		  fact=(p*p-p3*p3)/(q0*(gamma*p3-p0x)*2.0);
			 f1=pionEffect.F2pion(xpion,aq2)/2.00/xpion;
			 fact1=f1+fact*pionEffect.F2pion(xpion,aq2);
                      if(noWCut){goto g13;}

          else {if(IswcutforPionCloud(xpion,q2,p0x,p3,p,q0)){


                           g13:rdcosid[l]=-6.0*4.0*pi/xxpi3*2.0*DM*fact1;
} else{rdcosid[l]=0.0;}
}

                        }
                        resxxx= spectral.DRG20R(dcosmin,dcosmax,ncosi,rdcosid);
}
                     else{resxxx=0.0; }
                  }
      pionEffect.forvlvt(p0x, p, Dzero, FormFactor, DTILDEzero,FormFactorTILDE);

              aux=Dzero*Dzero*(f2/(xmpion*xmpion))*FormFactor*FormFactor*p*p;

                  res=res+resxxx*pi*drho*aux;
		  }
               else{res=0.0;}
               rpd[k]=res*p*p;
            }
            res= spectral.DRG20R(plow,plim,np,rpd);
            rrd[j]=res*r*r;
         }
         res= spectral.DRG20R(rlow,rlim,nr,rrd);

         return res/A;
        }

//Shadowing Effect
 //    Q2 is given in fm**(-2)

double dshadowingF1(double x, double q2,std::string setname, int imem,int iset,const PDF* pdf)
{
ShadowingEffect shadowEff;
ForShadowinfF1 shadowf1;
     double shadow;
     double q2gev2=q2*GeVtfm*GeVtfm;
	 //double qgev=sqrt(q2gev2);
      double deltR=shadowEff.deltaR2(x,q2gev2,setname,imem,iset,pdf);
      double f1free= shadowf1.f1nucshad( x, q2gev2, setname, imem,iset,pdf);
      if(noWCut) {goto g12;}

          else{if(IsWacceptforShad(x,q2gev2)){

   	g12:shadow=f1free*deltR;
}else{shadow=0.0;}
}
	return shadow;
}



