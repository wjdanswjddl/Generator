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

// For  Isoscalar target

// calculates structue function convulated with spectral function. q2 in fm.
void spectralF2_Iso(double& x,double& q2,const std::string& setname, int& imem,int& iset,const PDF* pdf,double& resaf);


// For Non Isoscalar target
// calculates structue function convulated with spectral function. q2 in fm.

void spectralF2_proton(double& x,double& q2,const std::string& setname, int& imem,int& iset,const PDF* pdf,double& resafp);


void spectralF2_neutron(double& x,double& q2,const std::string& setname, int& imem,int& iset,const PDF* pdf,double& resafn);

// for pion q2 in fm
double F2ApionEM(double x,double q2);

//Shadowing Functions q2 in fm
double dshadowingF2(double x, double q2,std::string setname, int imem,int iset,const PDF* pdf);

void spectralF2_Iso(double& x,double& q2,const std::string& setname, int& imem,int& iset,const PDF* pdf,double& resaf)
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
  double resa1,resa2,resa3,resa,spect,res0a,F2facNuc, F2facProton, F2facNeutron;
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
      D2(x, q2, p, r, p0,setname,imem,iset,pdf,F2facNuc, F2facProton, F2facNeutron);

      rp0ad[ip0]=spect*F2facNuc;

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
       D2(x, q2, p, r, p0,setname,imem,iset,pdf,F2facNuc,F2facProton,F2facNeutron);
	rp0ad[ip0]=spect*F2facNuc;
      }
      res0a= spectral.DRG20R(p0lim1,p0up,np0,rp0ad);
      D2(x, q2, p, r, dk0r+xm,setname,imem,iset,pdf,F2facNuc,F2facProton,F2facNeutron);
      double d2fac=F2facNuc;
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

      D2(x, q2, p, r, p0,setname,imem,iset,pdf,F2facNuc, F2facProton, F2facNeutron);
	rp0ad[ip0]=spect*F2facNuc;
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


void spectralF2_proton(double& x,double& q2,const std::string& setname, int& imem,int& iset,const PDF* pdf,double& resafp)
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
  double resa1,resa2,resa3,resa,spect,res0a,F2facNuc, F2facProton, F2facNeutron;
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
      D2(x, q2, p, r, p0,setname,imem,iset,pdf,F2facNuc, F2facProton, F2facNeutron);

      rp0ad[ip0]=spect*F2facProton;

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
       D2(x, q2, p, r, p0,setname,imem,iset,pdf,F2facNuc,F2facProton,F2facNeutron);
	rp0ad[ip0]=spect*F2facProton;
      }
      res0a= spectral.DRG20R(p0lim1,p0up,np0,rp0ad);
      D2(x, q2, p, r, dk0r+xm,setname,imem,iset,pdf,F2facNuc,F2facProton,F2facNeutron);
      double d2fac=F2facProton;
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

      D2(x, q2, p, r, p0,setname,imem,iset,pdf,F2facNuc, F2facProton, F2facNeutron);
	rp0ad[ip0]=spect*F2facProton;
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


void spectralF2_neutron(double& x,double& q2,const std::string& setname, int& imem,int& iset,const PDF* pdf,double& resafn)
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
  double resa1,resa2,resa3,resa,spect,res0a,F2facNuc, F2facProton, F2facNeutron;
  double rlow=0.0;
  double pflow=0.0;
  double pfup,der,dif,p0up;

  ///////// r loop
  spectral.DSG20R(rlow,rlim,nr,rd,nrp);
  for(int ir=1;ir<=nrp;++ir){
    r=rd[ir];
    drho=spectral.Densityn(r);
//    cout<<drho<<endl;
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
      D2(x, q2, p, r, p0,setname,imem,iset,pdf,F2facNuc, F2facProton, F2facNeutron);

      rp0ad[ip0]=spect*F2facNeutron;

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
       D2(x, q2, p, r, p0,setname,imem,iset,pdf,F2facNuc,F2facProton,F2facNeutron);
	rp0ad[ip0]=spect*F2facNeutron;
      }
      res0a= spectral.DRG20R(p0lim1,p0up,np0,rp0ad);
      D2(x, q2, p, r, dk0r+xm,setname,imem,iset,pdf,F2facNuc,F2facProton,F2facNeutron);
      double d2fac=F2facNeutron;
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

      D2(x, q2, p, r, p0,setname,imem,iset,pdf,F2facNuc, F2facProton, F2facNeutron);
	rp0ad[ip0]=spect*F2facNeutron;
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



//*******************************************************************

//aq2 in GeV
      double  F2ApionEM(double x,double aq2)
{

 SpectralFunction spectral;
 PionCloud pionEffect;
      double  rd[2000],pd[2000],dcosid[2000],p0d[2000];
      double rrd[2000],rpd[2000],rdcosid[2000],rp0d[2000];

//c We have an integral for r from 0 to infinity, we must put a
//c cutoff in rlim
         double rlim=12.0;
// we have another integral for p (momentum) from 0 to infinity,
//c and we must put another cut, let's take for instance plim=20.d0 and
//c let's see what happen.
         double plim=20.0;
         double dcosmin=0.0;
         double dcosmax=1.0;
//c bjorken limit
         double gamma=1.0;

         int nr=1;
         int np=4,npp,nrr,np00,ncosii;
         int ncosi=1;
         int np0=2;
         double p0,lim,p3,res,dcostheta,fact,aux,fac0,apq,fac2,p0x;
         double Dzero,FormFactor, DTILDEzero,FormFactorTILDE,drho,deltaIm;
         double xpion;
//        double Q2gev2=q2/(GeVtofm*GeVtofm);

	   double q2=aq2/(GeVtfm*GeVtfm); // in fm
	   double q0=q2/(2.0*DM*x);
	   double rlow=0.0;
         double plow=0.0;
         double xxpi3=pow((2.0*pi),3);

         spectral.DSG20R(rlow,rlim,nr,rd,nrr);
//c            p0limsup=plim
         spectral.DSG20R(plow,plim,np,pd,npp);

         int kk=0;
         for(int j=1;j<=nrr;++j){
            double r=rd[j];
            drho=spectral.Density(r);
            for(int k=1;k<=npp;k++){
               double p=pd[k];
               double p0limsup=(p*gamma)-(DM*x);
//            cout<<"p0limsup= "<<p0limsup<<endl;
	       if(p0limsup > 0.0){
                  spectral.DSG20R(plow,p0limsup,np0,p0d,np00);
                  for(int m=1;m<=np00;++m){
                      p0=p0d[m];
                     pionEffect.deltaimD2(p0,p,drho,deltaIm);
//                   cout<<"deltaIm= "<<deltaIm<<endl;
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

         fac0=q2/((gamma*q0*gamma*q0))*(p*p-p3*p3)/(2.0*mpi*mpi);
	 apq=p3*gamma-p0;
	fac2=(p3*q2/(gamma*q0*q0*apq)+1.0)*(p3*q2/(gamma*q0*q0*apq)+1.00);
	 fact=(fac0+(apq*apq*fac2)/(mpi*mpi))*(mpi/(p3*gamma-p0));



                         if(noWCut){goto g14;}
    else {if(IswcutforPionCloud(xpion,q2,p0,p3,p,q0)){

       g14:rdcosid[l]=-(6.0*4.0*pi*2.0*mpi*fact*pionEffect.F2pion(xpion,aq2))/xxpi3;
}
else{rdcosid[l]=0.0;}
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
                  p0x=p*p/(2.0*DM);
                  if(p0x < p0limsup){
                     dcosmin=(DM*x+p0x)/(gamma*p);
                     dcosmax=min(1.0,(DM+p0x)/(gamma*p));

                     if(dcosmin < dcosmax){
                        spectral.DSG20R(dcosmin,dcosmax,ncosi,dcosid,ncosii);
                        for(int l=1;l<=ncosii;++l){
                           dcostheta=dcosid[l];
                           p3=p*dcostheta;
                             xpion=DM*x/(gamma*p3-p0x);

   fac0=q2/((gamma*q0*gamma*q0))*(p*p-p3*p3)/(2.0*mpi*mpi);
	 apq=p3*gamma-p0x;
	fac2=(p3*q2/(gamma*q0*q0*apq)+1.0)*(p3*q2/(gamma*q0*q0*apq)+1.00);
	 fact=(fac0+(apq*apq*fac2)/(mpi*mpi))*(mpi/(p3*gamma-p0x));

		   if(noWCut){goto g12;}

          else {if(IswcutforPionCloud(xpion,q2,p0x,p3,p,q0)){
//             std::cout<<"I am here hum"<<std::endl;

            g12:rdcosid[l]=-6.0*4.0*pi/xxpi3*2.0*mpi*fact*pionEffect.F2pion(xpion,aq2);
}//std::cout<<" I am true "<<rdcosid[l]<<std::endl;}
else{rdcosid[l]=0.0;
//std::cout<<" I am here for false "<<rdcosid[l]<<std::endl;}
}
 }                       }
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

double dshadowingF2(double x, double q2,std::string setname, int imem,int iset,const PDF* pdf)
{
    ShadowingEffect shadowF2;
    NucleonNLO Nucleon;
  double f2proton,f2neutron,f2Nucleon,shadow;
  double q2gev2=q2*GeVtfm*GeVtfm; // from fm to GeV
	 //double qgev=sqrt(q2gev2);
      double deltR=shadowF2.deltaR2_F2(x,q2gev2,setname,imem,iset,pdf);
      Nucleon.f2(x,q2gev2,setname,imem,iset,pdf,f2proton,f2neutron,f2Nucleon);
        if(noWCut){goto g12;}

         else{if(IsWacceptforShad(x,q2gev2)){
             g12:shadow=f2Nucleon*deltR;}
              else{shadow=0.0;}
}
         return shadow;
}

