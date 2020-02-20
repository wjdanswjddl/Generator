#include "Nucleon.hh"
#include "NuclearEffect.hh"
//#include "param.hh"
using namespace std;
//using namespace Param;

// aq2 in GeV
// by xx means x is modified in Nucleus enviroment.
//q2 is in fm unless is sepecified in GeV


void NucleonNLO::f1(double xx,double& aq2,const std::string setname, int& imem,int& iset,const PDF* pdf,double& f1Nucleon,double& f1Proton,double& f1Neutron)
{
      double x=xx;
      double q=aq2;
      double f2proton,f2neutron,f2Nucleon,flNucleon,flProton,flNeutron;
      if ((x < 0.0) || (x > 1.0)){
          f1Nucleon=0.0;
           f1Proton=0.0;
           f1Neutron=0.0;
           }
else{
 //     double aq2=q*qgev;
      double gamma2=1.0;
      //double gamma2=1.0+4.0*Mn*Mn*x*x/(q);
      f2(x,q,setname,imem,iset,pdf,f2proton,f2neutron,f2Nucleon);
      Fl(x,q,setname,imem,iset,pdf,flNucleon,flProton,flNeutron);
//      cout<<x<<" q= "<<q<<" f2f= "<<f2f<<" flf= "<<Flf<<endl;
      f1Nucleon=(gamma2*f2Nucleon-flNucleon)/(2.0*x);
       f1Proton=(gamma2*f2proton-flProton)/(2.0*x);
       f1Neutron=(gamma2*f2neutron-flNeutron)/(2.0*x);
//    cout<<f2neutron<<"  "<<flNeutron<<" "<<f1Neutron<<endl;
}
  }


void  NucleonNLO::Fl(double& x, double& aq2,const std::string setname, int& imem,int& iset,const PDF* pdf,double& flNucleon,double& flProton,double& flNeutron)
{
      double q=aq2;
     double yd[2000],ryps2[2000],ryds[2000],rygd[2000],ryns2[2000],rygd2[2000];
     double ryps2P[2000],rydsP[2000],rygdP[2000],ryns2P[2000],rygd2P[2000];
     double ryps2N[2000],rydsN[2000],rygdN[2000],ryns2N[2000],rygd2N[2000];
     if(x<=0.0||x>=1.0){
           flNucleon=0.0;
           flProton=0.0;
           flNeutron=0.0;
          	}
       SpectralFunction spectral;

        //const double zeta2=std::pow(pi,2)/6.0;
     	double as,xsinglet,xnonsinglet,xgluon,xinsinglet,xinnonsinglet,xingluon;
      double xinsingneutron,xinnonsingneutron,xingluon_N,xinsingproton,xinnonsingproton,xingluon_P,xsingproton,xnonsingproton,xgluon_P,xsingneutron,xnonsingneutron,xgluon_N;

      xf2nuc(x,q,as,xsinglet,xnonsinglet,xgluon,setname, imem, iset,pdf);
      xf2proton(x,aq2,as,xsingproton,xnonsingproton,xgluon_P,setname, imem, iset, pdf);
  xf2neutron(x,aq2,as,xsingneutron,xnonsingneutron,xgluon_N,setname, imem, iset, pdf);
      double xclns2_noint=0.012*(esquared*xsinglet + xnonsinglet);
      double xclns2_nointP=0.012*(esquared*xsingproton + xnonsingproton);
      double xclns2_nointN=0.012*(esquared*xsingneutron + xnonsingneutron);

      int npp;
      double a=x;
      double b=1.0;
      int n=1;
      spectral.DSG20R(a,b,n,yd,npp);
      for(int i=1;i<=npp;i++){
         double y=yd[i];
         double xin=x/y;
          xf2nuc(xin,q,as,xinsinglet,xinnonsinglet,xingluon,setname, imem, iset,pdf);
//lo
         double xclq1=4.0*Cf*y;
         double xclg1=8.0*nf*y*(1.0-y);
          ryds[i] = xclq1*(esquared*xinsinglet+xinnonsinglet);
          rygd[i] = xclg1*xingluon*esquared;

//nlo
         double y1=1.0-y;
         double xl1=log(1.0-y);
         double xl0=log(y);

         double xclps2=nf*((15.94-5.212*y)*y1*y1*xl1+ (0.421+1.520*y)*xl0*xl0+ 28.09*y1*xl0 -(2.370/y - 19.27)*pow(y1,3));
         ryps2[i] = xclps2*xinsinglet*esquared;

         double xclns2=(128.0/9.0)*y*xl1*xl1 - 46.50*y*xl1 - 84.094*xl0*xl1-37.338 + 89.53*y + 33.82*y*y + y*xl0*(32.90 + 18.410*xl0) - (128.0/9.0)*xl0+(16.0/27.0)*nf*(6.0*y*xl1 - 12.0*y*xl0-25.0*y + 6.0);
         ryns2[i]=xclns2*(esquared*xinsinglet+xinnonsinglet);


         double xclg2=nf*((94.740-49.20*y)*y1*xl1*xl1 + 864.80*y1*xl1 + 1161.00*y*xl1*xl0 +60.06*y*xl0*xl0 + 39.660*y1*xl0 - 5.3330*(1.00/y-1.0));

         rygd2[i] = xclg2*xingluon*esquared;

// Proton
 xf2proton(xin,aq2,as,xinsingproton,xinnonsingproton,xingluon_P,setname, imem, iset, pdf);

         double xclq1P=4.0*Cf*y;
         double xclg1P=8.0*nf*y*(1.0-y);
          rydsP[i] = xclq1P*(esquared*xinsingproton+xinnonsingproton);
          rygdP[i] = xclg1P*xingluon_P*esquared;

//nlo

         double xclps2P=nf*((15.94-5.212*y)*y1*y1*xl1+ (0.421+1.520*y)*xl0*xl0+ 28.09*y1*xl0 -(2.370/y - 19.27)*pow(y1,3));
         ryps2P[i] = xclps2P*xinsingproton*esquared;

         double xclns2P=(128.0/9.0)*y*xl1*xl1 - 46.50*y*xl1 - 84.094*xl0*xl1-37.338 + 89.53*y + 33.82*y*y + y*xl0*(32.90 + 18.410*xl0) - (128.0/9.0)*xl0+(16.0/27.0)*nf*(6.0*y*xl1 - 12.0*y*xl0-25.0*y + 6.0);
         ryns2P[i]=xclns2P*(esquared*xinsingneutron+xinnonsingproton);


         double xclg2P=nf*((94.740-49.20*y)*y1*xl1*xl1 + 864.80*y1*xl1 + 1161.00*y*xl1*xl0 +60.06*y*xl0*xl0 + 39.660*y1*xl0 - 5.3330*(1.00/y-1.0));

         rygd2P[i] = xclg2P*xingluon_P*esquared;

//Neutron

  xf2neutron(xin,aq2,as,xinsingneutron,xinnonsingneutron,xingluon_N,setname, imem, iset, pdf);

         double xclq1N=4.0*Cf*y;
         double xclg1N=8.0*nf*y*(1.0-y);
          rydsN[i] = xclq1N*(esquared*xinsingneutron+xinnonsingneutron);
          rygdN[i] = xclg1N*xingluon_N*esquared;

//nlo

         double xclps2N=nf*((15.94-5.212*y)*y1*y1*xl1+ (0.421+1.520*y)*xl0*xl0+ 28.09*y1*xl0 -(2.370/y - 19.27)*pow(y1,3));
         ryps2N[i] = xclps2N*xinsingneutron*esquared;

         double xclns2N=(128.0/9.0)*y*xl1*xl1 - 46.50*y*xl1 - 84.094*xl0*xl1-37.338 + 89.53*y + 33.82*y*y + y*xl0*(32.90 + 18.410*xl0) - (128.0/9.0)*xl0+(16.0/27.0)*nf*(6.0*y*xl1 - 12.0*y*xl0-25.0*y + 6.0);
         ryns2N[i]=xclns2N*(esquared*xinsingneutron+xinnonsingneutron);


         double xclg2N=nf*((94.740-49.20*y)*y1*xl1*xl1 + 864.80*y1*xl1 + 1161.00*y*xl1*xl0 +60.06*y*xl0*xl0 + 39.660*y1*xl0 - 5.3330*(1.00/y-1.0));

         rygd2N[i] = xclg2N*xingluon_N*esquared;



        }
      double rydsInt = spectral.DRG20R(a,b,n,ryds);
      double rygdInt = spectral.DRG20R(a,b,n,rygd);
      double ryps2Int = spectral.DRG20R(a,b,n,ryps2);
      double ryns2Int = spectral.DRG20R(a,b,n,ryns2);
      double rygd2Int = spectral.DRG20R(a,b,n,rygd2);


      double rydsIntP = spectral.DRG20R(a,b,n,rydsP);
      double rygdIntP = spectral.DRG20R(a,b,n,rygdP);
      double ryps2IntP = spectral.DRG20R(a,b,n,ryps2P);
      double ryns2IntP = spectral.DRG20R(a,b,n,ryns2P);
      double rygd2IntP = spectral.DRG20R(a,b,n,rygd2P);


      double rydsIntN = spectral.DRG20R(a,b,n,rydsN);
      double rygdIntN = spectral.DRG20R(a,b,n,rygdN);
      double ryps2IntN = spectral.DRG20R(a,b,n,ryps2N);
      double ryns2IntN = spectral.DRG20R(a,b,n,ryns2N);
      double rygd2IntN = spectral.DRG20R(a,b,n,rygd2N);



      double fllo=(rydsInt+rygdInt)*as;
       double flloP=(rydsIntP+rygdIntP)*as;
      double flloN=(rydsIntN+rygdIntN)*as;


      double higher_terms=as*as*(ryps2Int+xclns2_noint+ryns2Int+rygd2Int);
      double higher_termsP=as*as*(ryps2IntP+xclns2_nointP+ryns2IntP+rygd2IntP);
      double higher_termsN=as*as*(ryps2IntN+xclns2_nointN+ryns2IntN+rygd2IntN);

      if(iset==3){
      flNucleon=fllo;
       flProton=flloP;
      flNeutron=flloN;}
       flNucleon=fllo+higher_terms;
        flProton=flloP+higher_termsP;
        flNeutron=flloN+higher_termsN;



//      cout<<"x= "<<x<<fllo<<"  "<<flfree<<endl;

}


void NucleonNLO::f3(double xx, double& aq2,const std::string setname, int& imem,int& iset,const PDF* pdf,double& f3Nucleon,double& f3Proton,double& f3Neutron)
{
SpectralFunction spectral;

	double ll0d[2000],ll1d[2000];

	double yd[2000],rryd[2000],rry2d[2000],rrydP[2000],rry2dP[2000],rrydN[2000],rry2dN[2000];
        double q=aq2;
         double x=xx;
        const double zeta2=pow(pi,2)/6.0;
       double as,f3pdfx,f3pdfxin,f3pdfxP,f3pdfxN,f3pdfxinP,f3pdfxinN;
	if((x<0.0)||(x>1.0)){
         f3Nucleon=0.0;
         f3Proton=0.0;
         f3Neutron=0.0;
		}
        f3pdf(x,q,as,f3pdfx,f3pdfxP,f3pdfxN,setname,imem,iset,pdf);
         double f30loc0=f3pdfx;
         double f30loc0P=f3pdfxP;

         double f30loc0N=f3pdfxN;


	  if(iset==3){
        f3Nucleon=f30loc0;
        f3Proton=f30loc0P;
         f3Neutron=f30loc0N;}

{
        double xdl1=log(1.-x);
     double f30loc1x=(-(9.0+4.0*zeta2)+4.0*xdl1*xdl1/2.0 -3.0*xdl1)*f3pdfx;
     double f30loc1xP=(-(9.0+4.0*zeta2)+4.0*xdl1*xdl1/2.0 -3.0*xdl1)*f3pdfxP;
    double f30loc1xN=(-(9.0+4.0*zeta2)+4.0*xdl1*xdl1/2.0 -3.0*xdl1)*f3pdfxN;
 	  int nintp;
      double a=x;
      double b=1.0;
 	int nint=1;
	if(q<2.10)nint=nint+2;
	if(q<1.40)nint=nint+3;
	if(x<1.4e-3)nint=nint+3;
	if(x<0.7e-4)nint=nint+15;
	 spectral.DSG20R(a,b,nint,yd,nintp);
      	for(int i=1;i<=nintp;++i){
      	   double yy=yd[i];
               double xin=x/yy;
               f3pdf(xin,q,as,f3pdfxin,f3pdfxinP,f3pdfxinN,setname, imem, iset,pdf);

	         //double f3pdfxind[i]=f3pdfxin;
               ll0d[i]=log(yy);
               ll1d[i]=log(1.0-yy);
	         rryd[i]=(-2.0*(1.0+yy)*(ll1d[i]-ll0d[i])-4.0*ll0d[i]/(1.0-yy)+4.0+2.0*yy)/yy*f3pdfxin;

               rry2d[i]=(4.0*ll1d[i]-3.0)/(1.0-yy)*(f3pdfxin/yy-f3pdfx);

//Proton

rrydP[i]=(-2.0*(1.0+yy)*(ll1d[i]-ll0d[i])-4.0*ll0d[i]/(1.0-yy)+4.0+2.0*yy)/yy*f3pdfxinP;

               rry2dP[i]=(4.0*ll1d[i]-3.0)/(1.0-yy)*(f3pdfxinP/yy-f3pdfxP);

// Neutron

rrydN[i]=(-2.0*(1.0+yy)*(ll1d[i]-ll0d[i])-4.0*ll0d[i]/(1.0-yy)+4.0+2.0*yy)/yy*f3pdfxinN;

               rry2dN[i]=(4.0*ll1d[i]-3.0)/(1.0-yy)*(f3pdfxinN/yy-f3pdfxN);

              }
       double rrydint  = spectral.DRG20R(a,b,nint,rryd);
       double rry2dint = spectral.DRG20R(a,b,nint,rry2d);

       double rrydintP  = spectral.DRG20R(a,b,nint,rrydP);
       double rry2dintP = spectral.DRG20R(a,b,nint,rry2dP);

       double rrydintN  = spectral.DRG20R(a,b,nint,rrydN);
       double rry2dintN = spectral.DRG20R(a,b,nint,rry2dN);

       double factor=as;

       double f30loc1=f30loc1x*factor*Cf;
       double f30r=rrydint*factor*Cf;
       double f30i=rry2dint*factor*Cf;

//Proton
       double f30loc1P=f30loc1xP*factor*Cf;
       double f30rP=rrydintP*factor*Cf;
       double f30iP=rry2dintP*factor*Cf;

//Neutron

        double f30loc1N=f30loc1xN*factor*Cf;
       double f30rN=rrydintN*factor*Cf;
       double f30iN=rry2dintN*factor*Cf;



    	  f3Nucleon=f30loc0+f30loc1+f30i+f30r;
        f3Proton=f30loc0P+f30loc1P+f30iP+f30rP;
        f3Neutron=f30loc0N+f30loc1N+f30iN+f30rN;

}
//        cout<<"f3pdfxin "<<f3pdfxin<<"f3pdfx "<<f3pdfx<<endl;

 }

// here aq2 is in Gev
void NucleonNLO::f3pdf(double& x,double& aq2,double& as,double& f3pdfx,double& f3pdfxP,double& f3pdfxN,const std::string& setname, int& imem,int& iset,const PDF* pdf)
{




        double q=aq2;
        const double u = pdf->xfxQ2(2, x, q);
        const double d = pdf->xfxQ2(1, x, q);
        const double s = pdf->xfxQ2(3, x, q);
        const double c = pdf->xfxQ2(4, x, q);
        const double ubar = pdf->xfxQ2(-2, x, q);
        const double dbar = pdf->xfxQ2(-1, x, q);
        const double sbar = pdf->xfxQ2(-3, x, q);
        const double cbar = pdf->xfxQ2(-4, x, q);
          as = pdf->alphasQ2(q)/(4.0*pi);

	   f3pdfx=(u-ubar+d-dbar+2.0*s-2.0*cbar)/x;
          f3pdfxP=2.0*(u+c-dbar-sbar)/x;
          f3pdfxN=2.0*(d+c-ubar-sbar)/x;

}




void NucleonNLO::f2(double xx, double& aq2,const std::string& setname, int& imem,int& iset,const PDF* pdf,double& f2proton, double& f2neutron,double& f2Nucleon){

  double yd[2000],rydns[2000],ryds[2000],rygd[2000],ryd1ns[2000], ryd1s[2000];
  double protonrydns[2000],protonryds[2000],protonrygd[2000],protonryd1ns[2000], protonryd1s[2000];
  double neutronrydns[2000],neutronryds[2000],neutronrygd[2000],neutronryd1ns[2000], neutronryd1s[2000];

    SpectralFunction spectral;
     double x=xx;
    // cout<<"x= "<<x<<endl;
  double /*q2=q*q,*/xin,y;
  double  xinsinglet,xinnonsinglet, xingluon,xinsingneutron,xinnonsingneutron,xingluon_N,xinsingproton,xinnonsingproton,xingluon_P;
  double as,xsinglet,xnonsinglet, xgluon,f2localorder0,xsingproton,xnonsingproton,xgluon_P,xsingneutron,xnonsingneutron,xgluon_N;
  double f2localorder0_proton,f2localorder0_neutron;
  const double zeta2=std::pow(pi,2)/6.0;
  //double f2free;
  if(x<0.0||x>1.00){
    f2Nucleon=0.0;
    f2proton=0.0;
    f2neutron=0.0;
  }
 else{
  xf2nuc(x,aq2,as,xsinglet,xnonsinglet,xgluon,setname, imem, iset, pdf);
  xf2proton(x,aq2,as,xsingproton,xnonsingproton,xgluon_P,setname, imem, iset, pdf);
  xf2neutron(x,aq2,as,xsingneutron,xnonsingneutron,xgluon_N,setname, imem, iset, pdf);

//   cout<<as<<endl;
    f2localorder0=xnonsinglet+esquared*xsinglet;
   f2localorder0_proton=xnonsingproton+esquared*xsingproton;

   f2localorder0_neutron=xnonsingneutron+esquared*xsingneutron;

//   f2free=f2localorder0;
  if(iset==3){

   f2Nucleon=f2localorder0;
   f2proton= f2localorder0_proton;
   f2neutron=f2localorder0_neutron;
      }
 else{
  double F2localorder1=Cf*as*(4.0/2.0*pow(log(1.0-x),2)-3.0*log(1.0-x)-(9.0+4.0*zeta2))*(xnonsinglet+esquared*xsinglet);

 double F2localorder1_proton=Cf*as*(4.0/2.0*pow(log(1.0-x),2)-3.0*log(1.0-x)-(9.0+4.0*zeta2))*(xnonsingproton+esquared*xsingproton);

 double F2localorder1_neutron=Cf*as*(4.0/2.0*pow(log(1.0-x),2)-3.0*log(1.0-x)-(9.0+4.0*zeta2))*(xnonsingneutron+esquared*xsingneutron);

  int n=1;
  double a=x;
  double b=1.0;

  int nyy;
  spectral.DSG20R(a,b,n,yd,nyy);
  for(int i=1;i<=nyy;++i){

    y=yd[i];
    xin=x/y;

    xf2nuc(xin,aq2,as,xinsinglet,xinnonsinglet, xingluon,setname, imem, iset,pdf);

   //  rydns[i]=(4.0*log(1.0-y)/(1.0-y)-3.0/(1.0-y))*(xinnonsinglet-xnonsinglet);
    rydns[i]=0.0;


    ryds[i]=(4.0*log(1.0-y)/(1.0-y)-3.0/(1.0-y))*(xinsinglet-xsinglet);

    ryd1ns[i]=(-2.0*(1.0+y)*(log(1.0-y)-log(y))-4.0*log(y)/(1.0-y)+6.0+4.0*y)*xinnonsinglet;

    ryd1s[i]=(-2.0*(1.0+y)*(log(1.0-y)-log(y))-4.0*log(y)/(1.0-y)+6.0+4.0*y)*xinsinglet;

    rygd[i]=xingluon*((2.0-4.0*y*(1.0-y))*(log(1.0-y)-log(y))-2.+16.0*y*(1.0-y));

// proton
    xf2proton(xin,aq2,as,xinsingproton,xinnonsingproton,xingluon_P,setname, imem, iset, pdf);

    protonrydns[i]=(4.0*log(1.0-y)/(1.0-y)-3.0/(1.0-y))*(xinnonsingproton-xnonsingproton);
    //rydns[i]=0.0;


    protonryds[i]=(4.0*log(1.0-y)/(1.0-y)-3.0/(1.0-y))*(xinsingproton-xsingproton);

    protonryd1ns[i]=(-2.0*(1.0+y)*(log(1.0-y)-log(y))-4.0*log(y)/(1.0-y)+6.0+4.0*y)*xinnonsingproton;

    protonryd1s[i]=(-2.0*(1.0+y)*(log(1.0-y)-log(y))-4.0*log(y)/(1.0-y)+6.0+4.0*y)*xinsingproton;

    protonrygd[i]=xingluon_P*((2.0-4.0*y*(1.0-y))*(log(1.0-y)-log(y))-2.+16.0*y*(1.0-y));

// neutron
      xf2neutron(xin,aq2,as,xinsingneutron,xinnonsingneutron,xingluon_N,setname, imem, iset, pdf);

     neutronrydns[i]=(4.0*log(1.0-y)/(1.0-y)-3.0/(1.0-y))*(xinnonsingneutron-xnonsingneutron);
    //rydns[i]=0.0;


    neutronryds[i]=(4.0*log(1.0-y)/(1.0-y)-3.0/(1.0-y))*(xinsingneutron-xsingneutron);

    neutronryd1ns[i]=(-2.0*(1.0+y)*(log(1.0-y)-log(y))-4.0*log(y)/(1.0-y)+6.0+4.0*y)*xinnonsingneutron;

    neutronryd1s[i]=(-2.0*(1.0+y)*(log(1.0-y)-log(y))-4.0*log(y)/(1.0-y)+6.0+4.0*y)*xinsingneutron;

    neutronrygd[i]=xingluon_N*((2.0-4.0*y*(1.0-y))*(log(1.0-y)-log(y))-2.+16.0*y*(1.0-y));



  }
  double rydnsInt = spectral.DRG20R(a,b,n,rydns);
  double rydsInt = spectral.DRG20R(a,b,n,ryds);
  double ryd1nsInt = spectral.DRG20R(a,b,n,ryd1ns);
  double ryd1sInt = spectral.DRG20R(a,b,n,ryd1s);
  double rygdInt = spectral.DRG20R(a,b,n,rygd);
//proton
 double protonrydnsInt = spectral.DRG20R(a,b,n,protonrydns);
  double protonrydsInt = spectral.DRG20R(a,b,n,protonryds);
  double protonryd1nsInt = spectral.DRG20R(a,b,n,protonryd1ns);
  double protonryd1sInt = spectral.DRG20R(a,b,n,protonryd1s);
  double protonrygdInt = spectral.DRG20R(a,b,n,protonrygd);
//neutron
 double neutronrydnsInt = spectral.DRG20R(a,b,n,neutronrydns);
  double neutronrydsInt = spectral.DRG20R(a,b,n,neutronryds);
  double neutronryd1nsInt = spectral.DRG20R(a,b,n,neutronryd1ns);
  double neutronryd1sInt = spectral.DRG20R(a,b,n,neutronryd1s);
  double neutronrygdInt = spectral.DRG20R(a,b,n,neutronrygd);


  double F2Dnonsinglet = rydnsInt*Cf*as;
  double F2Dsinglet = rydsInt*Cf*as*esquared;
  double F2regularnonsinglet=ryd1nsInt*Cf*as;
  double F2regularsinglet=ryd1sInt*Cf*as*esquared;

  double F2gluon=rygdInt*esquared*as*nf;
  f2Nucleon=f2localorder0+F2localorder1+F2Dnonsinglet+F2Dsinglet+F2regularnonsinglet+F2regularsinglet+F2gluon;
//proton
 double F2Dnonsinglet_proton = protonrydnsInt*Cf*as;
  double F2Dsinglet_proton  = protonrydsInt*Cf*as*esquared;
  double F2regularnonsinglet_proton =protonryd1nsInt*Cf*as;
  double F2regularsinglet_proton=protonryd1sInt*Cf*as*esquared;

  double F2gluon_proton=protonrygdInt*esquared*as*nf;
 f2proton=f2localorder0_proton+F2localorder1_proton+F2Dnonsinglet_proton+F2Dsinglet_proton+F2regularnonsinglet_proton+F2regularsinglet_proton+F2gluon_proton;
//neutron
 double F2Dnonsinglet_neutron = neutronrydnsInt*Cf*as;
  double F2Dsinglet_neutron  = neutronrydsInt*Cf*as*esquared;
  double F2regularnonsinglet_neutron =neutronryd1nsInt*Cf*as;
  double F2regularsinglet_neutron=neutronryd1sInt*Cf*as*esquared;

  double F2gluon_neutron=neutronrygdInt*esquared*as*nf;
  f2neutron=f2localorder0_neutron+F2localorder1_neutron+F2Dnonsinglet_neutron+F2Dsinglet_neutron+F2regularnonsinglet_neutron+F2regularsinglet_neutron+F2gluon_neutron;
}

}
}



void NucleonNLO::xf2nuc(double& x,double& aq2,double& as,double& xsinglet,double& xnonsinglet,double& xgluon, const std::string& setname, int& imem,int& iset,const PDF* pdf){
  //double aq2=q*q;
//if((x < 0.0) || (x>1.0)){
  //         xsinglet=0.0;
//	}
//else{
if((x < 0.0) || (x>1.0)){
           xsinglet=0.0;
           xnonsinglet=0.0;

	}
else{

  double q=aq2;
  const double u = pdf->xfxQ2(2, x, q);
  const double d = pdf->xfxQ2(1, x, q);
  const double s = pdf->xfxQ2(3, x, q);
  const double c = pdf->xfxQ2(4, x, q);
  const double ubar = pdf->xfxQ2(-2, x, q);
  const double dbar = pdf->xfxQ2(-1, x, q);
  const double sbar = pdf->xfxQ2(-3, x, q);
  const double cbar = pdf->xfxQ2(-4, x, q);
  const double g = pdf->xfxQ2(0, x, q);
  as = pdf->alphasQ2(q)/(4.0*pi);
  xgluon=g;
  xsinglet=(u+ubar+d+dbar+s+sbar+c+cbar);

  xnonsinglet=0.0;
}
//  cout<<as<<endl;
}

void NucleonNLO::xf2proton(double& x,double& aq2,double& as,double& xsingproton,double& xnonsingproton,double& xgluon, const std::string& setname, int& imem,int& iset,const PDF* pdf){
if((x < 0.0) || (x>1.0)){
           xsingproton=0.0;
           xnonsingproton=0.0;

	}
else{

  //double aq2=q*q;
  double q=aq2;
  const double u = pdf->xfxQ2(2, x, q);
  const double d = pdf->xfxQ2(1, x, q);
  const double s = pdf->xfxQ2(3, x, q);
  const double c = pdf->xfxQ2(4, x, q);
  const double ubar = pdf->xfxQ2(-2, x, q);
  const double dbar = pdf->xfxQ2(-1, x, q);
  const double sbar = pdf->xfxQ2(-3, x, q);
  const double cbar = pdf->xfxQ2(-4, x, q);
  const double g = pdf->xfxQ2(0, x, q);
  as = pdf->alphasQ2(q)/(4.0*pi);
  xgluon=g;
  xsingproton=u+ubar+d+dbar+s+sbar+c+cbar;

//  xnonsingproton=0.0;
  xnonsingproton=u-ubar+dbar-d+sbar-s+c-cbar;
}
//  cout<<as<<endl;
}

void NucleonNLO::xf2neutron(double& x,double& aq2,double& as,double& xsingneutron,double& xnonsingneutron,double& xgluon, const std::string& setname, int& imem,int& iset,const PDF* pdf){

if((x < 0.0) || (x>1.0)){
           xsingneutron=0.0;
           xnonsingneutron=0.0;

	}
else{

  //double aq2=q*q;
  double q=aq2;

  const double u = pdf->xfxQ2(2, x, q);
  const double d = pdf->xfxQ2(1, x, q);
  const double s = pdf->xfxQ2(3, x, q);
  const double c = pdf->xfxQ2(4, x, q);
  const double ubar = pdf->xfxQ2(-2, x, q);
  const double dbar = pdf->xfxQ2(-1, x, q);
  const double sbar = pdf->xfxQ2(-3, x, q);
  const double cbar = pdf->xfxQ2(-4, x, q);
  const double g = pdf->xfxQ2(0, x, q);
  as = pdf->alphasQ2(q)/(4.0*pi);
  xgluon=g;
  xsingneutron=u+ubar+d+dbar+s+sbar+c+cbar;

  xnonsingneutron=d-dbar+ubar-u+sbar-s+c-cbar;
}
//  cout<<as<<endl;
}
// q2 in fm
double TargetMassEffect::gamman(double x, double q2)
{
double Factor1= (4.0*DM*DM*x*x)/q2;
return sqrt(1.0+Factor1);
}
//q2 in fm
double TargetMassEffect::sigma(double x, double q2)
{
return 2.0*x/(1.0+gamman(x,q2));
}




void TargetMassEffect::f1tmc(double xx, double& q2gev2,const std::string setname, int& imem,int& iset,const PDF* pdf,double& f1tmc_Nuc,double& f1tmc_P,double& f1tmc_N)
{
            //double q2=q*q;
             NucleonNLO Nucleon;
           double f1Nucleon,f1Proton,f1Neutron;
           double x=xx;
           double q2=q2gev2/(GeVtfm*GeVtfm);

           Nucleon.f1(sigma(x,q2),q2gev2,setname,imem,iset,pdf,f1Nucleon,f1Proton,f1Neutron);

	      double fac1=(x/(sigma(x,q2)*gamman(x,q2)))*f1Nucleon;
            double fac1P=(x/(sigma(x,q2)*gamman(x,q2)))*f1Proton;
            double fac1N=(x/(sigma(x,q2)*gamman(x,q2)))*f1Neutron;

            double af1=2.0*(DM*DM/q2)*x*sigma(x,q2)/gamman(x,q2);

            double af3;
           if(sigma(x,q2)>=0.0){
              af3=pow((1.0-sigma(x,q2)),2);
            }
           else af3=0.0;

            f1tmc_Nuc=fac1*(1.0+af1*af3);
            f1tmc_P=fac1P*(1.0+af1*af3);
            f1tmc_N=fac1N*(1.0+af1*af3);
	}

void TargetMassEffect::f2tmc(double& xx,double& q2gev2,const std::string setname, int& imem,int& iset,const PDF* pdf,double& f2P_TMC,double& f2N_TMC,double& f2NucTMC)
{
	     //double q2 = q*q;
           //double x2 = x*x;

       NucleonNLO Nucleon;
       double x=xx;
       double af2,f2proton,f2neutron,f2Nucleon;
       double q2=q2gev2/(GeVtfm*GeVtfm);
       Nucleon.f2(sigma(x,q2),q2gev2,setname,imem,iset,pdf,f2proton,f2neutron,f2Nucleon);

          double fac1=(x*x)/(sigma(x,q2)*sigma(x,q2)*pow(gamman(x,q2),3));
       double facNucl=fac1*f2Nucleon;
        double facProton=fac1*f2proton;
        double facNeutron=fac1*f2neutron;
      double af1=6.0*DM*DM/q2*x*sigma(x,q2)/gamman(x,q2);
      if(sigma(x,q2)>0.0){
	af2=(1.0-sigma(x,q2))*(1.0-sigma(x,q2));
      }
      else  {af2=0.0;}
       f2NucTMC=facNucl*(1.0+af1*af2);
       f2P_TMC=facProton*(1.0+af1*af2);
       f2N_TMC=facNeutron*(1.0+af1*af2);
}

void  TargetMassEffect::f3tmc(double xx, double& q2gev2,const std::string setname, int& imem,int& iset,const PDF* pdf,double& f3tmc_Nuc,double& f3tmc_P,double& f3tmc_N)
{
	// double q2 = q*q;
     double x=xx;
     NucleonNLO Nucleon;

       double q2=q2gev2/(GeVtfm*GeVtfm);

	 double af1,af2,fac1,f3Nucleon,f3Proton,f3Neutron,fac1P,fac1N;

         Nucleon.f3(sigma(x,q2),q2gev2,setname,imem,iset,pdf,f3Nucleon,f3Proton,f3Neutron);

       fac1=(x/(sigma(x,q2)*gamman(x,q2)*gamman(x,q2)))*f3Nucleon;
       fac1P=(x/(sigma(x,q2)*gamman(x,q2)*gamman(x,q2)))*f3Proton;
       fac1N=(x/(sigma(x,q2)*gamman(x,q2)*gamman(x,q2)))*f3Neutron;

	 const double DM2=pow(DM,2);

       if(sigma(x,q2)>=0.0){

       af1=(DM2/q2)*x*sigma(x,q2)/gamman(x,q2);
	 af2=(1.0-sigma(x,q2))*log(sigma(x,q2));

	 //f3ptmc=fac1*(1.0-af1*af2);
      }
      else
      { af1=0.0;
       af2=0.0; }
        f3tmc_Nuc=fac1*(1.0-af1*af2);
        f3tmc_P=fac1P*(1.0-af1*af2);
        f3tmc_N=fac1N*(1.0-af1*af2);

	}

//Higher Twist



void HigherTwist4::pdfs_nucleon(double& x,double& aq2,double& as,double& xsinglet,double& xnonsinglet,double& xgluon, const std::string& setname, int& imem,int& iset,const PDF* pdf){
  double q=aq2;

if((x < 0.0) || (x>1.0)){
           xsinglet=0.0;
           xnonsinglet=0.0;

	}
else {
 // double aq2=q*q;
  const double u = pdf->xfxQ2(2, x, q);
  const double d = pdf->xfxQ2(1, x, q);
  const double s = pdf->xfxQ2(3, x, q);
  const double c = pdf->xfxQ2(4, x, q);
  const double ubar = pdf->xfxQ2(-2, x, q);
  const double dbar = pdf->xfxQ2(-1, x, q);
  const double sbar = pdf->xfxQ2(-3, x, q);
  const double cbar = pdf->xfxQ2(-4, x, q);
  const double g = pdf->xfxQ2(0, x, q);
  as = pdf->alphasQ2(q)/(4.0*pi);
  xgluon=g;
  xsinglet=(u+ubar+d+dbar+s+sbar+c+cbar);

  xnonsinglet=(u-ubar+d-dbar)/x;
  }
  //cout<<<<endl;
}

void HigherTwist4::pdfs_proton(double& x,double& aq2,double& as,double& xsingproton,double& xnonsingproton,double& xgluon_P, const std::string& setname, int& imem,int& iset,const PDF* pdf){
  double q=aq2;

if((x < 0.0) || (x>1.0)){
           xsingproton=0.0;
           xnonsingproton=0.0;

	}
else {
 // double aq2=q*q;
  const double u = pdf->xfxQ2(2, x, q);
  const double d = pdf->xfxQ2(1, x, q);
  const double s = pdf->xfxQ2(3, x, q);
  const double c = pdf->xfxQ2(4, x, q);
  const double ubar = pdf->xfxQ2(-2, x, q);
  const double dbar = pdf->xfxQ2(-1, x, q);
  const double sbar = pdf->xfxQ2(-3, x, q);
  const double cbar = pdf->xfxQ2(-4, x, q);
  const double g = pdf->xfxQ2(0, x, q);
  as = pdf->alphasQ2(q)/(4.0*pi);
  xgluon_P=g;
  xsingproton=(u+ubar+d+dbar+s+sbar+c+cbar);

  xnonsingproton=(u-ubar+d-dbar)/x;
  }
  //cout<<<<endl;
}
void HigherTwist4::pdfs_neutron(double& x,double& aq2,double& as,double& xsingneutron,double& xnonsingneutron,double& xgluon_N, const std::string& setname, int& imem,int& iset,const PDF* pdf){
  double q=aq2;

if((x < 0.0) || (x>1.0)){
           xsingneutron=0.0;
           xnonsingneutron=0.0;

	}
else {
 // double aq2=q*q;
  const double u = pdf->xfxQ2(2, x, q);
  const double d = pdf->xfxQ2(1, x, q);
  const double s = pdf->xfxQ2(3, x, q);
  const double c = pdf->xfxQ2(4, x, q);
  const double ubar = pdf->xfxQ2(-2, x, q);
  const double dbar = pdf->xfxQ2(-1, x, q);
  const double sbar = pdf->xfxQ2(-3, x, q);
  const double cbar = pdf->xfxQ2(-4, x, q);
  const double g = pdf->xfxQ2(0, x, q);
  as = pdf->alphasQ2(q)/(4.0*pi);
  xgluon_N=g;
  xsingneutron=(u+ubar+d+dbar+s+sbar+c+cbar);

  xnonsingneutron=(u-ubar+d-dbar)/x;
  }
  //cout<<<<endl;
}

void HigherTwist4::f2ht(double& xx, double& aq2,const std::string setname, int& imem,int& iset,const PDF* pdf,double& f2htNucleon, double& f2htProton, double& f2htNeutron)
{
             TargetMassEffect TMC;
            double x=xx;
            double f2P_TMC,f2N_TMC,f2NucTMC,htNucleon,htProton,htNeutron;
           //double qgev=sqrt(aq2);
          htns_nucleon(x,aq2,setname,imem,iset,pdf,htNucleon,htProton,htNeutron);

          double f2ht=htNucleon/aq2;
          double f2htp=htProton/aq2;
           double f2htn=htNeutron/aq2;
        TMC.f2tmc(x,aq2,setname,imem,iset,pdf,f2P_TMC,f2N_TMC,f2NucTMC);


           if(x <= 1.0 && x >= 0.0){

           f2htNucleon=f2NucTMC*(1.00+f2ht);
            f2htProton=f2P_TMC*(1.00+f2htp);
            f2htNeutron=f2N_TMC*(1.00+f2htn);}
           else{f2htNucleon=0.0;
             f2htProton=0.0;
              f2htNeutron=0.0; }
//            cout<<f1tmcE<<endl;
 }

void HigherTwist4::htns_nucleon(double& x, double& aq2,const std::string setname, int& imem,int& iset,const PDF* pdf,double& htNucleon,double& htProton,double& htNeutron)
{
          SpectralFunction spectral;

          double q=aq2;
           //double rdd[2000];
           double  xd[2000],rd[2000],rd1[2000];
           double  rdp[2000],rd1p[2000];
           double  rdn[2000],rd1n[2000];

             int np;
            int nm=1;

              double xsingneutron,xnonsingneutron=0.0,xgluon_N,xesingneutron,xenonsingneutron,xegluon_N,Pterm=0.0,Nterm=0.0;
             double xinsingneutron,xinnonsingneutron=0.0,xingluon_N,x1singneutron,x1nonsingneutron=0.0,x1gluon_N=0.0;

             double xsingproton,xnonsingproton=0.0,xgluon_P,xesingproton,xenonsingproton,xegluon_P;
             double xinsingproton,xinnonsingproton=0.0,xingluon_P,x1singproton,x1nonsingproton=0.0,x1gluon_P=0.0;

              double as,xsinglet,xnonsinglet=0.0,xgluon,xesinglet,xenonsinglet,xegluon,term=0.0;
             double xinsinglet,xinnonsinglet=0.0,xingluon,x1singlet,x1nonsinglet=0.0,x1gluon=0.0;

             double forth,disc,htheta,third,first,second,dirac,Pforth,Nforth,Pthird,Nthird,Pfirst,Nfirst,Psecond,Nsecond;
            if(x<1.0){

  pdfs_nucleon(x,q,as,xsinglet,xnonsinglet,xgluon,setname, imem, iset, pdf);
 pdfs_proton(x,q,as,xsingproton,xnonsingproton,xgluon_P,setname, imem, iset, pdf);
 pdfs_neutron(x,q,as,xsingneutron,xnonsingneutron,xgluon_N,setname, imem, iset, pdf);

          double epsi=1e-3;

          double xe=x+epsi;

  pdfs_nucleon(xe,q,as,xesinglet,xenonsinglet,xegluon,setname, imem, iset, pdf);
 pdfs_proton(xe,q,as,xesingproton,xenonsingproton,xegluon_P,setname, imem, iset, pdf);
 pdfs_neutron(xe,q,as,xesingneutron,xenonsingneutron,xegluon_N,setname, imem, iset, pdf);

            term=(xenonsinglet-xnonsinglet)/epsi;
            Pterm=(xenonsingproton-xnonsingproton)/epsi;
            Nterm=(xenonsingneutron-xnonsingneutron)/epsi;



             disc=0.0;
             dirac=0.0;
              htheta=0.0;

            double x1=0.999990;
 pdfs_nucleon(x1,q,as,x1singlet,x1nonsinglet,x1gluon,setname, imem, iset, pdf);
 pdfs_proton(x1,q,as,x1singproton,x1nonsingproton,x1gluon_P,setname, imem, iset, pdf);
 pdfs_neutron(x1,q,as,x1singneutron,x1nonsingneutron,x1gluon_N,setname, imem, iset, pdf);


            forth=-x1nonsinglet*dirac+xnonsinglet+term*x*(1.0+disc-2.0*htheta);
            Pforth=-x1nonsingproton*dirac+xnonsingproton+Pterm*x*(1.0+disc-2.0*htheta);
            Nforth=-x1nonsingneutron*dirac+xnonsingneutron+Nterm*x*(1.0+disc-2.0*htheta);
}
            else if(x==1.0){

            disc=1.0;
            htheta=0.5;
            dirac=1.0;

          forth=x1nonsinglet*dirac+xnonsinglet+term*x*(1.0+disc-2.0*htheta);
          Pforth=x1nonsingproton*dirac+xnonsingproton+Pterm*x*(1.0+disc-2.0*htheta);
          Nforth=x1nonsingneutron*dirac+xnonsingneutron+Nterm*x*(1.0+disc-2.0*htheta);


            }
            else if(x>1.0){
            forth=0.0;
            Pforth=0.0;
            Nforth=0.0;}



          third=-9.0*xnonsinglet;
          Pthird=-9.0*xnonsingproton;
          Nthird=-9.0*xnonsingneutron;

            double xm=1.0;
           spectral.DSG20R(x,xm,nm,xd,np);
           for(int iy=1;iy<=np;++iy){
           double y=xd[iy];
           double xin=x/y;

           pdfs_nucleon(xin,q,as,xinsinglet,xinnonsinglet, xingluon,setname, imem, iset,pdf);
            pdfs_proton(xin,q,as,xinsingproton,xinnonsingproton, xingluon_P,setname, imem, iset,pdf);
            pdfs_neutron(xin,q,as,xinsingproton,xinnonsingneutron, xingluon_N,setname, imem, iset,pdf);

           rd1[iy]=xinnonsinglet*2.0*(2.0+y+6.0*y*y)/y;
           rd1p[iy]=xinnonsingproton*2.0*(2.0+y+6.0*y*y)/y;
           rd1n[iy]=xinnonsingneutron*2.0*(2.0+y+6.0*y*y)/y;
           }

           double res1=spectral.DRG20R(x,xm,nm,rd1);
           double res1p=spectral.DRG20R(x,xm,nm,rd1p);
           double res1n=spectral.DRG20R(x,xm,nm,rd1n);
           double xup=1.0;
           spectral.DSG20R(x,xup,nm,xd,np);
            for(int iy=1;iy<=np;++iy){
           double y=xd[iy];
           double xin=x/y;

            pdfs_nucleon(xin,q,as,xinsinglet,xinnonsinglet, xingluon,setname, imem, iset,pdf);
            pdfs_proton(xin,q,as,xinsingproton,xinnonsingproton, xingluon_P,setname, imem, iset,pdf);
            pdfs_neutron(xin,q,as,xinsingneutron,xinnonsingneutron, xingluon_N,setname, imem, iset,pdf);
          rd[iy]=(xinnonsinglet/y-xnonsinglet)*(1.0/(1.0-y));
          rdp[iy]=(xinnonsingproton/y-xnonsingproton)*(1.0/(1.0-y));
          rdn[iy]=(xinnonsingneutron/y-xnonsingneutron)*(1.0/(1.0-y));
           }

           double res= spectral.DRG20R(x,xup,nm,rd);
           double resp= spectral.DRG20R(x,xup,nm,rdp);
           double resn= spectral.DRG20R(x,xup,nm,rdn);
           first=-4.0*(res+xnonsinglet*log(1.0-x));
            second=res1;

            Pfirst=-4.0*(resp+xnonsingproton*log(1.0-x));
            Psecond=res1p;

            Nfirst=-4.0*(resn+xnonsingneutron*log(1.0-x));
            Nsecond=res1n;
            double a2pp=0.0;
          if(aq2<10.0){
          a2pp=-0.15;}
          else if(aq2 == 10.0){
          a2pp=-0.2;}
          else if(aq2 > 10.0){
          a2pp=-0.19;}
           //double xnsv;
           htNucleon=(first+second+third+forth)*a2pp/xnonsinglet;
           htProton=(Pfirst+Psecond+Pthird+Pforth)*a2pp/xnonsingproton;
           htNeutron=(Nfirst+Nsecond+Nthird+Nforth)*a2pp/xnonsingneutron;


}



/////////////// F1 Higher Twist

void HigherTwist4::f1ht(double& xx, double& aq2,const std::string setname, int& imem,int& iset,const PDF* pdf,double& f1htNucleon, double& f1htProton, double& f1htNeutron)
{
             TargetMassEffect TMC;
            double h1tNucleon,h1tProton,h1tNeutron,f1tmc_Nuc,f1tmc_P,f1tmc_N;
            double x=xx;
           //double qgev=sqrt(aq2);
           htns_f1(x,aq2,setname,imem,iset,pdf,h1tNucleon,h1tProton,h1tNeutron);

           double f1_ns=x*h1tNucleon;
            double f1_ht=f1_ns/(2.0*x);

           double f1_nsP=x*h1tProton;
            double f1_htP=f1_nsP/(2.0*x);


            double f1_nsN=x*h1tNeutron;
            double f1_htN=f1_nsN/(2.0*x);



           double twist4_f1=f1_ht/aq2;
           double twist4_f1P=f1_htP/aq2;
           double twist4_f1N=f1_htN/aq2;


           TMC.f1tmc(x,aq2,setname,imem,iset,pdf,f1tmc_Nuc,f1tmc_P,f1tmc_N);

           if(x <= 1.0 && x >= 0.0){

           f1htNucleon=f1tmc_Nuc*(1.00+twist4_f1);
            f1htProton=f1tmc_P*(1.00+twist4_f1P);
             f1htNeutron=f1tmc_N*(1.00+twist4_f1N); }
           else{f1htNucleon=0.0;
            f1htProton=0.0;
              f1htNeutron=0.0;}
//            cout<<f1tmcE<<endl;
       }





void HigherTwist4::htns_f1(double x, double aq2,std::string setname, int imem,int iset,const PDF* pdf,double& h1tNucleon,double& h1tProton,double& h1tNeutron)
{
           double xd[2000];
           double  rd[2000],rd1[2000], rdp[2000],rd1p[2000], rdn[2000],rd1n[2000];
           //double xd1[2000],xd2[2000];
SpectralFunction spectral;

  double as,xsinglet,xnonsinglet=0.0,xgluon,xesinglet,xenonsinglet=0.0,xegluon,term=0.0,a2pp=0.0;
 double  xinsinglet,xinnonsinglet=0.0,xingluon,x1singlet,x1nonsinglet=0.0,x1gluon=0.0;

   double xsingneutron,xnonsingneutron=0.0,xgluon_N,xesingneutron,xenonsingneutron,xegluon_N,Pterm=0.0,Nterm=0.0;
             double xinsingneutron,xinnonsingneutron=0.0,xingluon_N,x1singneutron,x1nonsingneutron=0.0,x1gluon_N=0.0;

             double xsingproton,xnonsingproton=0.0,xgluon_P,xesingproton,xenonsingproton,xegluon_P;
             double xinsingproton,xinnonsingproton=0.0,xingluon_P,x1singproton,x1nonsingproton=0.0,x1gluon_P=0.0;
	double forth,disc,htheta,third,first,second,dirac,Pforth,Nforth,Pthird,Nthird,Pfirst,Nfirst,Psecond,Nsecond;

          // double forth,disc,htheta,third,first,second,dirac;

            double q=aq2;
            int np;
            int nm=1;


         if(x<1.0){

  pdfs_nucleon(x,q,as,xsinglet,xnonsinglet,xgluon,setname, imem, iset, pdf);
pdfs_proton(x,q,as,xsingproton,xnonsingproton,xgluon_P,setname, imem, iset, pdf);
 pdfs_neutron(x,q,as,xsingneutron,xnonsingneutron,xgluon_N,setname, imem, iset, pdf);
          double epsi=1e-3;

          double xe=x+epsi;

  pdfs_nucleon(xe,q,as,xesinglet,xenonsinglet,xegluon,setname, imem, iset, pdf);
pdfs_proton(xe,q,as,xesingproton,xenonsingproton,xegluon_P,setname, imem, iset, pdf);
 pdfs_neutron(xe,q,as,xesingneutron,xenonsingneutron,xegluon_N,setname, imem, iset, pdf);


           term=(xenonsinglet-xnonsinglet)/epsi;
            Pterm=(xenonsingproton-xnonsingproton)/epsi;
           Nterm=(xenonsingneutron-xnonsingneutron)/epsi;




             disc=0.0;
             dirac=0.0;
              htheta=0.0;

            double x1=0.999990;
        pdfs_nucleon(x1,q,as,x1singlet,x1nonsinglet,x1gluon,setname, imem, iset, pdf);
 pdfs_proton(x1,q,as,x1singproton,x1nonsingproton,x1gluon_P,setname, imem, iset, pdf);
 pdfs_neutron(x1,q,as,x1singneutron,x1nonsingneutron,x1gluon_N,setname, imem, iset, pdf);


           forth=-x1nonsinglet*dirac+xnonsinglet+term*x*(1.0+disc-2.0*htheta);
           Pforth=-x1nonsingproton*dirac+xnonsingproton+Pterm*x*(1.0+disc-2.0*htheta);
           Nforth=-x1nonsingneutron*dirac+xnonsingneutron+Nterm*x*(1.0+disc-2.0*htheta);


}
            else if(x==1.0){

            disc=1.0;
            htheta=0.5;
            dirac=1.0;

         forth=-x1nonsinglet*dirac+xnonsinglet+term*x*(1.0+disc-2.0*htheta);
           Pforth=-x1nonsingproton*dirac+xnonsingproton+Pterm*x*(1.0+disc-2.0*htheta);
           Nforth=-x1nonsingneutron*dirac+xnonsingneutron+Nterm*x*(1.0+disc-2.0*htheta);


            }


           else if(x>1.0){
            forth=0.0;
            Pforth=0.0;
            Nforth=0.0;}

           third=-5.0*xnonsinglet;
           Pthird=-5.0*xnonsingproton;
           Nthird=-5.0*xnonsingneutron;


            double xm=1.0;
           spectral.DSG20R(x,xm,nm,xd,np);
           for(int iy=1;iy<=np;++iy){
           double y=xd[iy];
           double xin=x/y;

      pdfs_nucleon(xin,q,as,xinsinglet,xinnonsinglet, xingluon,setname, imem, iset,pdf);
      pdfs_proton(xin,q,as,xinsingproton,xinnonsingproton, xingluon_P,setname, imem, iset,pdf);
      pdfs_neutron(xin,q,as,xinsingproton,xinnonsingneutron, xingluon_N,setname, imem, iset,pdf);

            rd1[iy]=xinnonsinglet*2.0*(2.0+y+2.0*y*y)/y;
            rd1p[iy]=xinnonsingproton*2.0*(2.0+y+2.0*y*y)/y;
            rd1n[iy]=xinnonsingneutron*2.0*(2.0+y+2.0*y*y)/y;



           }

          double res1=spectral.DRG20R(x,xm,nm,rd1);
          double res1p=spectral.DRG20R(x,xm,nm,rd1p);
          double res1n=spectral.DRG20R(x,xm,nm,rd1n);

           double xup=1.0;
           spectral.DSG20R(x,xup,nm,xd,np);
            for(int iy=1;iy<=np;++iy){
           double y=xd[iy];
           double xin=x/y;

          pdfs_nucleon(xin,q,as,xinsinglet,xinnonsinglet, xingluon,setname, imem, iset,pdf);
  pdfs_proton(xin,q,as,xinsingproton,xinnonsingproton, xingluon_P,setname, imem, iset,pdf);
 pdfs_neutron(xin,q,as,xinsingneutron,xinnonsingneutron, xingluon_N,setname, imem, iset,pdf);

         rd[iy]=(xinnonsinglet/y-xnonsinglet)*(1.0/(1.0-y));
          rdp[iy]=(xinnonsingproton/y-xnonsingproton)*(1.0/(1.0-y));
          rdn[iy]=(xinnonsingneutron/y-xnonsingneutron)*(1.0/(1.0-y));
           }

           double res= spectral.DRG20R(x,xup,nm,rd);
           double resp= spectral.DRG20R(x,xup,nm,rdp);
           double resn= spectral.DRG20R(x,xup,nm,rdn);



//       ____________________________________________________________

        first=-4.0*(res+xnonsinglet*log(1.0-x));
        Pfirst=-4.0*(resp+xnonsingproton*log(1.0-x));
        Nfirst=-4.0*(resn+xnonsingneutron*log(1.0-x));

            second=res1;
             Psecond=res1p;
             Nsecond=res1n;

          if(aq2<10.0){
          a2pp=-0.15;}
          else if(aq2 == 10.0){
          a2pp=-0.2;}
          else if(aq2 > 10.0){
          a2pp=-0.19;}
           //double cch=(first+second+third+forth)*a2pp/xnonsinglet;
//           cout<<cch<<endl;
           h1tNucleon=(first+second+third+forth)*a2pp/xnonsinglet;
           h1tProton=(Pfirst+Psecond+Pthird+Pforth)*a2pp/xnonsingproton;
            h1tNeutron=(Nfirst+Nsecond+Nthird+Nforth)*a2pp/xnonsingneutron;


}


/////// F3 Higher Twist

void HigherTwist4::f3ht(double& xx, double& aq2,const std::string setname, int& imem,int& iset,const PDF* pdf,double& f3htNucleon, double& f3htProton, double& f3htNeutron)
{
      TargetMassEffect TMC;
             double x=xx;
            double f3tmc_Nuc,f3tmc_P,f3tmc_N,h1tNucleon,h1tProton,h1tNeutron;
           //double qgev=sqrt(aq2);
            htns_f1(x,aq2,setname,imem,iset,pdf,h1tNucleon,h1tProton,h1tNeutron);

           double twist4Nuc=h1tNucleon/aq2;
           double twist4P=h1tProton/aq2;
           double twist4N=h1tNeutron/aq2;


         TMC.f3tmc(x,aq2,setname,imem,iset,pdf,f3tmc_Nuc,f3tmc_P,f3tmc_N);


           if(x <= 1.0 && x >= 0.0){

           f3htNucleon=f3tmc_Nuc*(1.00+twist4Nuc);
           f3htProton=f3tmc_P*(1.00+twist4P);
           f3htNeutron=f3tmc_N*(1.00+twist4N);

}
           else{f3htNucleon=0.0;
               f3htProton=0.0;
                f3htNeutron=0.0;}

}




double ForShadowinfF1::f1nucshad(double x,double q,std::string setname, int imem,int iset,const PDF* pdf)
{
   NucleonNLO Nucleon;
   double F1freshad,f2proton,f2neutron,f2Nucleon;
      if ((x < 0.0) || (x > 1.0)){
         F1freshad=0.0;
         return F1freshad;
      }
 //     double aq2=q*qgev;
      double gamma2=1.0+4.0*Mn*Mn*x*x/(q);
      Nucleon.f2(x,q,setname,imem,iset,pdf,f2proton,f2neutron,f2Nucleon);
      double Flf=Flongitudinal(x,q,setname,imem,iset,pdf);


      F1freshad=(gamma2*f2Nucleon-Flf)/(2.0*x);
      return F1freshad;
  }

double  ForShadowinfF1::Flongitudinal(double x, double q,std::string setname, int imem,int iset,const PDF* pdf)
{
SpectralFunction spectral;
NucleonNLO Nucleon;
      double yd[2000],ryd1[2000],ryd2[2000],rygd[2000];
       double as,xinsinglet,xinnonsinglet, xingluon, Flongi;
      if ((x < 0.0) || (x > 1.0)){
         Flongi=0.0;
         return Flongi;
      }

  //   double aq2=qgev*qgev;

//c Beginning of the Mellin Convolution for the singlet quark
//c distribution, the non-singlet quark distribution
//c and the gluon distribution
      int n=1;
      int nyy;
      double xup=1.0;
       double xin,y;
      spectral.DSG20R(x,xup,n,yd,nyy);
      for(int iy=1;iy<=nyy;++iy){
         y=yd[iy];
         xin=x/y;
    Nucleon.xf2nuc(xin,q,as,xinsinglet,xinnonsinglet, xingluon,setname, imem, iset,pdf);
         ryd1[iy]=Cf*4.0*y*xinsinglet;
         ryd2[iy]=Cf*4.0*y*xinnonsinglet;

         rygd[iy]=nf*8.0*y*(1.0-y)*xingluon;
      }

      double res1 = spectral.DRG20R(x,xup,n,ryd1);
      double res2 = spectral.DRG20R(x,xup,n,ryd2);
      double resglu = spectral.DRG20R(x,xup,n,rygd);

      double Flong_nonsingl_order1=res2*as;
      double Flong_singl_order1=res1*as*esquared;
     double  Flong_glu_order1=resglu*as*esquared;


    Flongi=Flong_singl_order1+Flong_glu_order1+Flong_nonsingl_order1;

      return Flongi;
}
