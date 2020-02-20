#include <iostream>
#include <string>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <algorithm>
#include "param.hh"
#include "NuclearEffect.hh"
#include <fstream>
using namespace std;
using namespace Param;

namespace onlyForRho0
{
  double rho0A;
//  double rho0p;
//   double rho0n;


}

namespace onlyForRho0p
{
  double rho0p;

}

namespace onlyForRho0n
{
  double rho0n;

}


namespace Densidad
{
  double DMDF;
  double dro;
}

namespace Densidadx
{
  double drox;
}

namespace Referd{
  double drefer;

}

namespace SigmaShad
{
  double sigmaT;
}


//     hole spectral function
// UNITS: fm
//input: DDDK0,DDDK, dro -> Energy, momentum, density
// output:  DSh1, dmu -> spectral function, chemical potential



void SpectralFunction::DSh1int(double& DDDK0,double& DDDK,double& drho,double& DMU, double& spect)
{
  //    EXPRESION DE LA ENERGIA DE FERMI

  //double drox=Density(r);
  using namespace Densidadx;
  using namespace Densidad;
  using namespace Referd;
  Densidad::dro=drho;
  Densidadx::drox=drho;


  double aux=pow((3.0/2.0*pi*pi*drho),1.0/3.0);
  Densidad::DMDF = aux;


  double DMUx = DMDF*DMDF/(2.0*DM);
  double drefer1=reasi_int(DMUx,aux);
  Referd::drefer=drefer1;
  DMU=DMUx+drefer;
  //   double DIMA1;

  double DIMA1,DTARDE,arg1,aux1;
  double deno=sqrt(DDDK*DDDK + DM*DM);
  double fac=DM/deno;

  if(DDDK0 < DMU)
    {
      arg1=DDDK0-drefer;
      DTARDE=DIMASM(arg1,DDDK)*fac;
      aux1=DDDK0-DDDK*DDDK/2.0/DM - reasi_int(arg1,DDDK);
      DIMA1=DTARDE/(DTARDE*DTARDE + aux1*aux1);

    }
  else {
    DIMA1=0.0;
  }

  spect = DIMA1/pi;
}

void SpectralFunction::DSh1int_proton(double& DDDK0,double& DDDK,double& drho,double& DMU,double& spect)
{
  //    EXPRESION DE LA ENERGIA DE FERMI

  //double drox=Density(r);
  using namespace Densidadx;
  using namespace Densidad;
  using namespace Referd;
  Densidad::dro=drho;
  Densidadx::drox=drho;


  double aux=pow((3.0*pi*pi*drho),1.0/3.0);
  Densidad::DMDF = aux;


  double DMUx = DMDF*DMDF/(2.0*DM);
  double drefer1=reasi_int(DMUx,aux);
  Referd::drefer=drefer1;
  DMU=DMUx+drefer;
  //   double DIMA1;

  double DIMA1,DTARDE,arg1,aux1;
  double deno=sqrt(DDDK*DDDK + DM*DM);
  double fac=DM/deno;

  if(DDDK0 < DMU)
    {
      arg1=DDDK0-drefer;
      DTARDE=DIMASM(arg1,DDDK)*fac;
      aux1=DDDK0-DDDK*DDDK/2.0/DM - reasi_int(arg1,DDDK);
      DIMA1=DTARDE/(DTARDE*DTARDE + aux1*aux1);

    }
  else {
    DIMA1=0.0;
  }

  spect = DIMA1/pi;
}

void SpectralFunction::DSh1int_neutron(double& DDDK0,double& DDDK,double& drho,double& DMU,double& spect)
{
  //    EXPRESION DE LA ENERGIA DE FERMI

  //double drox=Density(r);
  using namespace Densidadx;
  using namespace Densidad;
  using namespace Referd;
  Densidad::dro=drho;
  Densidadx::drox=drho;


  double aux=pow((3.0*pi*pi*drho),1.0/3.0);
  Densidad::DMDF = aux;


  double DMUx = DMDF*DMDF/(2.0*DM);
  double drefer1=reasi_int(DMUx,aux);
  Referd::drefer=drefer1;
  DMU=DMUx+drefer;
  //   double DIMA1;

  double DIMA1,DTARDE,arg1,aux1;
  double deno=sqrt(DDDK*DDDK + DM*DM);
  double fac=DM/deno;

  if(DDDK0 < DMU)
    {
      arg1=DDDK0-drefer;
      DTARDE=DIMASM(arg1,DDDK)*fac;
      aux1=DDDK0-DDDK*DDDK/2.0/DM - reasi_int(arg1,DDDK);
      DIMA1=DTARDE/(DTARDE*DTARDE + aux1*aux1);

    }
  else {
    DIMA1=0.0;
  }

  spect = DIMA1/pi;
}


//      FUNCIONES DE LINHARD

void SpectralFunction::ULIND(double& QZR,double& Q,double& XKF,std::complex<double>& CUFUN)
{
  double DCOCI,RUN,YMN,B,RO,A,AP,AM,TDIR,TCROS,TERS,TERF,SQS,QFR,GAMH,FAC;
  std::complex<double> CUNUC,CUDEL,CAC,CAPC,CDIR,CCROS;
  //double Q=q;
  double XMN=6.73480;
  double XMD=8.8315;
  double WRES=XMD-XMN;
  double FNS=1.0;
  double FDS=4.520;
  std::complex<double> CYI (0.0,1.0);
  double QA=Q/XKF;
  double QZ=abs(QZR);
  double QZA=abs(QZ*XMN/pow(XKF,2));
   double pi2=pow(pi,2);
  //https://www.programiz.com/c-programming/c-break-continue-statement

  if(Q==0.0)
    {
      RUN=FNS*XMN*XKF/(pi*pi)*2.0/3.0*pow(QA/QZA,2);
    }
  //while(Q!=0.0)
  else
    {
      AM=((abs(1.0- (QZA/QA-QA/2.0)) < 0.000010) ? 1.000010 :  QZA/QA-QA/2.0) ;
      AP=((abs(1.0- (QZA/QA+QA/2.0)) < 0.000010) ? 1.000010 :  QZA/QA+QA/2.0) ;
      TERF=abs((1.0+AM)/(1.0-AM))+1.e-15;
      TERS=abs((1.0+AP)/(1.0-AP))+1.e-15;
      if(QZA-0.000010 <= 0.0)
	{
	  RUN=FNS*XMN*XKF/pi2*(-1.0+(1.0-AM*AM)/
				     (2.0*QA)*log(TERF)-(1.0-AP*AP)/(2.0*QA)*log(TERS));
	}
      else
	{
	  DCOCI=abs(QA/QZA);
	  RUN=(( DCOCI-0.1 <= 0.) ? FNS*XMN*XKF/(pi*pi)*2.0/3.0*pow(QA/QZA,2) :
	       FNS*XMN*XKF/pi2*(-1.0+(1.0-AM*AM)/
				      (2.0*QA)*log(TERF)-(1.0-AP*AP)/(2.0*QA)*log(TERS)) );
	}
    }
  YMN=0.0;
  if( ((QA*QA/2.0+QA)>=QZA) && (QZA>=(QA*QA/2.0-QA)) && (QA>=2.0) )
    {
      YMN=-2.0*XMN*XKF*FNS/(4.0*pi*QA)*(1.0-pow(AM,2));
    }
  else if((QA < 2.0) && ((QA+QA*QA/2.0) >= QZA) && (QZA >= (QA-QA*QA/2.0)))
    {
      YMN=-2.0*XMN*XKF*FNS/(4.0*pi*QA)*(1.0-pow(AM,2));
    }
  else if((QA < 2.0) && (0.0 <= QZA) && (QZA <= (QA-QA*QA/2.0)))
    {
      YMN=-2.0*XMN*XKF*FNS/(4.0*pi*QA)*2.0*QZA;
    }

  CUNUC = RUN + CYI*YMN;
  B=Q/XMD;
  RO=pow(XKF,3)*2.0/(3.0*pi*pi);
  FAC= ((Q == 0. ) ? 0.: pow((4.0/3.0*XKF/(2.0*pi)),2)/pow(B,3)*FDS );
  if((abs(QZ)-1.0) <= 0.)
    {
      A=(QZ-WRES-pow(Q,2)/(2.0*XMD))/XKF;
      AP=(-QZ-WRES-pow(Q,2)/(2.0*XMD))/XKF;
      if((abs(B/A)-0.10) <= 0.)
	{
	  TDIR=4.0/9.0*FDS*RO/(A*XKF);
	}
      else
	{
	  TERF=abs((A+B)/(A-B))+1.e-15;
	  TDIR=FAC*(B*A+(pow(B,2)-pow(A,2))/2.0*log(TERF));
	}
      //--------------
      if((abs(B/AP)-0.10) <= 0.)
	{
	  TCROS=4.0/9.0*FDS*RO/(AP*XKF);
        CUDEL=TDIR+TCROS;
        goto ll;
	}
      else
	{
	  TERS=abs((AP-B)/(AP+B))+1.e-15;
	  TCROS=FAC*(B*AP-(B*B-AP*AP)/2.0*log(TERS));
	  CUDEL=TDIR+TCROS;
        goto ll;

	}
      SQS=XMN+abs(QZ);
    }
  else
    {
      SQS=XMN+abs(QZ);
    }
  QFR=sqrt((pow(QZ,2)-1.0)/(1.0+QZ/XMN));
  GAMH=1.0/(3.0*4.0*pi)*FDS*XMN/SQS*pow(QFR,3);

  CAC=((QZ <= 0.) ? (QZ-WRES-pow(Q,2)/(2.0*XMD))/XKF : (QZ-WRES+CYI*GAMH-pow(Q,2)/(2.0*XMD))/XKF);
  CAPC=((QZ <= 0.) ? (-QZ-WRES+CYI*GAMH-pow(Q,2)/(2.0*XMD))/XKF : (-QZ-WRES-pow(Q,2)/(2.0*XMD))/XKF);

  CDIR=(((abs(B/CAC)-0.10) <= 0. ) ? 4.0/9.0*FDS*RO/(CAC*XKF) :
	FAC*(B*CAC+(pow(B,2)-pow(CAC,2))/2.0*log((CAC+B)/(CAC-B))));

  CCROS=( (abs(B/CAPC)-0.10 <= 0. ) ? 4.0/9.0*FDS*RO/(CAPC*XKF) :
	  FAC*(B*CAPC-(pow(B,2)-pow(CAPC,2))/2.0*log((CAPC-B)/(CAPC+B))) );
  CUDEL=CDIR+CCROS;
ll:
  CUFUN=CUNUC + CUDEL;
}



double SpectralFunction::UIM(double W00,double q,double K1,double K2,const double DM)
{
  double W,L2,NUM,UIM1;

  W=abs(W00);
  L2=pow((DM*W/q-q/2),2);
  NUM=K2*K2-2*DM*W;
  if((L2 <= NUM) && (NUM <= pow(K1,2))){
    UIM1=-DM/(2*pi*q)*(K1*K1 - K2*K2 +2*DM*W);}
  else if((NUM <= L2) && (L2 < pow(K1,2))){
    UIM1=-DM/(2*pi*q)*(K1*K1 - L2);}
  else{
    UIM1=0.0;
  }
  return UIM1;
}


double SpectralFunction::DIMASM(double DK0,double DK)
{
  double DID[100000],DIFD[100000];
  using namespace Densidad;
  double DTDTDT,DI;

  if(DK0 >= 0.0){
    DTDTDT=-DK+sqrt(2.0*DM*DK0);}
  else{
    DTDTDT=0.0;
  }
  double add=-DMDF+DK;
  double aux=std::max(0.0,add);
  double D1 = std::max(aux,DTDTDT);
  double D2 = DK+DMDF;
  STGS48(D1,D2,DID);
  for(int I=1;I<=48; ++I){
    DI=DID[I];
    DIFD[I]=DI*DI*DFINM(DI, DK0, DK);
  }
//  double resu=reGS48(D1,D2,DIFD);
  //  cout<<resu<<endl;
  return reGS48(D1,D2,DIFD);
}


double SpectralFunction::DFINM(double DQ,double DK0,double DK)//hhbbi,hbb??
{
  double DIX[100000],DIFX[100000];
  double DGPRIM=0.70;
  using namespace Densidad;

  double Dm2=DM*DM;
  double Dk2=DK*DK;

  double DS=pow((2.0*DM+Dk2/2.0/DM+3.0/5.0*DMDF*DMDF/2.0/DM),2)-Dk2;

  double dlamb=1690000.0*hbbi*hbbi ;
  double ddff4=pow((dlamb+Dk2/2.0)/(dlamb+DM*DK0),4);

  double DI, DSLIM, DSNN, DENER, DSNN0, D0PROP, DCRO, DFRO2, DMESON, DVT,RRR;
  std::complex<double> CRES,CRESUL;
  DSLIM=4.0*Dm2;
  if((DS <= DSLIM ) && (DS > 0.0)){
    DENER=sqrt(DSLIM)*hbb;
    SELA(DENER,DSNN0);
    DSNN=DSNN0;
  }
  else if(DS < 0.0){
    DSNN=0.0;
  }
  else{
    DENER=sqrt(DS)*hbb;
    SELA(DENER,DSNN);
  }
  //  double DI;
  double DMPI=139.35*hbbi;
  double mdr=1400.0*hbbi;
  double mdrpi=770.0*hbbi;
  double mdr2=pow(mdr,2);
  double mdrp2=pow(mdrpi,2);
  double DCONS=-1.0/pi*DSNN*ddff4/Dm2;

  double almit=-.50*(Dk2+DQ*DQ-DMDF*DMDF)/(DK*DQ);
  double D1=std::min(1.0,almit);
  double d2limit=-DM*((Dk2+DQ*DQ)/(2.0*DM)-DK0)/(DK*DQ);
  double D2=std::max(-1.0,d2limit);
  //  cout<<D1<<" D2= "<<D2<<endl;
  STGS48(D2,D1,DIX);
  for(int I=1;I<=48; ++I)
    {
      DI=DIX[I];
      D0PROP=DK0-(Dk2+DQ*DQ)/(2.0*DM)-DQ*DK*DI/DM;
      DCRO=3.94;
      //double deno= mdr2 - D0PROP*D0PROP+DQ*DQ;
      //double nemo=mdr2-mdrp2;
//      cout<<deno<<" nemo= "<<nemo<<endl;
      double  DFaux=(mdr2-mdrp2)/(mdr2 - D0PROP*D0PROP+DQ*DQ);
      DFRO2=DFaux*DFaux;

      DMESON=1.0/(D0PROP*D0PROP - DQ*DQ - 770.0*770.0*hbbi*hbbi);
      double acons=0.080*4.0*pi/(DMPI*DMPI);
      double anum=DQ*DQ*DMESON*DFRO2*DCRO+DGPRIM;
      //      double anum=DQ*DQ*DMESON*DFRO2;
      DVT=acons*anum;
//            cout<<" DVT= "<<DVT<<endl;
      double ULINDaux1=(DK0-Dk2/(2.0*DM)-DQ*DQ/(2.0*DM)-DQ*DK*DI/DM)/DMPI;
      double ULINDaux2=DQ/DMPI;
      double ULINDaux3=DMDF/DMPI;
//       cout<<"a1= "<<ULINDaux1<<" a2= "<<ULINDaux2<<" a3= "<<ULINDaux3<<endl;
      ULIND(ULINDaux1,ULINDaux2,ULINDaux3,CRESUL);
//      cout<<CRESUL<<endl;
      CRES=DMPI*DMPI*CRESUL;
      double UIMaux=DK0-Dk2/(2.*DM)-DQ*DQ/(2.*DM)-DQ*DK*DI/DM;
      RRR=UIM(UIMaux,DQ,DMDF,DMDF,DM)/(std::abs(1.0-DVT*CRES)*std::abs(1.0-DVT*CRES));
      DIFX[I]=DCONS*RRR;
    }
  reGS48(D2,D1,DIFX);

  return reGS48(D2,D1,DIFX);
}

//     ELASTIC CROSS-SECTIONS
//     FIT BY J.VANDERMEULEN
//     PROMEDIO DE LA SECCION EFICAZ PARA PP(O NN) Y PN .5*(PP+PN)
//     E=CM ENERGY (IN MEV):SEL IN FM+2(PARA ELLO MULTIPLICO POR .1)
void SpectralFunction::SELA(double& E,double& SEL)
{
  double PLAB=E*sqrt(E*E-3.52e6)/1876.6;
  double P1=0.0010*PLAB;
  if(PLAB < 800.0){
    SEL=(23.50+1000.0*pow((P1-0.70),4)+33.0+196.0*sqrt(pow(abs(P1-0.950),5)))*.050;
    //  3rd July 2019 Rafi
    return;
  }
  else if(PLAB > 800.0 && PLAB < 2000.0) {
    SEL=(1250.0/(50.0+P1)-4.0*pow((P1-1.30),2)+31.0/sqrt(P1))*.050;
    //  3rd July 2019 Rafi
    return;
  }
  else{
    SEL=77.0/(P1+1.50)*.10;
    //  3rd July 2019 Rafi
    return;
  }
}

double SpectralFunction::reasi_int(double a,double b)
{
     // cout<<"a= "<<a<<" b= "<<b<<endl;
//  using namespace forcrho;
  double fd [15][200][200];
  int i,j,k;
  double  drorel,res1;
  double daux;
  static int ipri=0;
  //   cout<<" a= "<<a<<" b= "<<b<<endl;
  using namespace Densidad;
    if(ipri != 13233)
    {
      ifstream myfile;
      // TODO: make the data file to use here configurable
      // Right now it will unconditionally open $GENIE/data/dissf/realsigma.txt
      std::string data_file_path = std::getenv( "GENIE" );
      data_file_path += "/data/dissf/realsigma.txt";
      myfile.open( data_file_path );
      // TODO: add check that file exists and is readable
      // TODO: add input error handling
      int k1,k2,k3;
      double aux;
      for (int ii=1; ii <=120000 ; ++ii)
	{
	  myfile >> k1 >> k2 >> k3 >> aux ;

	  fd[k1][k2][k3]=aux;
	}

      myfile.close();

      ipri=13233;
    }

  double fac=DM/sqrt(b*b+DM*DM);
  drorel=(dro -0.0010)/0.170;
  double xidro= drorel*10.0+1.0;
  i= int(xidro);
  double s=xidro-i;
  if (i==0){
    i=1;
    s=0.0;
  }
  // energy
  double deltae=5.01705380;
  double eini=-1.013546210*4.0;
  double xiene=99.0*(a-eini)/deltae+1.0;
  j=int(xiene);
  double t=xiene-j;
  // momentum
  double deltam=5.017053730;
  double ximom=(99.0*b)/deltam+1.0;
  k=int(ximom);
  double u=ximom-k;


  if(i > 11 || j > 99 || k > 99){
    cout<<" sorry!! out of range! hashasdioasd"<<endl;
    return 0;
  }

  if(i < 1 || j < 1 || k < 1){
    cout<<" sorry! out of range!"<<endl;
    return 0;
  }
  daux = fd[i][j][k]*(1.0 - s)*(1.0 - t)*(1.0 - u) + fd[i + 1][j][k]*s*(1.0 - t)*(1.0 - u) + fd[i] [j + 1] [k]*(1.0 - s)*t*(1.0 - u) + fd[i] [j] [k + 1]*(1.0 - s)*(1.0 - t)*u + fd[i + 1] [j + 1] [k]*s*t*(1.0 - u) + fd[i + 1] [j] [k + 1]*s*(1.0 - t)*u + fd[i] [j + 1] [k + 1]*(1.0 - s)*t*u + fd[i + 1] [j + 1][ k + 1]*s*t*u;

//  double crho=22.00;
 if(isIso){
  res1=daux*fac+crhoIso*drorel/197.33;}
else{
//std::cout<<"me here for non Iso"<<std::endl;
res1=daux*fac+crhoNonIso*drorel/197.33;}
  return   res1;

}



double SpectralFunction::findzero(double DK, double drho, double DMU)
{
  double fd,fu,fc,xc;
  double xlimd=-1.0;
  double xlimu=DMU;
  fd=ff(xlimd,DK,drho);
  fu=ff(xlimu,DK,drho);

 ll:xc=xlimd-fd*(xlimu-xlimd)/(fu-fd);
  fc=ff(xc,DK,drho);
  double apow=pow(10,-6);
  if(abs(fc) < apow){
    return xc;
  }
  else if(fc*fd < 0.0){
    fd=fc;
    xlimd=xc; }
  else{
    fu=fc;
    xlimu=xc;
  }
  goto ll;
}




double SpectralFunction::ff(double DK0,double DK,double drho)
{
  using namespace Referd;
  double arg1=DK0-drefer;
  return  DK0-DK*DK/2.0/DM-reasi_int(arg1,DK);
}


double SpectralFunction::Deriv(double DK0,double DK,double drho)
{
  using namespace Referd;
  double dh=pow(10,-5);
  double aux1=DK0+dh-drefer;
  double aux2=DK0-drefer;

//  cout<<"Deriv "<<DMDF<<"  "<<dro<<" "<<drox<<" "<<drefer<<endl;

  double ff1=(reasi_int(aux1,DK)-reasi_int(aux2,DK))/dh;
  return sqrt((1.0-ff1)*(1.0-ff1));
}

//integration
void SpectralFunction::DSG20R(double& A, double& B, int& N,  double X[2000], int& NP)
{
  double Y[11] = {0.0,0.9931285991,0.9639719272,0.9122344282,0.8391169718,
		  0.7463319064,0.6360536807,0.5108670019,0.3737060887,0.2277858511,0.0765265211};
  NP = 20*N;
  //cout<<"NP = "<<NP<<endl;
  double XN = N;
  double DINT = B - A;
  DINT = DINT/XN;
  double DELT = DINT*0.5;
  double ORIG=A-DELT;
  double I1=-20;

  for(int i =1; i<= N; ++i)
    {
      ORIG=ORIG+DINT;
      double DORIG=ORIG+ORIG;
      I1=I1+20;
      double I2=I1+21;

      for(int j =1; j<= 10; ++j)
	{
	  int J1=I1+j;
	  int J2=I2-j;
	  X[J1]=ORIG-DELT*Y[j];
	  X[J2]=DORIG-X[J1];
	}
    }
}

double SpectralFunction::DRG20R(double A, double B, int N,double CF[2000] )
{
  double W[11] = {0.0,0.0176140071,0.0406014298,0.0626720483,0.0832767415,0.1019301198,
		  0.1181945319,0.1316886384,0.1420961093,0.1491729864,0.1527533871};
  double CR = 0.0;

  double I1=-20;

  for(int i =1; i<= N; ++i)
    {
      I1=I1+20;
      double I2=I1+21;
      for(int j =1; j<= 10; ++j)
	{
	  int J1=I1+j;
	  int J2=I2-j;
	  CR=CR+W[j]*(CF[J1]+CF[J2]);
	}
    }


  return CR*0.5*(B-A)/N;
}

void SpectralFunction::STGS48(double& X1,double& X2,double X[48])
{
  double Y[25]={0.,.998771007252426118601,.993530172266350757548,.984124583722826857745,.970591592546247250461,.952987703160430860723,.931386690706554333114,.905879136715569672822,.876572020274247885906,.843588261624393530711,.807066204029442627083,.767159032515740339254,.724034130923814654674,.677872379632663905212,.628867396776513623995,.577224726083972703818,.523160974722233033678,.466902904750958404545,.408686481990716729916,.348755886292160738160,.287362487355455576736,.224763790394689061225,.161222356068891718056,.097004699209462698930,.032380170962869362033};
  double DORIG=X1+X2;
  double ORIG=DORIG*0.5;
  double DELT=(X2-X1)*0.5;
  for(int i=1;i<=24;++i){
    int L=49-i;
    X[i]=-Y[i]*DELT+ORIG;
    X[L]=DORIG-X[i];

  }

}

double SpectralFunction::reGS48(double X1,double X2,double CF[48])
{
  double W[25]={0.,.003153346052305838633,.007327553901276262102,.011477234579234539490,.015579315722943848728,.019616160457355527814,.023570760839324379141,.027426509708356948200,.031167227832798088902,.034777222564770438893,.038241351065830706317,.041545082943464749214,.044674560856694280419,.047616658492490474826,.050359035553854474958,.052890189485193667096,.055199503699984162868,.057277292100403215705,.059114839698395635746,.060704439165893880053,.062039423159892663904,.063114192286254025657,.063924238584648186624,.06446616443590082207,.064737696812683922503};

  double CRES=0.;
  for(int i=1;i<=24;++i){
    int L=49-i;
    CRES=CRES+W[i]*(CF[i]+CF[L]);
  }
  return CRES*(X2-X1)*0.5;
}


//void startIron(double& rho0p, double& rho0N,double& rho0A)
/*
void  SpectralFunction::startNucleus(double& rho0A,double& rho0p,double& rho0n)
{
SpectralFunction spectral;

  double  xd[2000],fd[2000],fdp[2000],fdn[2000];

  double a=0.0;
  int np=1;
  int npp;
  double rlimA=12.0;
  //double rho0A;
  spectral.DSG20R(a,rlimA,np,xd,npp);
  for(int i=1;i<=npp;++i){
    double r=xd[i];
    fd[i]=r*r*(1.0/(1.0+exp((r-aden)/th)));
    fdp[i]=r*r*(1.0/(1.0+exp((r-aproton)/thp)));
    fdn[i]=r*r*(1.0/(1.0+exp((r-aneutron)/thn)));
  }
  double fdIntA=spectral.DRG20R(a,rlimA,np,fd);
  double fdIntp=spectral.DRG20R(a,rlimA,np,fdp);
  double fdIntn=spectral.DRG20R(a,rlimA,np,fdn);

  double fdresA=fdIntA*4.0*pi;
  double fdresp=fdIntp*4.0*pi;
  double fdresn=fdIntn*4.0*pi;

  onlyForRho0::rho0A =A/fdresA;
  onlyForRho0::rho0p =Z/fdresp;
  onlyForRho0::rho0n =N/fdresn;

//  return onlyForRho0::rho0A;

}
*/

double SpectralFunction::startNucleus()
{
SpectralFunction spectral;

  double  xd[2000],fd[2000];

  double a=0.0;
  int np=1;
  int npp;
  double rlimA=12.0;
  //double rho0A;
  spectral.DSG20R(a,rlimA,np,xd,npp);
  for(int i=1;i<=npp;++i){
    double r=xd[i];
    fd[i]=r*r*(1.0/(1.0+exp((r-aden)/th)));
      }
  double fdIntA=spectral.DRG20R(a,rlimA,np,fd);


  double fdresA=fdIntA*4.0*pi;

  onlyForRho0::rho0A =A/fdresA;

  return onlyForRho0::rho0A;

}

double SpectralFunction::startProton()
{
SpectralFunction spectral;

  double  xd[2000],fdp[2000];

  double a=0.0;
  int np=1;
  int npp;
  double rlimA=2.5*aproton;
  //double rho0A;
  spectral.DSG20R(a,rlimA,np,xd,npp);
  for(int i=1;i<=npp;++i){
    double r=xd[i];
       fdp[i]=r*r*(1.0/(1.0+exp((r-aproton)/thp)));

  }
   double fdIntp=spectral.DRG20R(a,rlimA,np,fdp);

   double fdresp=fdIntp*4.0*pi;
    onlyForRho0p::rho0p =Z/fdresp;

  return  onlyForRho0p::rho0p;

}

double SpectralFunction::startNeutron()
{
SpectralFunction spectral;

  double  xd[2000],fdn[2000];

  double a=0.0;
  int np=1;
  int npp;
  double rlimA=2.5*aneutron;
  //double rho0A;
  spectral.DSG20R(a,rlimA,np,xd,npp);
  for(int i=1;i<=npp;++i){
    double r=xd[i];
     fdn[i]=r*r*(1.0/(1.0+exp((r-aneutron)/thn)));
  }
   double fdIntn=spectral.DRG20R(a,rlimA,np,fdn);

   double fdresn=fdIntn*4.0*pi;

   onlyForRho0n::rho0n =N/fdresn;
  return onlyForRho0n::rho0n;

}



double  SpectralFunction::Density(double r)
{
  return onlyForRho0::rho0A/(1.0+exp((r-aden)/th));
}

double  SpectralFunction::Densityp(double r)
{
  return onlyForRho0p::rho0p/(1.0+exp((r-aproton)/thp));
}

double  SpectralFunction::Densityn(double r)
{
  return onlyForRho0n::rho0n/(1.0+exp((r-aproton)/thn));
}

///////////////////////////////////////////Pion Part////////////////
void PionCloud::forvlvt(double& p_zero, double& p_mom, double& Dzero, double& FormFactor, double&  DTILDEzero, double& FormFactorTILDE)
{
  Dzero=1.0/(p_zero*p_zero-p_mom*p_mom-xmpion*xmpion);
  double FormFacNum= xlambda*xlambda-xmpion*xmpion;
  double deno=(xlambda*xlambda)+(p_mom*p_mom);
  FormFactor=FormFacNum/deno;
DTILDEzero=1.0/((p_zero*p_zero)-(p_mom*p_mom+qc*qc)-(xmpion*xmpion));


FormFactorTILDE=FormFacNum/(xlambda*xlambda-(p_zero*p_zero-(p_mom*p_mom+qc*qc)));
//double Fdeno=xlambda*xlambda-(p_zero*p_zero-(p_mom*p_mom+qc*qc));
//cout<<"p_zero= "<<p_zero<<" p_mom= "<<p_mom<<endl;
}


void PionCloud::forvlvtrho(double& p_zero, double& p_mom, double& Dzerorho, double& FormFactorRho, double&  DTILDEzerorho, double& FormFactorTILDErho)
{
  Dzerorho=1.0/(p_zero*p_zero-p_mom*p_mom-xmrho*xmrho);
  double FormFacNumrho= xlambdarho*xlambdarho-xmrho*xmrho;
  FormFactorRho=FormFacNumrho/((xlambdarho*xlambdarho)-(0.0*p_zero*p_zero-p_mom*p_mom));

DTILDEzerorho=1.0/(p_zero*p_zero-(p_mom*p_mom+qc*qc)-xmrho*xmrho);

FormFactorTILDErho=FormFacNumrho/(xlambdarho*xlambdarho-(p_zero*p_zero-(p_mom*p_mom+qc*qc)));

}

void PionCloud::VlVtprima(double& p_zero, double& p_mom, double& Vlprima, double& Vtprima)
{
double Dzero,FormFactor,DTILDEzero,FormFactorTILDE;
double  Dzerorho,FormFactorRho,DTILDEzerorho, FormFactorTILDErho;
forvlvt( p_zero, p_mom, Dzero,FormFactor,DTILDEzero,FormFactorTILDE);
forvlvtrho( p_zero, p_mom, Dzerorho,FormFactorRho,DTILDEzerorho, FormFactorTILDErho);
//all pow
double pm2 = pow(p_mom,2);
double FFac2 = pow(FormFactor,2);
double FFacT2 = pow(FormFactorTILDE,2);
double FFacTrho2=pow(FormFactorTILDErho,2);
double  FFacRho2 = pow(FormFactorRho,2);
double D0Trho2=pow(DTILDEzerorho,2);
Vlprima=pm2*Dzero*FFac2-pm2*DTILDEzero*FFacT2-(1.0/3.0)*qc*qc*DTILDEzero*FFacT2-(2.0/3.0)*qc*qc*DTILDEzerorho*xCrho*FFacTrho2;

Vtprima=pm2*Dzerorho*xCrho*FFacRho2-(1.0/3.0)*qc*qc*DTILDEzero*FFacT2-(pm2+(2.0/3.0)*qc*qc)*D0Trho2*xCrho*FFacTrho2;

//     cout<<"Vtprima= "<<Vtprima<<endl;
}



void PionCloud::lindhard(double& q_zero,double& q_mod,double& drho,double& k_fermi,std::complex<double>& lind)
{

         const double m= DM; //mass of nucleon in fermi
          double q_mod2 =   pow(q_mod,   2);
         double rq_zero;
         const std::complex<double> iN(0,1.0);
         std::complex<double> zz;
         //double nb=pow(10.0,-36);
         double img2=imag(q_zero)*imag(q_zero);
         if(img2<1e-36)
{
          rq_zero = real(q_zero);
}


     //   z  = m /(q_mod*k_fermi) *( q_zero - q_mod**2/(2.d0*m))
   //     zp = m /(q_mod*k_fermi) *(-q_zero - q_mod**2/(2.d0*m))

 std::complex<double> z(DM/(q_mod*k_fermi)*(q_zero-q_mod2/(2.0*DM)));
 std::complex<double> zp(DM/(q_mod*k_fermi)*(-q_zero-q_mod2/(2.0*DM)));
//cout<<<<endl;

  std::complex<double> pzeta(0.0);
  if(abs(z) > 100.0){
    pzeta = 2.0/(3.0*z)+2.0/(15.0*z*z*z);
  }else if(abs(z) < pow(10.0,-2)){
    pzeta = 2.0*z-2.0/3.0*z*z*z-iN*pi/2.0*(1.0-z*z);
  }else{
    pzeta = z + (1.0-z*z) * log((z+1.0)/(z-1.0))/2.0;
  }

  std::complex<double> pzetap(0.0);
  if(abs(zp) > 100.0){
    pzetap = 2.0/(3.0*zp)+2.0/(15.0*zp*zp*zp);
//   std::cout<<pzetap<<std::endl;
  }else if(abs(zp) < pow(10.0,-2)){
    pzetap = 2.0*zp-2.0/3.0*zp*zp*zp-iN*pi/2.0*(1.0-zp*zp);

  }else{
    zz=(zp+1.0)/(zp-1.0);
    if(imag(zz)<0.0){
      zz=real(zz)+abs(imag(zz));
}
    else{zz=real(zz)+imag(zz);}

    pzetap = zp + (1.0-zp*zp)* (std::log(zz))/2.0;

 //  std::cout<<pzetap<<std::endl;
   }
 //    else{
// pzetap = zp + (1.0-zp*zp)* log(zz)/2.0;}



  //}



/*
  if(abs(zp) > 100.0){
    pzetap = 2.0/(3.0*zp)+2.0/(15.0*zp*zp*zp);
  }else if(abs(zp) < pow(10.0,-2)){
    pzetap = 2.0*zp-2.0/3.0*zp*zp*zp-iN*pi/2.0*(1.0-zp*zp);
  }else{
    pzetap = zp+ (1.0-zp*zp) * log((zp+1.0)/(zp-1.0))/2.0;
  }
 */

        lind= 3.0/2.0 * drho * m/(q_mod*k_fermi) * (pzeta +pzetap);
//        cout<<m<<" "<<pzeta<<" "<<pzetap<<endl;

  }



void PionCloud::delta_lind(double q_zero, double q_mod, double drho, double k_fermi,std::complex<double>& dlt_lind)
{
  //double q_zero = q0/fhbarc;
  //double q_mod = dq/fhbarc;
  //double k_fermi = kF/fhbarc;
  //Divide by hbarc in order to use natural units (rho is already in the correct units)

  //m = 939/197.3, md = 1232/197.3, mpi = 139/197.3
//  double m = 4.7592; DM
//  double md = 6.2433;
//  double mpi = 0.7045;
   const double m= DM;
   const double fdel_f = 2.13;
   const double wr = md-m;
    double gamma = 0;
    double gammap = 0;

  double q_zero2 =  pow(q_zero,  2);
  double q_mod2 =   pow(q_mod,   2);
  double k_fermi2 = pow(k_fermi, 2);

  double m2 =       pow(m,       2);
  double m4 =       pow(m,       4);
  double mpi2 =     pow(mpi,     2);
  double mpi4 =     pow(mpi,     4);

  double fdel_f2 =  pow(fdel_f,  2);
  double rq_zero;
  //For the current code q_zero is always real
  //If q_zero can have an imaginary part then only the real part is used
  //until z and zp are calculated
    double nb=pow(10.0,-36);
    double img2=imag(q_zero)*imag(q_zero);
    if(img2<nb)
{
          rq_zero = real(q_zero);
}

  double s = m2+q_zero2-q_mod2+
    2.0*q_zero *sqrt(m2+3.0/5.0*k_fermi2);

  if(s>pow(m+mpi,2)){
    double srot = sqrt(s);
    double qcm = sqrt(pow(s,2)+mpi4+m4-2.0*(s*mpi2+s*m2+
	      	mpi2*m2)) /(2.0*srot);
    gamma = 1.0/3.0 * 1.0/(4.0*pi) * fdel_f2*pow(qcm,3)/srot*(m+sqrt(m2+pow(qcm,2)))/mpi2;
  }
  double sp = m2+q_zero2-q_mod2-2.0*q_zero *sqrt(m2+3.0/5.0*k_fermi2);


  if(sp >pow(m+mpi,2)){
    double srotp = sqrt(sp);
    double qcmp = sqrt(pow(sp,2)+mpi4+m4-2.0*(sp*mpi2+sp*m2+
		 mpi2*m2))/(2.0*srotp);
    gammap = 1.0/3.0 * 1.0/(4.0*pi) * fdel_f2*pow(qcmp,3)/srotp*(m+sqrt(m2+pow(qcmp,2)))/mpi2;
  }
  //}//End if statement
  const std::complex<double> iNum(0,1.0);

  std::complex<double> z(md/(q_mod*k_fermi)*(q_zero-q_mod2/(2.0*md)
                         -wr +iNum*gamma/2.0));
  std::complex<double> zp(md/(q_mod*k_fermi)*(-q_zero-q_mod2/(2.0*md)
                          -wr +iNum*gammap/2.0));

  std::complex<double> pzeta(0.0);
  if(abs(z) > 50.0){
    pzeta = 2.0/(3.0*z)+2.0/(15.0*z*z*z);
  }else if(abs(z) < pow(10.0,-2)){
    pzeta = 2.0*z-2.0/3.0*z*z*z-iNum*pi/2.0*(1.0-z*z);
  }else{
    pzeta = z + (1.0-z*z) * log((z+1.0)/(z-1.0))/2.0;
  }

  std::complex<double> pzetap(0);
  if(abs(zp) > 50.0){
    pzetap = 2.0/(3.0*zp)+2.0/(15.0*zp*zp*zp);
  }else if(abs(zp) < pow(10.0,-2)){
    pzetap = 2.0*zp-2.0/3.0*zp*zp*zp-iNum*pi/2.0*(1.0-zp*zp);
  }else{
    pzetap = zp+ (1.0-zp*zp) * log((zp+1.0)/(zp-1.0))/2.0;
  }

  //Multiply by hbarc^2 to give answer in units of GeV^2
  dlt_lind=2.0/3.0 * drho * md/(q_mod*k_fermi) * (pzeta +pzetap) * fdel_f2;
//cout<<dlt_lind<<endl;

//      cout<<"z= "<<pzeta<<"  "<<pzetap<<endl;
}



 void PionCloud::ULIND2(double& QZR,double& q,double& drho,std::complex<double>& CUFUN)
{
      //complexNumber dlt_lind,lind,sum;
      std::complex<double> lind, dlt_lind;

//         using namespace Densidad;
          double q_zero=real(QZR);
         double xkf=pow((3.0/2.0*pi*pi*drho),(1.0/3.0));
        double q_zerol = q_zero;

            lindhard(q_zero,q,drho,xkf,lind);

           delta_lind(q_zerol,q,drho,xkf,dlt_lind);
           CUFUN=lind+dlt_lind;
//          sum=addCN(lind,dlt_lind);

	if(abs(imag(CUFUN))<1.e-30){
	CUFUN=real(CUFUN);
}
//cout<<lind<<" "<<dlt_lind<<endl;
}



/*************************************************
c Here we are going to define D(p)-D_0(p). We will need the imaginary
c part of this function for the contribution of the pions to the
c nuclear structure function.
***********************************************/
   std::complex<double> PionCloud::CDminusD0(double p_zero,double p_mom,double drho)
   {

	double pzr=p_zero;
	double pmom=p_mom;
	double xrho=drho;
	 std::complex<double> CUFUN;
      double  Dzero, FormFactor,DTILDEzero,FormFactorTILDE;
      double  Vlprima,Vtprima;

     forvlvt(p_zero, p_mom, Dzero, FormFactor, DTILDEzero,FormFactorTILDE);

    VlVtprima(p_zero, p_mom, Vlprima,Vtprima);

	 ULIND2(pzr,pmom,xrho,CUFUN);


     return Dzero*Dzero*((f2/xmpion*xmpion)*FormFactor*FormFactor*p_mom*p_mom*CUFUN)/(1.0-(f2/xmpion*xmpion)*Vlprima*CUFUN);

}
/*
***************************************************************
c here we are going to define Drho(p)-D_0rho(p).We will need the
c imaginary part of this function for the contribution of the rho
c to the nuclear structure function.
****************************************************************/
std::complex<double> PionCloud::CDminusD0rho(double p_zero,double p_mom,double drho)
{

	double pzr=p_zero;
	double pmom=p_mom;
	double xrho=drho;
      std::complex<double> CUFUN;

      double Dzerorho,FormFactorRho,DTILDEzerorho,FormFactorTILDErho, Vlprima,Vtprima;


    forvlvtrho(p_zero,p_mom, Dzerorho,FormFactorRho,DTILDEzerorho,FormFactorTILDErho);
	 VlVtprima(p_zero, p_mom, Vlprima,Vtprima);
	  ULIND2(pzr,pmom,xrho,CUFUN);



	return Dzerorho*Dzerorho*((f2/xmpion*xmpion)*xCrho*FormFactorRho*FormFactorRho*p_mom*p_mom*CUFUN)/(1.0-(f2/xmpion*xmpion)*Vtprima*CUFUN);

}

void PionCloud::deltaimD2(double& p_zero,double& p_mom,double& drho,double& deltaIm)
{


	double pzr=p_zero;
	double pmom=p_mom;
	double xrho=drho;

   //   cout<<"pzr= "<<pzr<<" pmom= "<<pmom<<endl;
      std::complex<double> CUFUN,CUFUN1,CDminusD1,cdeno,check;
       double  Vlprima,Vtprima;
      double  Dzero, FormFactor,DTILDEzero,FormFactorTILDE;
      forvlvt(p_zero, p_mom, Dzero, FormFactor, DTILDEzero,FormFactorTILDE);
      double D02=pow(Dzero,2);
      double xmpi2=pow(xmpion,2);
      double FormFac2=pow(FormFactor,2);
      double p_mom2=pow(p_mom,2);

	 ULIND2(pzr,pmom,xrho,CUFUN);
       VlVtprima(p_zero, p_mom, Vlprima,Vtprima);

      double aux=D02*(f2/xmpi2)*FormFac2*p_mom2;
// checked      cout<<"aux= "<<aux<<"  "<<FormFactor<<endl;
     // cout<<" Vlprima= "<<Vlprima<<" Vtprima= "<<Vtprima<<endl;
      std::complex<double> caux2=aux*CUFUN/(1.0-(f2/xmpi2)*Vlprima*CUFUN);

//      cout<<CUFUN<<endl;
       double num=pow(10.0,-20);

      ULIND2(pzr,pmom,num,CUFUN1);
        cdeno=pow((1.0-(f2/xmpi2)*Vlprima*CUFUN1),2);

       check=num*aux*CUFUN1;
	CDminusD1=caux2-drho/num*aux*CUFUN1/cdeno;
       deltaIm=imag(CDminusD1);
//       cout<<deltaIm<<endl;
    }

/*
**************************************************************
*************          FILE RICHARDSON.F         *************
**************************************************************




*************************************************************
*********   SUBROUTINE TO CALCULATE THE DERIVATIVE    *******
*********             OF A FUNCTION                   *******
*************************************************************

c     We use the Richardson method to calculate the derivative.
c     If we calculate f(x+h)-f(x-h) in powers of h
c     f(x+h)-f(x-h)=2h*f^(1)(x)+(2/3!)h^3*f^(3)(x)+(2/5!)h^5*f^(5)(x)+...
c     where f^(i)(x) indicates the derivative of order i. Calling D=f^(1)(x),
c     from this expression we can see that

c                     D=D0(h)+a2*h^2+a4*h^4+a6*h^6+...

c     Changing h --> h/2 we can remove the a2 coefficient and put

c                     D=D1(h)+b4*h^4+b6*h^6+...

c     where D1(h)=(4/3)*D0(h/2)-(1/3)*D0(h). Repeating this process, we
c     arrive to a recursive formula

c              D(n,k)=(4^k/(4^k-1))*D(n,k-1)-(1/((4^k-1))*D(n-1,k-1))

c     with k=1,2,...,m and n=k,k+1,...,m (m is the number of extrapolations).

c     The method consists of obtaining D(l,0)=D0(h/2^l), 0<= l <= n, and
c     then D(n,k).


*************************************************
***        SUBROUTINE THAT CALCULATE D0       ***
***          FOR THE FIRST DERIVATIVE         ***
*************************************************
c Variables de entrada
c     'h' es el paso que tomamos inicialmente
c     'rho' es el punto (densidad) donde se desea calcular la derivada
c	'p_zero' es la energía de la que depende la funcion
c	'p_mom' es el momento del que depende la funcion


c Variables de salida
c     'D0' es D0(h), o sea, la primera aproximacion a la derivada. A partir
c de este D0(h), D0(h/2),D0(h/4)... se calculan las extrapolaciones de Richardson
*/
      void PionCloud::Difirstandrho(double& h,double& drho,double& p_zero,double& p_mom,double& D0,double& D0rho)
{


      double xplush=drho+h;
      double xminush=drho-h;


//c Calcularemos la función en x+h y en x-h y la pondremos por pantalla
//c para controlar el proceso

//c      D0=(funcion(xplush)-funcion(xminush))/(2.d0*h)
      D0=(1.0/(2.0*h))*(imag(CDminusD0(p_zero,p_mom,xplush))-imag(CDminusD0(p_zero,p_mom,xminush)));
      D0rho=(1.0/(2.0*h))*(imag(CDminusD0rho(p_zero,p_mom,xplush))-imag(CDminusD0rho(p_zero,p_mom,xminush)));

      }

/*
*************************************************
***     SUBROUTINE THAT CALCULATE D(n,k)      ***
***         FOR THE FIRST DERIVATIVE          ***
*************************************************
c Variables de entrada
c     'm' es el número de extrapolaciones de Richardson
c     'N' es un entero mayor que 'm' que permite dimensionar la matriz
c adecuadamente en el programa principal y en la subrutina
c     'h' es el paso que tomamos inicialmente
c     'rho' es el punto donde se desea calcular la derivada
c	'p_zero' es la energía de la que depende la funcion
c	'p_mom' es el momento del que depende la funcion


c Variables de salida
c     'D' es la matriz donde estan contenidas las extrapolaciones de Richardson
c Los elementos de la diagonal de dicha matriz son las extrapolaciones de Richardson
c D0(h),D1(h),D2(h),...,Dm(h)
*/
      void PionCloud::RCH1(int& m,double& h,double& drho,double& p_zero,double& p_mom,double D[100][100])
{
       //D[N][N];
      double hi,D0,D0rho;
      for(int i=0;i<=m;++i){
         hi=h/(pow(2.0,i));
         Difirstandrho(hi,drho,p_zero,p_mom,D0,D0rho);
         D[i][0]=D0;
       }

      for(int k=1;k<=m;++k){
        for(int j=k;j<=m;++j){
            D[j][k]=(pow(4.0,k)/(pow(4.0,k)-1.0))*D[j][k-1]-(1.0/(pow(4.0,k)-1.0))* D[j-1][k-1];
         }
      }
}

/**************************************************************
**************************************************************
**************************************************************
*************************************************
***     SUBROUTINE THAT CALCULATE D(n,k)      ***
***         FOR THE FIRST DERIVATIVE          ***
*************************************************
c Variables de entrada
c     'm' es el número de extrapolaciones de Richardson
c     'N' es un entero mayor que 'm' que permite dimensionar la matriz
c adecuadamente en el programa principal y en la subrutina
c     'h' es el paso que tomamos inicialmente
c     'rho' es el punto donde se desea calcular la derivada
c	'p_zero' es la energía de la que depende la funcion
c	'p_mom' es el momento del que depende la funcion


c Variables de salida
c     'D' es la matriz donde estan contenidas las extrapolaciones de Richardson
c Los elementos de la diagonal de dicha matriz son las extrapolaciones de Richardson
c D0(h),D1(h),D2(h),...,Dm(h)
*/
      void PionCloud::RCH1rho(int& m,double& h,double& drho,double& p_zero,double& p_mom,double D[100][100])
{
  //    double D[N][N];
      double hi,D0,D0rho;
      for(int i=0;i<=m;++i){
         hi=h/(pow(2.0,i));
         Difirstandrho(hi,drho,p_zero,p_mom,D0,D0rho);
         D[i][0]=D0;
       }

      for(int k=1;k<=m;++k){
        for(int j=k;j<=m;++j){
            D[j][k]=(pow(4.0,k)/(pow(4.0,k)-1.0))*D[j][k-1]-(1.0/(pow(4.0,k)-1.0))* D[j-1][k-1];
         }
      }
}



/********************************************************************
********************************************************************
********************************************************************
******************************************************
******************************************************
******************************************************
c This subroutine is for calculating the derivative of
c the imaginary part of D(p)-D0(p) with respect to density
c at density zero. It is also useful for any extrapolation
*/
/*
*********************************************************
*********************************************************
c Input variables

c     rho1 is the density of the first point
c     rho2 is the density (abscisa) of the second point
c     f1 is the value of the function at rho1
c     f2 is the value of the function at rho2
c     rho is the point where we want to know the value of
c              the extrapolation, in our case at rho=0.d0


c Output variable

c     f is the value of the function in the point of the
c              extrapolation, point rho
**********************************************************
**********************************************************
*/

//      void extrapol(double& drho1,double& drho2,double& f1,double& f2,double& drho,double& f)
  // {
    //  f=f1+(f2-f1)/(drho2-drho1)*(drho-drho1);
  // }

/*
******************************************************************
******************************************************************
******************************************************************
c Here we are going to define the subroutine deltaimD in terms of
c the imaginary part of CDminusD0 and the derivative with respect
c to density (at density=0.d0) of the same imaginary part.

****************************************************************
****************************************************************
c Input variables
c p_zero is the energy
c p_mom is the momentum
c rho is the density

c Output variables
c deltaIm is the output
****************************************************************
****************************************************************
*/

	void PionCloud::deltaimD(double& p_zero,double& p_mom,double& drho,double& deltaIm)
{
//	int N=100;
	double D[100][100]={};
	double D1[100][100]={};
/*
c First of all, we are going to evaluate the derivative with
c an extrapolation. For it, we are going to need to evaluate
c the derivative in two points near zero (rho1=0.001d0 and
c rho2=0.002d0) and then extrapolate to rho0=0.d0

c D is a matrix where we are going to put the values of the
c derivatives calculated with 3 Richardson extrapolations (3
c is a value for which we have seen that everything is OK). We
c will only need the value D(m,m) with m=3 (the number of
c Richardson extrapolations)
*/
	int m=3;
// Initial step for the derivative h=0.0001d0
	double h=0.00010;
// Values for the density before the extrapolation.
// These values must be close to zero and close between them
	double drho1=0.001;
	double drho2=0.002;
  //    double der1,der2;
//c Now we call to the subroutine that does the calculation of
//c the first derivative at rho=rho1

	RCH1(m,h,drho1,p_zero,p_mom,D);
//	 D=der1;


//c Call to the subroutine for rho2

       RCH1(m,h,drho2,p_zero,p_mom,D1);
   //     D1=der2;

//c Point of the extrapolation
 //     double drho0=0.0;
//c Call to the subroutine of the extrapolation to know
//c the value of the derivative at rho0=0.d0

       double der;
  //     extrapol(drho1,drho2,D,D1,drho0,der);
          der=D[m][m]+(D1[m][m]-D[m][m])/(drho2-drho1)*(drho-drho1);

//c Now we can define the output variable

	deltaIm=imag(CDminusD0(p_zero,p_mom,drho))-drho*der;
}

/*
**************************************************************
**************************************************************
**************************************************************
c Here we are going to define the subroutine deltaImDrho2
c This subroutine calculates deltaImDrho(p) making use of the
c low density theorem for calculating the derivative of the imaginary
c part of Drho(p)-D_0rho(p)
*/
void PionCloud::deltaImDrho2(double& p_zero,double& p_mom,double& drho,double& deltaImrho)
{
      std::complex<double>CUfunc;
//c The first term of deltaImrho is the imaginary part
//c of Drho(p)-D_0rho(p)

	double firstterm=imag(CDminusD0rho(p_zero,p_mom,drho));

//c Now the second term with the derivative with respect to rho of
//c Drho(p)-D_0rho(p)
      double Dzerorho,FormFactorRho,DTILDEzerorho,FormFactorTILDErho,Vlprima,Vtprima;

     forvlvtrho(p_zero,p_mom,Dzerorho,FormFactorRho,DTILDEzerorho,FormFactorTILDErho);
      double D0rho2=pow(Dzerorho,2);

      double xmpi2=pow(xmpion,2);
      double FormFacRho2=pow(FormFactorRho,2);
      double p_mom2=pow(p_mom,2);

	double aux=D0rho2*(f2/xmpi2)*xCrho*FormFacRho2*p_mom2;

//c Now we call to ULIND2 to obtain U(rho) at rho-->0
	double pzr=p_zero;
	double pmom=p_mom;
	double xrho=pow(10,-20);
	ULIND2(pzr,pmom,xrho,CUfunc);
// Now CUfunc is the Lindhard function at density 1.d-20
       VlVtprima(p_zero, p_mom, Vlprima,Vtprima);

	std::complex<double> cderiv=(aux*(CUfunc/xrho))/((1.0-(f2/xmpi2)*Vtprima*CUfunc)*(1.0-(f2/xmpi2)*Vtprima*CUfunc));

//c Now the second term, which is rho times the imaginary part
//c of the derivative
	double secondterm=drho*imag(cderiv);
	  deltaImrho=firstterm-secondterm;

}

/*
******************************************************************
******************************************************************
******************************************************************
c Here we are going to define the subroutine deltaimDrho in terms of
c the imaginary part of CDminusD0rho and the derivative with respect
c to density (at density=0.d0) of the same imaginary part.

****************************************************************
****************************************************************
c Input variables
c p_zero is the energy
c p_mom is the momentum
c rho is the density

c Output variables
c deltaImrho is the output
****************************************************************
****************************************************************
*/
/*c First of all, we are going to evaluate the derivative with
c an extrapolation. For it, we are going to need to evaluate
c the derivative in two points near zero (rho1=0.001d0 and
c rho2=0.002d0) and then extrapolate to rho0=0.d0

c D is a matrix where we are going to put the values of the
c derivatives calculated with 3 Richardson extrapolations (3
c is a value for which we have seen that everything is OK). We
c will only need the value D(m,m) with m=3 (the number of
c Richardson extrapolations)
*/

	void PionCloud::deltaimDrho(double& p_zero,double& p_mom,double& drho,double& deltaImrho)
{
	double D[100][100];
	double D1[100][100];


       int m=3;
//c Initial step for the derivative h=0.0001d0
	double h=0.0001;
//c Values for the density before the extrapolation.
//c These values must be close to zero and close between them
     double	drho1=0.001;
	double drho2=0.002;

//c Now we call to the subroutine that does the calculation of
//c the first derivative at rho=rho1

	 RCH1rho(m,h,drho1,p_zero,p_mom,D);
//	der1=D(m,m)


//c Call to the subroutine for rho2

       RCH1rho(m,h,drho2,p_zero,p_mom,D1);
//      der2=D(m,m)

//c Point of the extrapolation
//      rho0=0.d0
//c Call to the subroutine of the extrapolation to know
//c the value of the derivative at rho0=0.d0

//      call extrapol(rho1,rho2,der1,der2,rho0,der)
	 double der=D[m][m]+(D1[m][m]-D[m][m])/(drho2-drho1)*(drho-drho1);
//c Now we can define the output variable

	deltaImrho=imag(CDminusD0rho(p_zero,p_mom,drho))-drho*der;

}


/**************************************************************
*************************************************************
*************************************************************
c Parameterizations of Pionic Parton Distribution Functions from a work
c by Glueck, Reya and Vogt, published in Z.Phys.C-Particles and Fields 53, 651-655 (1992)

c For each distribution, you must give x (Bjorken variable) and Q2 (scale in GeV**2)
c the ranges are:
c xmu2LO <= Q2 <= 10**8 GeV**2
c 10**(-5) <= x < 1
c These are the parameterizations for LO parton distributions
c Always x times the distribution is returned

c VALENCE DISTRIBUTION
*/
	double PionCloud::xvalencepion(double x,double Q2)
{
//c Parameters for the parameterization
	double xmu2LO=0.25;
//c Definition of s
	double s=log(log(Q2/(0.2320*0.2320))/log(xmu2LO/(0.2320*0.2320)));
//      cout<<s<<"  "<<log(Q2/(0.2320*0.2320))<<"  "<<Q2<<endl;
//c Definition of parameters in terms of s
	double c=0.5190+0.1800*s-0.0110*s*s;
	double a=0.4990-0.0270*s;
	double Aprima=0.3810-0.4190*s;
      double D=0.3670+0.5630*s;
      double  xvalenpion;
	if ((x > 0.0) && (x < 1.0)){
	 xvalenpion=c*pow(x,a)*(1.0+Aprima*sqrt(x))*pow((1.0-x),D);
       }
	else
       {
		xvalenpion=0.0;
	}

return xvalenpion;
}

//c SEA QUARK DISTRIBUTIONS


	double PionCloud::xseapion(double x,double Q2)
{
	double xmu2LO=0.25,xsea;
//c Definition of s
   double s=log(log(Q2/(0.2320*0.2320))/(log(xmu2LO/(0.2320*0.2320))));
//   cout<<s<<endl;
//c Definition of parameters in terms of s
	double sqbar=0.0;
	double alpha=0.55;
	double beta=0.56;
	double a=2.5380-0.7630*s;
	double Aprima=-0.7480;
	double B=0.3130+0.9350*s;
      double D=3.3590;
	double E=4.4330+1.3010*s;
	double Eprima=9.300-0.8870*s;

	if ((x > 0.0) && (x < 1.0)) {
	xsea= (pow((s-sqbar),alpha))/(pow((log(1.0/x)),a))*(1.0+Aprima*sqrt(x)+B*x)*pow((1.0-x),D)*exp(-E+sqrt(Eprima*pow(s,beta)*log(1.0/x)));
}
	else {xsea= 0.0;}
return xsea;
}

/*********************************************************************
c function F2pion at LO in terms of xvalencepion and xseapion
c remember that F2pion=2*x*v_pion +6*x*qbar_pion
*/
	double PionCloud::F2pion(double x,double Q2)
{
       double f2pi;
	if ((x > 0.0) && (x < 1.0)){
	   f2pi= 2.0*xvalencepion(x,Q2)+6.0*xseapion(x,Q2);
	   }
	else{f2pi= 0.0;}
//     cout<<xvalencepion(x,Q2)<<"   "<<xseapion(x,Q2)<<endl;
return f2pi;
}

void ShadowingEffect::c2A_c3A_ct1A(double& x,double& q2, std::complex<double>& C2A,std::complex<double> & C3A, std::complex<double>& CT1A)
{
  using namespace SigmaShad;
  const std::complex<double> iN(0,1.0);
//    RA is given on Page-184/Appendix/NPA/KP
  double RA=sqrt(5.0/3.0)*rms;
  double RA3=pow(RA,3);
  double rho0=A/(4.0*pi/3.0*RA3);
    const double sigma0=27.0;
const double sigma1=0.0;
const double Q02=1.39;

   std::complex<double> ca0T,cy,cy3,cy4,cy5,cphi2,cphi3;
// Complex amplitude a0T
//    ca0T is the forward scattering amplitude on Page 8/PRD76/KP
    double sigT=(sigma1+(sigma0-sigma1)/(1.0+q2/Q02))*1.e-1;
    SigmaShad::sigmaT=sigT;
      ca0T=(iN+alpha0T)*sigmaT/2.0;
//c Real variable kL
//c    kL is given on Page 144 of KP NPA765
      double xkL=DM*x*(1.0+xmeff*xmeff/q2);

//c Complex adimensional variable y
//c    cy is given on Page-184/Appendix-B/NPA/KP
      cy=2.0*iN*(rho0*ca0T-xkL)*RA;
      cy3=pow(cy,3);
      cy4=pow(cy,4);
       cy5=pow(cy,5);
     cphi2=(6.0-3.0*cy*cy-2.0*cy3+6.0*(cy-1.0)*exp(cy))/cy4;

//c Definition of the function
//c    c2a is given on Page-184/Appendix-B/NPA/KP
      C2A=A*rho0*RA*cphi2;

//      c This is Eq.B.2b of KP NPA765
      cphi3=12.0*(-4.0+cy*cy+cy3/3.0+(2.0-cy)*(2.0-cy)*exp(cy))/cy5;

//c Definition of the function
//c    c3a is given on Page-184/Appendix-B/NPA/KP
      C3A=A*(rho0*RA)*(rho0*RA)*cphi3;

    CT1A=2.0*iN*ca0T*C2A-ca0T*ca0T*C3A;
//cout<<"C2A ="<<C2A<<endl;
}

double ShadowingEffect::deltaR2(double x,double q2,std::string setname, int imem,int iset,const PDF* pdf)
{
     using namespace SigmaShad;
//c Q2 in GeV**2
      const std::complex<double> iN(0,1.0);
       std::complex<double> C2A,C3A,CT1A,ca0t;
       //double qgev=sqrt(q2);
        c2A_c3A_ct1A(x,q2,C2A,C3A,CT1A);
        ca0t=(iN + alpha0T)*SigmaShad::sigmaT/2.0;

         double sigma0T=2.0*imag(ca0t);
       double dRtransF1=sigma0T*real((iN+alpha0T)*(iN+alpha0T)*C2A)/(2.0*A);
//       cout<<ca0t<<"  "<<sigma0T<<"  "<<dRtransF1<<endl;
	return dRtransF1;
}

double ShadowingEffect::deltaR2_F2(double x,double q2,std::string setname, int imem,int iset,const PDF* pdf)
{
//c Q2 in GeV**2
//      NucleonNLO Nucleon;
     ForShadowinfF1 shadowf1;
     using namespace SigmaShad;
      const std::complex<double> iN(0,1.0);
       std::complex<double> C2A,C3A,CT1A;
       double ratio,f1f,Flf;
       double qgev=sqrt(q2);
      if(iset==3){
	ratio=4.0*Mn*Mn*x*x/(qgev*qgev);}
        else{
       f1f=shadowf1.f1nucshad(x,q2,setname,imem,iset,pdf);
       Flf=shadowf1.Flongitudinal(x,q2,setname,imem,iset,pdf);
        ratio=Flf/(2.0*x*f1f);}
        c2A_c3A_ct1A(x,q2,C2A,C3A,CT1A);

       double dRtrans=sigmaT*real((iN+alpha0T)*(iN+alpha0T)*C2A)/(2.0*A);
       //double dR=(1.0+ratio*ratio)/(1.0+ratio)*dRtrans;
//       cout<<"ratio= "<<ratio<<"   "<<dR<<endl;
//c definition of the function
//     cout<<ratio<<"   "<<dR<<endl;
	return (1.0+ratio*ratio)/(1.0+ratio)*dRtrans;
}

//c This is formula 45b in second reference hep-ph/0703033 for
//c the shadowing correction to the F3(nu+nubar) structure function
//c Q2 in GeV**2
double ShadowingEffect::deltaR_F3(double x,double q2,std::string setname, int imem,int iset,const PDF* pdf)
{
       std::complex<double> C2A,C3A,CT1A;
      const double alpha0delta=1.15;
       const std::complex<double> iN(0,1.0);
         c2A_c3A_ct1A(x,q2,C2A,C3A,CT1A);
//c     alpha0delta is taken from Table-I/Page-12/ PRD76/KP
//c    deltaRdelta is taken from Eq.45b/Page-10/ PRD76/KP
//c       deltaRdelta=dimag((i+alpha0delta)*cT1A(x,Q2))/A
//    	dRdelta=imag((iN+alpha0delta)*CT1A(x,q2))/A
      return imag((iN+alpha0delta)*CT1A)/A;
}
