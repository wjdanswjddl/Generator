#ifndef _Param_H_
#define _Param_H_
/*#include <cmath>*/
//  aq2 is in GeV. unless specified.
// q2 is in fm. unless specified.

namespace Param
{
// At Nucleon level
//   const bool IsWcut=true;
   const bool noWCut=false;
   const bool isNLO=false;
  const bool isTMC_NLO=true;
  const bool isTMC_NLO_HT=false;
// At Nucleus level 
  const bool isIso=true;


  double const wlimit=2.0; // in GeV
  const double A=56.0;
  const double Z=26.0;
  const double N=30.0;
//For Non-ISo density parameters
//presently fot Fe-56
  const double aproton=3.971;
  const double thp=0.5935;
  const double aneutron=4.05;
  const double thn=0.5935;
// For Iso density parameters
  const double aden=4.106;
  const double th=0.519;
// Fixed parameter on the basis of binding energy of the nucleus 
  const double crhoIso=22.0;
  const double crhoNonIso=-22.0;
//spectral function Normalization
  const double Normal_Iso=57.74176;
  const double Normal_NonIso=59.07;//27.40+31.671;

  
  //end here
  const double pi=3.141592653589793;
  const double Mn=0.93826; //mass of Nucleon in Gev
  const double DM=4.75599250;//mass of nucleon in fm
  const double md =6.2433;// mass of Delta in fm 1232./197.3
  const double mpi = 0.7045; // mass of pion in fm mpi = 139./197.3 
  const double hbb=197.33; //for MeV to fm
  const double hbbi=0.0050676; //1/hbb
  const double GeVtfm=0.197330; // for convertin in to fm or Gev
  const double esquared=1.0;// charge square =1 for weak interactions 
  const double Cf=4.0/3.0;
  const int nf=4; 
/// For Pion all in fm
const double f2=1.00531;
const double qc=3.95277;
const double xCrho=3.94;
const double xmrho=3.90716;
const double xmpion=0.704404;
const double xlambda=5.06765; //in fm
const double xlambdarho=5.06765; //in fm

//Parameters for Shadowing

const double alpha0T=-0.20;
const double xmeff=0.8451;
const double rms=3.721; //for Fe
}


#endif
