//____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Jesús González-Rosa <jgrosa \at us.es>
 University of Seville

 Steven Gardiner <gardiner \at fnal.gov>
 Fermi National Accelerator Laboratory
*/
//____________________________________________________________________________

// Standard library includes
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

// ROOT includes
#include "TSystem.h"

// GENIE includes
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/Resonance/XSection/SuSAv2InelPXSec.h"

namespace {

  constexpr size_t VECTOR_SIZE = 1002000u;

}

//____________________________________________________________________________
genie::SuSAv2InelPXSec::SuSAv2InelPXSec()
  : XSecAlgorithmI( "genie::SuSAv2InelPXSec" )
{
}
//____________________________________________________________________________
genie::SuSAv2InelPXSec::SuSAv2InelPXSec( std::string config )
  : XSecAlgorithmI( "genie::SuSAv2InelPXSec", config )
{
}
//____________________________________________________________________________
genie::SuSAv2InelPXSec::~SuSAv2InelPXSec()
{
}
//____________________________________________________________________________
double genie::SuSAv2InelPXSec::XSec( const genie::Interaction* interaction,
  genie::KinePhaseSpace_t kps ) const
{
  //if ( !this->ValidProcess    (interaction) ) return 0.;
  //if ( !this->ValidKinematics (interaction) ) return 0.;

  const genie::InitialState& init_state = interaction->InitState();
  const genie::ProcessInfo& proc_info = interaction->ProcInfo();
  const genie::Target& target = init_state.Tgt();

  // Get kinematical parameters
  const genie::Kinematics& kinematics = interaction->Kine();

  // Hadronic invariant mass
  double W = kinematics.W();

  // Final lepton kinetic energy
  double Tl = kinematics.GetKV( genie::kKVTl );

  // Final lepton scattering cosine
  double ctl = kinematics.GetKV( genie::kKVctl );

  // Get target nucleus information
  int Anumber = target.A();
  int Znumber = target.Z();

  // Check whether this is an EM or weak (CC only for now) process
  bool is_em = interaction->ProcInfo().IsEM();

  // Choose the appropriate minimum Q^2 value based on the interaction
  // mode (this is important for EM interactions since the differential
  // cross section blows up as Q^2 --> 0)
  double Q2min = genie::controls::kMinQ2Limit; // CC/NC limit
  if ( is_em ) Q2min = genie::utils::kinematics::electromagnetic::kMinQ2Limit;

  // Get a sign dependent on whether we're working with neutinos or
  // antineutrinos
  // TODO: revisit for electrons
  int  probepdgc = init_state.ProbePdg();

  int inu = 0;
  if ( genie::pdg::IsNeutrino(probepdgc) ) inu = 1;
  else if ( genie::pdg::IsAntiNeutrino(probepdgc) ) inu = -1;
  else {
    // TODO: revisit this!
    LOG( "SuSAv2InelPXSec", pFATAL ) << "Bad inu value for probe!";
    std::exit( 1 );
  }

  // Final-state lepton mass
  double xmil = interaction->FSPrimLepton()->Mass();

  // Probe energy
  double ener1 = init_state.ProbeE( genie::kRfLab );

  // Energy transfer
  double w = ener1 - xmil - Tl;

  // Final-state lepton momentum
  double pl2 = std::pow( xmil + Tl, 2 ) - xmil*xmil;
  double pl = std::sqrt( std::max(0., pl2) );

  // Probe momentum
  TLorentzVector* p4_probe = init_state.GetProbeP4( genie::kRfLab );
  double pv = p4_probe->P();
  delete p4_probe;

  // Momentum transfer
  double q =std::sqrt(std::max(0., pv*pv + pl2 - 2.*pv*pl*ctl));

//////////////////////////////////////////////////////////////
// BEGIN CODE FROM JESUS

  // Constants
  double pi = std::acos( -1.0 );
  double GF = genie::constants::kGF; //1.16639E-5; //Cte Fernu, en GeV**-2
  double rmn = genie::constants::kNucleonMass; //(rmProton + rmNeutron)/2.0;

  // Fermi momentum
  double etaF = pF / rmn;
  double epsF = std::sqrt( 1.0 + etaF*etaF );
  double xiF = epsF - 1.0;

  //KINEMATICS
  double ener2 = ener1 - w; // Energy of the lepton
  double xkp = std::sqrt( ener2*ener2 - xmil*xmil ); // Momentum of the lepton
  double Q2 = w*w - q*q; //Four-momentum
  double vo = 4*ener1*ener2 + Q2;
  
  //std::cout<<" Q2 " << Q2 << std::endl;

  // Enforce the global GENIE minimum limit on Q^2. Note that the variable
  // above has negative values (it's really q2). We therefore add a minus
  // sign in the check.
  // TODO: change this to avoid confusion
  if ( -Q2 < Q2min ) {
    //LOG( "SuSAv2Inel", pWARN ) << "Q2 = " << -Q2 << ", Q2min = " << Q2min;
    return 0.;
  }

  // Scaling and Adimensional variables
  double xk = q / ( 2 * rmn );
  double xlambda = w / ( 2 * rmn );
  double tau = std::pow( xk, 2 ) - std::pow( xlambda, 2 );

   /* LEPTONIC KINEMATICAL FACTOR

   Evaluating the expressions for the leptonic tensor*/

   double xtg=-Q2/vo;
   double xdelta=xmil/sqrt(-Q2);
   double xnu=xlambda/xk;
   double xrho=1 - pow(xnu,2);
   double xrho_p=q/(ener1 + ener2);

  // cout << "tg " << xtg << ", delta " << xdelta <<", nu " << xnu <<", rho "<<xrho <<", rhop "<<xrho_p <<endl;

   double Vcc= 1 - xtg*xdelta*xdelta;
   double Vcl=-(xnu + xtg*xdelta*xdelta/xrho_p);
   double Vll=xnu*xnu + xtg*(1 + 2*xnu/xrho_p + xrho*xdelta*xdelta)*xdelta*xdelta;
   double Vt=xrho/2 + xtg -xtg*(xnu + xrho*xrho_p*xdelta*xdelta/2)*xdelta*xdelta/xrho_p;
   double Vt_p=(1 - xnu*xrho_p*xdelta*xdelta)*xtg/xrho_p;
  // double Vl=Vcc + 2*Vcl*xlambda/xk + Vll*pow(xlambda/xk,2);

  // std::cout << "Vl " << Vl << " Vcc " << Vcc << " Vcl " << Vcl << " Vll "<< Vll << " Vt "<< Vt << " Vt_p "<< Vt_p <<std::endl;


  /* INELASTIC RESPONSES
     We look in the files the value of Q² that is closest to the one which we have according
     to exchange momentum and exchange energy. */
    
   // std::cout<<" W "<< W <<" "<< rmn + w<<std::endl;

     if(W > rmn + w) return 0;

    int j=0;
     for (int m=0; m<1001000; ++m)
     {
        if( (-Q2)>=Q2vec[m] &&  (-Q2)<=Q2vec[m+1000]) break;
        
         //  std::cout << (-Q2) << " " << Q2vec[m] <<" "<<Q2vec[m + 1000] << " " << j << " " << Q2vec[j]/(4*rmn*rmn)<< std::endl;
           j=j+1;
     }

       int n_step=j;

   // We search the value of the invariant mass that it is the closest to the upper limit that the kinematics gave us.

    int f=0;
    for ( int i = 0; i < 1000; ++i )
    {
       if ( W >= Wvec[i + n_step] && W <= Wvec[i + n_step + 1] ) break;
       f=f+1;
    }
    int i_step = f;

   //    cout << n_step <<" " << i_step << endl;

        double xmux=Wvec[i_step + n_step]/rmn;
        double w1_p=w1pvec[i_step + n_step];
        double w1_n=w1nvec[i_step + n_step];
        double w2_p=w2pvec[i_step + n_step];
        double w2_n=w2nvec[i_step + n_step];
        double w3_p=w3pvec[i_step + n_step];
        double w3_n=w3nvec[i_step + n_step];

      //  std::cout <<i_step<<" "<<n_step<<std::endl;
      //  std::cout << "xmux (table) "<<xmux << " tau (table) " << Q2vec[n_step]/(4*rmn*rmn) << std::endl;
      //  std::cout<<"xmux "<<W/rmn <<" tau " <<tau << " w1p " << w1_p << " w1n " << w1_n << std::endl;
       // std::cout << "w2p " << w2_p << " w2n " << w2_n << " w3p " << w3_p << " w3n "<< w3_n<< std::endl;

       // We define the inelastic psi variable, the inelasticity paraemter and the Bjorken variable

       double rhox=1 + (xmux*xmux -1)/(4*tau);
      // double xxx=1/rhox; //The x-Bjorken variable.
       // double rnu= ( pow(xmux*rmn,2) -rmn*rmn + tau*rmn*rmn*4)/(2*rmn) //The y-Bjorken variable
       double psix=(1/sqrt(xiF))*(xlambda - tau*rhox)/sqrt(tau*(1 + xlambda*rhox) + xk*sqrt(tau*(1+ tau*rhox*rhox)) );

      // cout << "rhox " << rhox << " x "<< xxx << " psix " << psix << endl;

       //w_1
       double gmmp=sqrt(w1_p/tau);
       double gmmn=sqrt(w1_n/tau);
       double gmisos=gmmp + gmmn;
       double gmisov=gmmp-gmmn;
       double ww1isov=tau*gmisov*gmisov/2;
       double ww1isos=tau*gmisos*gmisos/2;

     //  cout << gmmp << " , " <<gmmn << " , " << gmisos << " , " << gmisov << " , "<<ww1isov << ", "<< ww1isos <<endl;

       //w_2
       double geep2geen2=(1 + tau)*(w2_p + w2_n) - (w1_p + w1_n);
       double rootgeepgeen=sqrt( ( (1 +tau)*w2_p -w1_p)*( (1 + tau)*w2_n -w1_n) );
       double ww2isov=(geep2geen2   - 2*rootgeepgeen + tau*gmisov*gmisov) /(2*(1+tau))  ;
       double ww2isos=(geep2geen2 + 2*rootgeepgeen + tau*gmisos*gmisos)/(2*(1+tau));

     //  cout << geep2geen2 << " , " << rootgeepgeen<< " , " << ww2isov << " , " << ww2isos << endl;

       //w_3
       double w3geep2geen2=(1 + tau)*(w3_p + w3_n) -(w1_p + w1_n);
       double w3rootgeepgeen=sqrt( ( (1+ tau)*w3_p - w1_p )*( (1 + tau)*w3_n - w1_n ) );
       double ww3isov=(w3geep2geen2 - 2*w3rootgeepgeen +tau*gmisov*gmisov)/(2*(1 + tau));
       double ww3isos=(w3geep2geen2 + 2*w3rootgeepgeen +tau*gmisos*gmisos)/(2*(1 + tau));

      // cout << w3geep2geen2 << " , " << w3rootgeepgeen << " , " << ww3isov << " , " << ww3isos <<endl;



       //Adding factor 18/5 in the case of neutrinos

       double ww1s=ww1isos*(18.0/5.0);
       double ww1v=ww1isov*(18.0/5.0);
       double ww2s=ww2isos*(18.0/5.0);
       double ww2v=ww2isov*(18.0/5.0);
       double ww3s=ww3isos*(18.0/5.0);
       double ww3v=ww3isov*(18.0/5.0);


      // std::cout << "w1s " << ww1s << " w1v " << ww1v << " w2s " << ww2s << " w2v " << ww2v <<std::endl;
      // std::cout << "w3s "<<ww3s << " w3v " << ww3v <<std::endl;

        //SINGLE-NUCLEON HADRONIC TENSOR

         double xD=xiF*(1 - psix*psix)*(1 + xiF*psix*psix - xlambda*psix*sqrt(xiF*(2 + xiF*psix*psix))/xk +
         tau*xiF*(1 - psix*psix)/(3*xk*xk));

        // cout <<"xD " <<xD << endl;

         //Responses for isovector

         double GCinelv=xk*xk*(ww2v*(1 + tau*rhox*rhox) - ww1v + ww2v*xD)/tau;
         double GLLinelv=GCinelv*pow(xlambda/xk,2);
         double GCLinelv=GCinelv*xlambda/xk;
         double GTinelv=2*ww1v + xD*ww2v;
         double GTpinelv=ww3v*tau*(xlambda*rhox + 1 + xiF*(1 + psix*psix)/2)/xk;

        // cout << "GCv "<<GCinelv << " GLLv " << GLLinelv << " GCLv "<<GCLinelv << " GTv "<<GTinelv << " GTpv " << GTpinelv << endl;


         //Responses for isoscalar
         double GCinels=xk*xk*( ww2s*(1 + tau*rhox*rhox) -ww1s + ww2s*xD )/tau;
         double GLLinels=GCinels*pow(xlambda/xk,2);
         double GCLinels=GCinels*xlambda/xk;
         double GTinels=2*ww1s + ww2s*xD;
         double GTpinels=ww3s*tau*(xlambda*rhox + 1 + xiF*(1 + psix*psix)/2)/xk;

        //cout << "GCs "<<GCinels << " GLLs " << GLLinels << " GCLs "<<GCLinels << " GTs "<<GTinels << " GTps " << GTpinels << endl;

        /* Defining SuSAv2-inelastic scaling function. G. D. Megias, PhD Thesis (2017)*/

        // A) SuSAv2-RMF

         //double shifte=-5; //shift inelastic parametrization for 12C
         // q en MeV to calculate Eshift_RMF and Eshift_RPWIA
         double xq=q*1000;

         // fL isovector (e, e') RMF q = 650 MeV
        /*  double a1=0.89225; // +/- 0.01183 (1.326 %)
         double a2=0.657214; // +/- 0.01588 (2.416 %)
         double a3=0.170801; // +/- 0.02409 (14.11 %)
         double a4=-0.750098; //+/- 0.08646 (11.53 %) */

         // Eshift_RMF
         double Eshift_RMF=-15.506 + 0.0548*xq + shifte;
         if (Eshift_RMF<0 && shifte>0) Eshift_RMF=0;
         if(shifte<0 && Eshift_RMF<shifte) Eshift_RMF=shifte;
          Eshift_RMF=Eshift_RMF/1000;

          //Definition of scaling variable and Psi and Psi'
          double xlambdap_RMF=(w - Eshift_RMF)/(2*rmn);
          double taup_RMF=xk*xk - xlambdap_RMF*xlambdap_RMF;
          double rhoxp_RMF=1 + (xmux*xmux -1)/(4*taup_RMF);

          double psip_RMF=(xlambdap_RMF - taup_RMF*rhoxp_RMF)/
          (sqrt(xiF)*sqrt(taup_RMF*(1 + xlambdap_RMF*rhoxp_RMF) + xk*sqrt(taup_RMF*(1 + taup_RMF*rhoxp_RMF*rhoxp_RMF))));
          double xlambdap_rmf_neg=-xlambdap_RMF;
          double psip_neg_RMF=(xlambdap_rmf_neg - taup_RMF*rhoxp_RMF)/
          (sqrt(xiF)*sqrt(taup_RMF*(1 + xlambdap_rmf_neg*rhoxp_RMF) + xk*sqrt(taup_RMF*(1 + taup_RMF*rhoxp_RMF*rhoxp_RMF))));

         // cout << "Psip " << psip_RMF << " Psip Neg "<< psip_neg_RMF << endl;

          //Defining SuSA and mirror SusA
          double fsusal=2*(a1/a2)*exp(-(psip_RMF -a3)/a2)*exp(-exp(-(psip_RMF - a3)/a2))/( 1 + exp(-a4*(psip_RMF -a3)/a2));
          double fsusal_neg=2*(a1/a2)*exp(-(psip_neg_RMF -a3)/a2)*exp(-exp(-(psip_neg_RMF - a3)/a2))/( 1 + exp(-a4*(psip_neg_RMF -a3)/a2));
          //Defining scaling function RMF
          double fLT1rmf=fsusal - fsusal_neg;
          if(fLT1rmf<0) fLT1rmf=0.0;

         // cout << "fSuSAl " <<fsusal << " fSuSAl neg " <<fsusal_neg << " fLT1rmf "<<fLT1rmf << endl;

          // B) SuSAv2-RPWIA

          /*double b1=-0.892196; // +/- 0.05081 (5.694 %)
          double b2=520.898; // +/- 1.733E4 (3327 %)
          double b3=-2906.94; // +/- 9.651E4 (3320 %)
          double b4=6475.57; //+/- 5.112E4 (789.4 %)
          double b5=1.74049; //+/- 1.767 (101.5 %)
          double b6=0.64559; // +/- 0.3952 (61.21 %) */


          // Eshift_RPWIA
          double Eshift_RPWIA=25.164 + 0.0112*xq +shifte;
          if(xq>1500) Eshift_RPWIA=25.164 + 0.0112*1500 + shifte;
          Eshift_RPWIA=Eshift_RPWIA/1000;

          //Definition of scaling variable and Psi and Psi'
          double xlambda_RPWIA=(w - Eshift_RPWIA)/(2*rmn);
          double taup_RPWIA=xk*xk - xlambda_RPWIA*xlambda_RPWIA;
          double rhoxp_RPWIA=1 + (xmux*xmux - 1)/(4*taup_RPWIA);

          //cout << "xlambda RPWIA " << xlambda_RPWIA << " taup " << taup_RPWIA << " rhox " << rhoxp_RPWIA <<endl;

          double psip_RPWIA=(xlambda_RPWIA - taup_RPWIA*rhoxp_RPWIA)/
          (sqrt(xiF)*sqrt(taup_RPWIA*(1 + xlambda_RPWIA*rhoxp_RPWIA) + xk*sqrt(taup_RPWIA*(1 + taup_RPWIA*rhoxp_RPWIA*rhoxp_RPWIA))));
          double xlambda_RPWIA_neg=-xlambda_RPWIA;
          double psip_RPWIA_neg=(xlambda_RPWIA_neg - taup_RPWIA*rhoxp_RPWIA)/
          (sqrt(xiF)*sqrt(taup_RPWIA*(1 + xlambda_RPWIA_neg*rhoxp_RPWIA) + xk*sqrt(taup_RPWIA*(1 + taup_RPWIA*rhoxp_RPWIA*rhoxp_RPWIA))));

          //cout << "Psip " << psip_RPWIA << " Psip neg " << psip_RPWIA_neg << endl;

          //Defining SuSA and mirror SusA
          double fLT1=(0.75/0.8)*(2*b4*exp(-pow(psip_RPWIA - b5,2)/b6))/(1 + exp(-b3*(psip_RPWIA - b1)/b2));
          double fLT1_neg=(0.75/0.8)*(2*b4*exp(-pow(psip_RPWIA_neg - b5,2)/b6))/(1 + exp(-b3*(psip_RPWIA_neg - b1)/b2));
          double fLT1rpwia=fLT1 - fLT1_neg;
          if(fLT1rpwia<0) fLT1rpwia=0;

         // cout << "fLT1 " << fLT1 << " fLT1 neg "<<fLT1_neg << " fLT1rpwia " << fLT1rpwia << endl;

         // DEFINITION of RMF + RPWIA SCALING FUNCTIONS
        /* double qi0=533.989; // +/- 122.5 (110.7 %)
         double qi1=0.651644; // +/- 0.8131 (31.87 %)
         double qi00=494.439; // +/- 122.5 (110.7 %)
         double qi11=0.706034; // +/- 0.8231 (31.87 %) */

         double q0imev=0.0;
         if(xq<1231.3858)  q0imev=qi0 + qi1*xq;
         if(xq>1231.3858)  q0imev=qi00 + qi11*xq;

         double q0i=q0imev/1000;
         //double w0=0.2;
         double arg=(pi/2)*( 1 - 1/(1 + exp((q - q0i)/w0)));
         double fscaling=pow(cos(arg),2)*fLT1rmf + pow(sin(arg),2)*fLT1rpwia;

         //cout << "q0 " << q0i << " arg " << arg << " fscaling " << fscaling <<endl;

         double cte=xiF/(xk*pow(etaF,3));  //Part of the expression of the response

         //Responses

         double Ginel_C=cte*( (Anumber - Znumber)*GCinels*fscaling + Znumber*GCinelv*fscaling );
         double Ginel_T=cte*( (Anumber - Znumber)*GTinels*fscaling + Znumber*GTinelv*fscaling );
         double Ginel_LL=cte*( (Anumber - Znumber)*GLLinels*fscaling + Znumber*GLLinelv*fscaling );
         double Ginel_CL=cte*( (Anumber - Znumber)*GCLinels*fscaling + Znumber*GCLinelv*fscaling );
         double Ginel_Tp=cte*( (Anumber - Znumber)*GTpinels*fscaling + Znumber*GTpinelv*fscaling );

         //cout << "cte " << cte << " GC " << Ginel_C << " GT " << Ginel_T << " GLL "<<Ginel_LL << endl;
         //cout << Ginel_CL << " GTp " << Ginel_Tp << endl


      /* INELASTIC CROSS SECTION */

      // Tensor contraction

      double Fx2=Vcc*Ginel_C + Vt*Ginel_T + Vll*Ginel_LL + 2*Vcl*Ginel_CL + inu*2*Vt_p*Ginel_Tp;

      //Definition of Sigma0 and units
      //double unit=hc*hc*1E-26; //to transform to cm^2/GeV/str dsdomegadKp
      double sigma0=vo*(GF*GF*cos_cabibbo*cos_cabibbo*xkp*xkp)/(8*pi*pi*ener1*ener2);

      // cout << "Fx2 " << Fx2 << " unit " << unit << " sigma0 " << sigma0 << endl;


      //CROSS SECTION

    //double   CS_inel=sigma0*Fx2*unit;

    double CS_inel = sigma0 * Fx2;

   //cout <<"d^2 sigma/d Omega d k{^prime} (GeV^{-2}/GeV/str) "
   //  <<CS_inel << endl;


   // If the hit nucleon PDG code is set, scale by N/A (neutron) or Z/A
   // so that the inclusive prediction is recovered for the sum of the
   // hit-nucleon-specific xsecs
   if ( target.HitNucIsSet() ) {
     int hit_nuc_pdg = target.HitNucPdg();
     if ( genie::pdg::IsNeutron(hit_nuc_pdg) ) {
       CS_inel *= static_cast< double >( Anumber - Znumber ) / Anumber;
     }
     else if ( genie::pdg::IsProton(hit_nuc_pdg) ) {
       CS_inel *= static_cast< double >( Znumber ) / Anumber;
     }
     else {
       LOG( "SuSAv2Inel", pERROR ) << "Unrecognized hit nucleon PDG code";
       CS_inel = 0.;
     }
   }

   if ( std::isnan(CS_inel) || CS_inel < 0. ) CS_inel = 0.0;

   return CS_inel;
}
// END CODE FROM JESUS
//////////////////////////////////////////////////////////////
//____________________________________________________________________________
double genie::SuSAv2InelPXSec::Integral(
  const genie::Interaction* interaction ) const
{
  const genie::InitialState& init_state = interaction->InitState();
   const genie::Target& target = init_state.Tgt();
  double xmil = interaction->FSPrimLepton()->Mass();
  double ev = init_state.ProbeE( genie::kRfLab );
  double pi = std::acos( -1.0 );
  double rmn = genie::constants::kNucleonMass;
  double W_min = rmn + genie::constants::kPionMass;
  double GF = genie::constants::kGF;
    // Fermi momentum
  double etaF = pF / rmn;
  double epsF = std::sqrt( 1.0 + etaF*etaF );
  double xiF = epsF - 1.0;
  
  int Anumber = target.A();
  int Znumber = target.Z();
  int  probepdgc = init_state.ProbePdg();
 int inu = 0;
 
  if ( genie::pdg::IsNeutrino(probepdgc) ) inu = 1;
  else if ( genie::pdg::IsAntiNeutrino(probepdgc) ) inu = -1;
  else {
    // TODO: revisit this!
    LOG( "SuSAv2InelPXSec", pFATAL ) << "Bad inu value for probe!";
    std::exit( 1 );
  }

  double Suma_Tl = 0.0;
  for ( int i = 0; i < 500; i++ ) {
    double Tl = i * ( ev - xmil ) / 500;
    double Suma_cos = 0.0;
    for ( int j = 0; j < 500; j++ ) {
      double costh = -1.00 + j*2.00/500;
      double Suma_W = 0.0;
      
      
   //KINEMATICS
   double ener2=Tl + xmil; //Energy of the lepton
   double xkp=sqrt(ener2*ener2 - xmil*xmil); //Momentum of the lepton
   double w=ev-ener2;
   double q=sqrt(ev*ev + xkp*xkp - 2*xkp*ev*costh);
   double Q2=w*w - q*q; //Four-momentum
  double vo=4*ev*ener2 + Q2;
  // double vo=pow((ener1 + ener2),2) - pow(q,2);
  // Scaling and Adimensional variables


   double xk=q/(2*rmn); 
   double xlambda=w/(2*rmn);
   double tau= pow(xk,2) - pow(xlambda,2);
   // double tau=-Q2/(4*rmn*rmn);
   //double psi=(1/sqrt(xiF))*(xlambda*tau)/sqrt( (1+ xlambda)*tau + xk*sqrt(tau*(tau + 1)))

   //cout << xk << ", " << xlambda << " , " << tau << endl;

   /* LEPTONIC KINEMATICAL FACTOR

   Evaluating the expressions for the leptonic tensor*/

   double xtg=-Q2/vo;
   double xdelta=xmil/sqrt(-Q2);
   double xnu=xlambda/xk;
   double xrho=1 - pow(xnu,2);
   double xrho_p=q/(ev + ener2);

  // cout << "tg " << xtg << ", delta " << xdelta <<", nu " << xnu <<", rho "<<xrho <<", rhop "<<xrho_p <<endl;

   double Vcc= 1 - xtg*xdelta*xdelta;
   double Vcl=-(xnu + xtg*xdelta*xdelta/xrho_p);
   double Vll=xnu*xnu + xtg*(1 + 2*xnu/xrho_p + xrho*xdelta*xdelta)*xdelta*xdelta;
   double Vt=xrho/2 + xtg -xtg*(xnu + xrho*xrho_p*xdelta*xdelta/2)*xdelta*xdelta/xrho_p;
   double Vt_p=(1 - xnu*xrho_p*xdelta*xdelta)*xtg/xrho_p;
  // double Vl=Vcc + 2*Vcl*xlambda/xk + Vll*pow(xlambda/xk,2);
   
   double W_max =genie::constants::kNeutronMass + ev - Tl - xmil;
   
     /* INELASTIC RESPONSES
     We look in the files the value of Q² that is closest to the one which we have according
     to exchange momentum and exchange energy. */
    
    int t=0;
     for (int m=0; m<1001000; ++m)
     {           
        if( abs(Q2)>=Q2vec[m]&&abs(Q2)<=Q2vec[m+1000]) break;
           t=t+1;
 
 
     }

       int n_step=t;
   

          if ( W_max > W_min ) {

         double W_step = ( W_max - W_min ) / 500;
      for( int z = 0; z < 500; z++ ) {
        double W = W_min + W_step*z;
   
      int c=0;
     for (int r=0; r<1000; ++r)
     {   
        if( W>=Wvec[r + n_step]&&W<=Wvec[r + n_step +1]) break;
        c=c+1;
     }
       int f_step=c;
   
        double xmux=Wvec[f_step + n_step]/rmn;
        double w1_p=w1pvec[f_step + n_step];
        double w1_n=w1nvec[f_step + n_step];
        double w2_p=w2pvec[f_step + n_step];
        double w2_n=w2nvec[f_step + n_step];
        double w3_p=w3pvec[f_step + n_step];
        double w3_n=w3nvec[f_step + n_step];    

 // We define the inelastic psi variable, the inelasticity paraemter and the Bjorken variable

       double rhox=1 + (xmux*xmux -1)/(4*tau);
      // double xxx=1/rhox; //The x-Bjorken variable.
       // double rnu= ( pow(xmux*rmn,2) -rmn*rmn + tau*rmn*rmn*4)/(2*rmn) //The y-Bjorken variable
       double psix=(1/sqrt(xiF))*(xlambda - tau*rhox)/sqrt(tau*(1 + xlambda*rhox) + xk*sqrt(tau*(1+ tau*rhox*rhox)) );
       
      // cout << "rhox " << rhox << " x "<< xxx << " psix " << psix << endl;

       //w_1
       double gmmp=sqrt(w1_p/tau);
       double gmmn=sqrt(w1_n/tau);
       double gmisos=gmmp + gmmn;
       double gmisov=gmmp-gmmn;
       double ww1isov=tau*gmisov*gmisov/2;
       double ww1isos=tau*gmisos*gmisos/2;

     //  cout << gmmp << " , " <<gmmn << " , " << gmisos << " , " << gmisov << " , "<<ww1isov << ", "<< ww1isos <<endl;

       //w_2
       double geep2geen2=(1 + tau)*(w2_p + w2_n) - (w1_p + w1_n);
       double rootgeepgeen=sqrt( ( (1 +tau)*w2_p -w1_p)*( (1 + tau)*w2_n -w1_n) );
       double ww2isov=(geep2geen2   - 2*rootgeepgeen + tau*gmisov*gmisov) /(2*(1+tau))  ;
       double ww2isos=(geep2geen2 + 2*rootgeepgeen + tau*gmisos*gmisos)/(2*(1+tau));

     //  cout << geep2geen2 << " , " << rootgeepgeen<< " , " << ww2isov << " , " << ww2isos << endl;

       //w_3
       double w3geep2geen2=(1 + tau)*(w3_p + w3_n) -(w1_p + w1_n);
       double w3rootgeepgeen=sqrt( ( (1+ tau)*w3_p - w1_p )*( (1 + tau)*w3_n - w1_n ) );
       double ww3isov=(w3geep2geen2 - 2*w3rootgeepgeen +tau*gmisov*gmisov)/(2*(1 + tau));
       double ww3isos=(w3geep2geen2 + 2*w3rootgeepgeen +tau*gmisos*gmisos)/(2*(1 + tau));

      // cout << w3geep2geen2 << " , " << w3rootgeepgeen << " , " << ww3isov << " , " << ww3isos <<endl;



       //Adding factor 18/5 in the case of neutrinos

       double ww1s=ww1isos*(18.0/5.0);
       double ww1v=ww1isov*(18.0/5.0);
       double ww2s=ww2isos*(18.0/5.0);
       double ww2v=ww2isov*(18.0/5.0);
       double ww3s=ww3isos*(18.0/5.0);
       double ww3v=ww3isov*(18.0/5.0);


     //  cout << "w1s " << ww1s << " w1v " << ww1v << " w2s " << ww2s << " w2v " << ww2v <<endl;
     //  cout << "w3s "<<ww3s << " w3v " << ww3v <<endl;

        //SINGLE-NUCLEON HADRONIC TENSOR

         double xD=xiF*(1 - psix*psix)*(1 + xiF*psix*psix - xlambda*psix*sqrt(xiF*(2 + xiF*psix*psix))/xk + 
         tau*xiF*(1 - psix*psix)/(3*xk*xk));

        // cout <<"xD " <<xD << endl;

         //Responses for isovector 

         double GCinelv=xk*xk*(ww2v*(1 + tau*rhox*rhox) - ww1v + ww2v*xD)/tau;
         double GLLinelv=GCinelv*pow(xlambda/xk,2);
         double GCLinelv=GCinelv*xlambda/xk;
         double GTinelv=2*ww1v + xD*ww2v;
         double GTpinelv=ww3v*tau*(xlambda*rhox + 1 + xiF*(1 + psix*psix)/2)/xk;

        // cout << "GCv "<<GCinelv << " GLLv " << GLLinelv << " GCLv "<<GCLinelv << " GTv "<<GTinelv << " GTpv " << GTpinelv << endl;


         //Responses for isoscalar
         double GCinels=xk*xk*( ww2s*(1 + tau*rhox*rhox) -ww1s + ww2s*xD )/tau;
         double GLLinels=GCinels*pow(xlambda/xk,2);
         double GCLinels=GCinels*xlambda/xk;
         double GTinels=2*ww1s + ww2s*xD;
         double GTpinels=ww3s*tau*(xlambda*rhox + 1 + xiF*(1 + psix*psix)/2)/xk;

        //cout << "GCs "<<GCinels << " GLLs " << GLLinels << " GCLs "<<GCLinels << " GTs "<<GTinels << " GTps " << GTpinels << endl;      

        /* Defining SuSAv2-inelastic scaling function. G. D. Megias, PhD Thesis (2017)*/

        // A) SuSAv2-RMF

       //  double shifte=-5; //shift inelastic parametrization for 12C
         // q en MeV to calculate Eshift_RMF and Eshift_RPWIA
         double xq=q*1000;

         // fL isovector (e, e') RMF q = 650 MeV
       /*  double a1=0.89225; // +/- 0.01183 (1.326 %)
         double a2=0.657214; // +/- 0.01588 (2.416 %)
         double a3=0.170801; // +/- 0.02409 (14.11 %)
         double a4=-0.750098; //+/- 0.08646 (11.53 %) */

         // Eshift_RMF
         double Eshift_RMF=-15.506 + 0.0548*xq + shifte;
         if ((Eshift_RMF<0)&(shifte>0)) Eshift_RMF=0;
         if((shifte<0)&(Eshift_RMF<shifte)) Eshift_RMF=shifte;
          Eshift_RMF=Eshift_RMF/1000;

          //Definition of scaling variable and Psi and Psi'
          double xlambdap_RMF=(w - Eshift_RMF)/(2*rmn);
          double taup_RMF=xk*xk - xlambdap_RMF*xlambdap_RMF;
          double rhoxp_RMF=1 + (xmux*xmux -1)/(4*taup_RMF); 
          
          double psip_RMF=(xlambdap_RMF - taup_RMF*rhoxp_RMF)/
          (sqrt(xiF)*sqrt(taup_RMF*(1 + xlambdap_RMF*rhoxp_RMF) + xk*sqrt(taup_RMF*(1 + taup_RMF*rhoxp_RMF*rhoxp_RMF))));
          double xlambdap_rmf_neg=-xlambdap_RMF;
          double psip_neg_RMF=(xlambdap_rmf_neg - taup_RMF*rhoxp_RMF)/
          (sqrt(xiF)*sqrt(taup_RMF*(1 + xlambdap_rmf_neg*rhoxp_RMF) + xk*sqrt(taup_RMF*(1 + taup_RMF*rhoxp_RMF*rhoxp_RMF))));

         // cout << "Psip " << psip_RMF << " Psip Neg "<< psip_neg_RMF << endl;

          //Defining SuSA and mirror SusA
          double fsusal=2*(a1/a2)*exp(-(psip_RMF -a3)/a2)*exp(-exp(-(psip_RMF - a3)/a2))/( 1 + exp(-a4*(psip_RMF -a3)/a2));
          double fsusal_neg=2*(a1/a2)*exp(-(psip_neg_RMF -a3)/a2)*exp(-exp(-(psip_neg_RMF - a3)/a2))/( 1 + exp(-a4*(psip_neg_RMF -a3)/a2));
          //Defining scaling function RMF
          double fLT1rmf=fsusal - fsusal_neg;
          if(fLT1rmf<0) fLT1rmf=0.0;
           
         // cout << "fSuSAl " <<fsusal << " fSuSAl neg " <<fsusal_neg << " fLT1rmf "<<fLT1rmf << endl;

          // B) SuSAv2-RPWIA

         /* double b1=-0.892196; // +/- 0.05081 (5.694 %)
          double b2=520.898; // +/- 1.733E4 (3327 %)
          double b3=-2906.94; // +/- 9.651E4 (3320 %)
          double b4=6475.57; //+/- 5.112E4 (789.4 %)
          double b5=1.74049; //+/- 1.767 (101.5 %)
          double b6=0.64559; // +/- 0.3952 (61.21 %) */

          // Eshift_RPWIA
          double Eshift_RPWIA=25.164 + 0.0112*xq +shifte;
          if(xq>1500) Eshift_RPWIA=25.164 + 0.0112*1500 + shifte;
          Eshift_RPWIA=Eshift_RPWIA/1000;

          //Definition of scaling variable and Psi and Psi'
          double xlambda_RPWIA=(w - Eshift_RPWIA)/(2*rmn);
          double taup_RPWIA=xk*xk - xlambda_RPWIA*xlambda_RPWIA;
          double rhoxp_RPWIA=1 + (xmux*xmux - 1)/(4*taup_RPWIA);

          //cout << "xlambda RPWIA " << xlambda_RPWIA << " taup " << taup_RPWIA << " rhox " << rhoxp_RPWIA <<endl;
        
          double psip_RPWIA=(xlambda_RPWIA - taup_RPWIA*rhoxp_RPWIA)/
          (sqrt(xiF)*sqrt(taup_RPWIA*(1 + xlambda_RPWIA*rhoxp_RPWIA) + xk*sqrt(taup_RPWIA*(1 + taup_RPWIA*rhoxp_RPWIA*rhoxp_RPWIA))));
          double xlambda_RPWIA_neg=-xlambda_RPWIA;
          double psip_RPWIA_neg=(xlambda_RPWIA_neg - taup_RPWIA*rhoxp_RPWIA)/
          (sqrt(xiF)*sqrt(taup_RPWIA*(1 + xlambda_RPWIA_neg*rhoxp_RPWIA) + xk*sqrt(taup_RPWIA*(1 + taup_RPWIA*rhoxp_RPWIA*rhoxp_RPWIA))));

          //cout << "Psip " << psip_RPWIA << " Psip neg " << psip_RPWIA_neg << endl;

          //Defining SuSA and mirror SusA
          double fLT1=(0.75/0.8)*(2*b4*exp(-pow(psip_RPWIA - b5,2)/b6))/(1 + exp(-b3*(psip_RPWIA - b1)/b2)); 
          double fLT1_neg=(0.75/0.8)*(2*b4*exp(-pow(psip_RPWIA_neg - b5,2)/b6))/(1 + exp(-b3*(psip_RPWIA_neg - b1)/b2)); 
          double fLT1rpwia=fLT1 - fLT1_neg;
          if(fLT1rpwia<0) fLT1rpwia=0;

         // cout << "fLT1 " << fLT1 << " fLT1 neg "<<fLT1_neg << " fLT1rpwia " << fLT1rpwia << endl;

         // DEFINITION of RMF + RPWIA SCALING FUNCTIONS
        /* double qi0=533.989; // +/- 122.5 (110.7 %)
         double qi1=0.651644; // +/- 0.8131 (31.87 %)
         double qi00=494.439; // +/- 122.5 (110.7 %)
         double qi11=0.706034; // +/- 0.8231 (31.87 %) */
        
         double q0imev=0.0; 
         if(xq<1231.3858)  q0imev=qi0 + qi1*xq;
         if(xq>1231.3858)  q0imev=qi00 + qi11*xq;

         double q0i=q0imev/1000;
       //  double w0=0.2;
         double arg=(pi/2)*( 1 - 1/(1 + exp((q - q0i)/w0)));
         double fscaling=pow(cos(arg),2)*fLT1rmf + pow(sin(arg),2)*fLT1rpwia;

         //cout << "q0 " << q0i << " arg " << arg << " fscaling " << fscaling <<endl;  

         double cte=xiF/(xk*pow(etaF,3));  //Part of the expression of the response

         //Responses

         double Ginel_C=cte*( (Anumber - Znumber)*GCinels*fscaling + Znumber*GCinelv*fscaling );
         double Ginel_T=cte*( (Anumber - Znumber)*GTinels*fscaling + Znumber*GTinelv*fscaling );
         double Ginel_LL=cte*( (Anumber - Znumber)*GLLinels*fscaling + Znumber*GLLinelv*fscaling );         
         double Ginel_CL=cte*( (Anumber - Znumber)*GCLinels*fscaling + Znumber*GCLinelv*fscaling );
         double Ginel_Tp=cte*( (Anumber - Znumber)*GTpinels*fscaling + Znumber*GTpinelv*fscaling );

         //cout << "cte " << cte << " GC " << Ginel_C << " GT " << Ginel_T << " GLL "<<Ginel_LL << endl;
         //cout << Ginel_CL << " GTp " << Ginel_Tp << endl         



          
      /* INELASTIC CROSS SECTION */

      // Tensor contraction

      double Fx2=Vcc*Ginel_C + Vt*Ginel_T + Vll*Ginel_LL + 2*Vcl*Ginel_CL + inu*2*Vt_p*Ginel_Tp;
      double sigma0=vo*(GF*GF*cos_cabibbo*cos_cabibbo*xkp*xkp)/(8*pi*pi*ev*ener2);

     double   CS_inel=sigma0*Fx2; 
     
 if ( target.HitNucIsSet() ) {
     int hit_nuc_pdg = target.HitNucPdg();
     if ( genie::pdg::IsNeutron(hit_nuc_pdg) ) {
       CS_inel *= static_cast< double >( Anumber - Znumber ) / Anumber;
     }
     else if ( genie::pdg::IsProton(hit_nuc_pdg) ) {
       CS_inel *= static_cast< double >( Znumber ) / Anumber;
     }
     else {
       LOG( "SuSAv2Inel", pERROR ) << "Unrecognized hit nucleon PDG code";
       CS_inel = 0.;
     }
   }

   if ( std::isnan(CS_inel) || CS_inel < 0. ) CS_inel = 0.0;
     
        Suma_W = Suma_W + 2*pi*W*CS_inel*W_step/pow( rmn,2 );
        }
        } else Suma_W=0.0;
      Suma_cos = Suma_cos + Suma_W*2.00/500;
    }
    Suma_Tl = Suma_Tl + Suma_cos*( ev - xmil ) / 500;
  }
  double xsec_tot = Suma_Tl;
  return xsec_tot;
}
//____________________________________________________________________________
bool genie::SuSAv2InelPXSec::ValidProcess(
  const genie::Interaction* interaction ) const
{
  if ( interaction->TestBit(kISkipProcessChk) ) return true;

  const genie::InitialState& init_state = interaction->InitState();
  const genie::ProcessInfo&  proc_info  = interaction->ProcInfo();
  //const genie::XclsTag&      xcls       = interaction->ExclTag();

  if ( !proc_info.IsResonant() ) return false;
  //if ( !xcls.KnownResonance() )  return false;

  int  hitnuc = init_state.Tgt().HitNucPdg();
  bool is_pn = ( pdg::IsProton(hitnuc) || pdg::IsNeutron(hitnuc) );

  if ( !is_pn ) return false;

  int  probe   = init_state.ProbePdg();
  bool is_weak = proc_info.IsWeak();
  bool is_em   = proc_info.IsEM();
  bool nu_weak = ( pdg::IsNeutralLepton(probe) && is_weak );
  bool l_em    = ( pdg::IsChargedLepton(probe) && is_em );

  if ( !nu_weak && !l_em ) return false;

  return true;
}
//____________________________________________________________________________
void genie::SuSAv2InelPXSec::Configure(const genie::Registry & config)
{
  genie::Algorithm::Configure( config );
  this->LoadConfig();
}
//____________________________________________________________________________
void genie::SuSAv2InelPXSec::Configure( std::string config )
{
  genie::Algorithm::Configure( config );
  this->LoadConfig();
}
//____________________________________________________________________________
void genie::SuSAv2InelPXSec::LoadConfig()
{
  //// Cross section scaling factors
  //this->GetParam( "RES-CC-XSecScale", fXSecScaleCC ) ;
  //this->GetParam( "RES-EM-XSecScale", fXSecScaleEM ) ;

  // TODO: read in the Fermi momentum in the usual GENIE way
  this->GetParam( "FermiMomentum", pF );
  //this->GetParam( "FermiMomentumTable", fKFTable );

  //this->GetParam( "Eshift",esep );

  this->GetParam( "Eshift-inel", shifte );

  this->GetParam( "param-RMF-1", a1 );
  this->GetParam( "param-RMF-2", a2 );
  this->GetParam( "param-RMF-3", a3 );
  this->GetParam( "param-RMF-4", a4 );

  this->GetParam( "param-RPWIA-1", b1 );
  this->GetParam( "param-RPWIA-2", b2 );
  this->GetParam( "param-RPWIA-3", b3 );
  this->GetParam( "param-RPWIA-4", b4 );
  this->GetParam( "param-RPWIA-5", b5 );
  this->GetParam( "param-RPWIA-6", b6 );

  this->GetParam( "param-transit-q-1", qi0 );
  this->GetParam( "param-transit-q-2", qi1 );
  this->GetParam( "param-transit-q-3", qi00 );
  this->GetParam( "param-transit-q-4", qi11 );
  this->GetParam( "param-transit-w", w0 );

  double ang_cabibbo;
  this->GetParam( "CabibboAngle", ang_cabibbo );
  cos_cabibbo = std::cos( ang_cabibbo );

  std::string struc_func_file_name;
  this->GetParam( "InelStrucFuncFile", struc_func_file_name );

  this->LoadStructureFunctions( struc_func_file_name );

  // Load the cross section integrator
  fXSecIntegrator = dynamic_cast< const genie::XSecIntegratorI* >(
    this->SubAlg("XSec-Integrator") );
  assert( fXSecIntegrator );
}
//____________________________________________________________________________
void genie::SuSAv2InelPXSec::LoadStructureFunctions(
  const std::string& input_file_name )
{
  // Clear out old contents of the structure function table. Also reserve an
  // appropriate amount of memory for the new table elements.
  Q2vec.clear();
  Q2vec.reserve( VECTOR_SIZE );

  Wvec.clear();
  Wvec.reserve( VECTOR_SIZE );

  w1pvec.clear();
  w1pvec.reserve( VECTOR_SIZE );

  w1nvec.clear();
  w1nvec.reserve( VECTOR_SIZE );

  w2pvec.clear();
  w2pvec.reserve( VECTOR_SIZE );

  w2nvec.clear();
  w2nvec.reserve( VECTOR_SIZE );

  w3pvec.clear();
  w3pvec.reserve( VECTOR_SIZE );

  w3nvec.clear();
  w3nvec.reserve( VECTOR_SIZE );

  // Temporary storage for reading in the table
  double Q2, W, w1p, w1n, w2p, w2n, w3p, w3n;

  std::string full_file_name = std::string( gSystem->Getenv("GENIE") )
    + "/data/evgen/nucl/" + input_file_name;

  std::ifstream in_file( full_file_name );

  // TODO: add error handling for a missing file
  while ( in_file >> Q2 >> W >> w1p >> w1n >> w2p >> w2n >> w3p >> w3n ) {
    Q2vec.push_back( Q2 );
    Wvec.push_back( W );
    w1pvec.push_back( w1p );
    w1nvec.push_back( w1n );
    w2pvec.push_back( w2p );
    w2nvec.push_back( w2n );
    w3pvec.push_back( w3p );
    w3nvec.push_back( w3n );
  }

}
