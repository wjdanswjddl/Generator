//____________________________________________________________________________
/*
 Copyright (c) 2003-2021, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Syrian Truong <struong@fnal.gov>
         Steven Gardiner <gardiner \at fnal.gov>
         Fermi National Acclerator Laboratory

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Math/IFunction.h"
#include "Math/Integrator.h"

#include "Physics/QuasiElastic/XSection/QELUtils.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Physics/NuclearState/NuclearModelI.h"
#include "Physics/NuclearState/NuclearUtils.h"

#include "FortranWrapperQELPXSec.h"
#include "FortranWrapperXSecIntegrator.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;
using namespace genie::utils;

//____________________________________________________________________________
extern"C"
{
void diracmatrices_(double *xmn_in);
}

//____________________________________________________________________________
extern"C"
{
void cc1_(double *xq, double *w, double *wt, double *xk, double *xp, double *ee0, double *theta, int *ig, double *xsec, double *nuphi);
}

//____________________________________________________________________________
FortranWrapperQELPXSec::FortranWrapperQELPXSec() :
XSecAlgorithmI("genie::FortranWrapperQELPXSec")
{

}
//____________________________________________________________________________
FortranWrapperQELPXSec::FortranWrapperQELPXSec(string config) :
XSecAlgorithmI("genie::FortranWrapperQELPXSec", config)
{

}
//____________________________________________________________________________
FortranWrapperQELPXSec::~FortranWrapperQELPXSec()
{

}
//____________________________________________________________________________
double FortranWrapperQELPXSec::XSec(const Interaction* interaction,
  KinePhaseSpace_t kps) const
{
  // Apply the global Q^2 cut imposed by GENIE
  TLorentzVector* tempProbeP4 = interaction->InitState().GetProbeP4();
  TLorentzVector probeP4 = *tempProbeP4;
  delete tempProbeP4;

  const TLorentzVector& fsLeptonP4 = interaction->Kine().FSLeptonP4();
  // 4-momentum transfer
  TLorentzVector qP4 = probeP4 - fsLeptonP4;
  double Q2 = -1. * qP4.Mag2();

  interaction->KinePtr()->SetQ2( Q2 );

  if ( Q2 < genie::utils::kinematics::electromagnetic::kMinQ2Limit ) return 0.;

  // This is the main "missing piece" for interfacing with Noemi's code.
  // Given an Interaction object which will be pre-loaded with the correct
  // variables, we need to extract the inputs to Noemi's function, call that
  // function, and then return a differential cross section in natural units.
  // You can ignore the kps variable for now. Just work with the output
  // of Noemi's cc1 subroutine here.

  // Input 1 of 11: xmn_in

   double pm_GeV = genie::PDGLibrary::Instance()->Find(2212)->Mass(); // Proton mass

   double nm_GeV = genie::PDGLibrary::Instance()->Find(2112)->Mass(); // Neutron mass

   double xmn_in = 0.5 * 1000 * (pm_GeV + nm_GeV) / 197.327053; // (proton_{mass} + neutron_{mass})/(2*h_bar*c) in units of 1/fm*c^2 because (MeV/c^2)*[1/(MeV*fm)]

   std::cout << "\n1 of 11) xmn_in = " << xmn_in << " 1/fm*c^2\n"; // double xmn_in1 = 4.7581860861217038;

  // Input 2 of 11: xq

   genie::InitialState* initstate = interaction->InitStatePtr(); // initial state functions

   TLorentzVector* ilep4 = initstate->GetProbeP4(genie::kRfLab); // initial lepton 4-momenta

   genie::Kinematics* kinmat = interaction->KinePtr(); // final state functions

   TLorentzVector flep4 = kinmat->FSLeptonP4(); // final lepton 4-momenta

   double il_Px = ilep4->Px()/(197.3269804 * 0.001); // initial lepton momentum in x, from GeV to 1/fm

   double il_Py = ilep4->Py()/(197.3269804 * 0.001); // initial lepton momentum in y, from GeV to 1/fm

   double il_Pz = ilep4->Pz()/(197.3269804 * 0.001); // initial lepton momentum in z, from GeV to 1/fm

   double fl_Px = flep4.Px()/(197.3269804 * 0.001); // final lepton momentum in x, from GeV to 1/fm

   double fl_Py = flep4.Py()/(197.3269804 * 0.001); // final lepton momentum in y, from GeV to 1/fm

   double fl_Pz = flep4.Pz()/(197.3269804 * 0.001); // final lepton momentum in z, from GeV to 1/fm

   double xq = sqrt( (il_Px - fl_Px)*(il_Px - fl_Px) + (il_Py - fl_Py)*(il_Py - fl_Py) + (il_Pz - fl_Pz)*(il_Pz - fl_Pz) ); // double xq = 2.2944884204923253;

   std::cout << "2 of 11) xq = " << xq << " 1/fm\n"; // three-momentum transfer in inverse fm (i.e. magnitude of |incident electron 3-momentum - final electron 3-momentum|)

  // Input 3 of 11: w

   double il_E = ilep4->E()*1000; // initial lepton energy, converted from GeV to MeV

   double fl_E = flep4.E()*1000; // final lepton energy, converted from GeV to MeV

   double w = sqrt((il_E - fl_E) * (il_E - fl_E)); // double w = 37.500000000000000;

   std::cout << "3 of 11) w = " << w << " MeV\n"; // lepton energy loss in MeV (i.e. |incident lepton energy - final lepton energy|)

  // Input 4 of 11: wt

   // w defined above

   double pm_MeV = 1000 * pm_GeV; // (proton_{mass} + neutron_{mass})/2 in units of (MeV/c^2), but note Noemi's code uses (proton_{mass} + neutron_{mass})/2 in units of (MeV/c^2), however, the difference is likely negligible, so we'll use proton mass here

   genie::Target inuc4 = initstate->Tgt(); // initial target functions for initial nucleon 4-momentum

   TLorentzVector* inuc4hit = inuc4.HitNucP4Ptr(); // initial hit functions for initial nucleon 4-momentum

   double in_E_MeV = inuc4hit->E()*1000; // initial nucleon energy, from GeV to MeV

   double in_Px_MeV = inuc4hit->Px()*1000; // initial nucleon momentum in x, from GeV to MeV

   double in_Py_MeV = inuc4hit->Py()*1000; // initial nucleon momentum in y, from GeV to MeV

   double in_Pz_MeV = inuc4hit->Pz()*1000; // initial nucleon momentum in z, from GeV to MeV

   double xk_MeV = sqrt(in_Px_MeV * in_Px_MeV + in_Py_MeV * in_Py_MeV + in_Pz_MeV * in_Pz_MeV); // |initial nucleon (proton) 3-momenta|

   double wt = w + in_E_MeV - sqrt(pm_MeV * pm_MeV + xk_MeV * xk_MeV); // |initial lepton (electron) energy - final lepton (electron) energy)| + initial nucleon (proton) energy - \sqrt{ [nucleon (proton) mass]^2 + |initial nucleon (proton) 3-momenta|^2 } = w + M - E - E_{\mathbf{p}} = E_{\mathbf{p'}} - E_{\mathbf{p}}

   std::cout << "4 of 11) wt = " << wt << " MeV\n"; // Need to make sure units are correct, but everything should be in MeV I believe as long as definitions for hbar and c are consistent

  // Input 5 of 11: xk

   // initial target functions defined above

   // initial hit functions defined above

   double in_Px = inuc4hit->Px()/(197.3269804 * 0.001); // initial nucleon momentum in x, from GeV to 1/fm

   double in_Py = inuc4hit->Py()/(197.3269804 * 0.001); // initial nucleon momentum in y, from GeV to 1/fm

   double in_Pz = inuc4hit->Pz()/(197.3269804 * 0.001); // initial nucleon momentum in z, from GeV to 1/fm

   double xk = sqrt(in_Px * in_Px + in_Py * in_Py + in_Pz * in_Pz); // double xk = 1.0642230591666515;

   std::cout << "5 of 11) xk = " << xk << " 1/fm\n"; // initial nucleon three-momenta in inverse fm (i.e. magnitude of |initial nucleon 3-momenta|)

  // Input 6 of 11: xp

   TLorentzVector fhad4 = kinmat->HadSystP4(); // final hadron 4-momenta

   // Check outgoing nucleon momentum against kF. If it is Pauli-blocked, then
   // just return zero without bothering to do the rest of the calculation.
   double pNf = fhad4.P();
   double kF = fNuclModel->LocalFermiMomentum( inuc4, inuc4.HitNucPdg(),
     inuc4hit->Vect().Mag() );
   if ( pNf < kF ) return 0.;

   double fh_E = fhad4.E()*1000; // final hadron E, converted from GeV to MeV

   double fh_Px = fhad4.Px()/(197.3269804 * 0.001); // final hadron momentum in x, from GeV to 1/fm

   double fh_Py = fhad4.Py()/(197.3269804 * 0.001); // final hadron momentum in y, from GeV to 1/fm

   double fh_Pz = fhad4.Pz()/(197.3269804 * 0.001); // final hadron momentum in z, from GeV to 1/fm

   double xp = sqrt(fh_Px * fh_Px + fh_Py * fh_Py + fh_Pz * fh_Pz); // double xp = 1.3112528672942025;

   std::cout << "6 of 11) xp = " << xp << " 1/fm\n"; // final nucleon (or hadron) three-momenta in inverse fm (i.e. magnitude of |final hadron 3-momenta|)

  // Input 7 of 11: ee0

   double ee0 = il_E; // double ee0 = 730.00000000000000;

   std::cout << "7 of 11) ee0 = " << ee0 << " MeV\n"; // incident lepton energy in MeV

  // Input 8 of 11: theta

   // initial lepton momentum defined above in components

   // final lepton momentum defined above in components

   double xil = sqrt(il_Px * il_Px + il_Py * il_Py + il_Pz * il_Pz); // Magnitude of initial lepton 3-momentum

   double xfl = sqrt(fl_Px * fl_Px + fl_Py * fl_Py + fl_Pz * fl_Pz); // Magnitude of final lepton 3-momentum

   double theta = acos( (il_Px * fl_Px + il_Py * fl_Py + il_Pz * fl_Pz) / ( xil * xfl) ); // double theta = 0.64577182323790194;

   std::cout << "8 of 11) theta = " << theta << " rad (or " << theta * 180 / M_PI << " degrees)\n"; // angle between the initial lepton 3-momentum and final lepton 3-momentum

  // Input 9 of 11: ig

   int ig = 2;

   std::cout << "9 of 11) ig = " << ig << " Setting? (need to figure out what this actually is in Noemi's code)\n"; // variable for setting purposes in Noemi's code?

  // Input 10 of 11: xsec

   double xsec;

   std::cout << "10 of 11) xsec is passed in without a defined value, then given one via cc1 below:" << '\n';

  // Input 11 of 11: nuphi

  // initial nucleon or hadron momentum defined above in components

  // final nucleon or hadron momentum defined above in components

  // magnitude of initial nucleon or hadron 3-momentum defined above

  // magnitude of final nucleon or hadrom 3-momentum defined above

  double nuphi = acos( (in_Px * fh_Px + in_Py * fh_Py + in_Pz * fh_Pz) / ( xk * xp) ); // angle between the initial nucleon 3-momentum and final nucleon 3-momentum

   std::cout << "11 of 11) nuphi = " << nuphi << " rad (or " << nuphi * 180 / M_PI << " degrees)\n"; // angle between the initial nucleon 3-momentum and final nucleon 3-momentum 

  // Set up dirac matrices via calling diracmatrices Fortran90 subroutine

   diracmatrices_(&xmn_in);

  // Calculate sigma via calling cc1 Fortran90 subroutine

   cc1_(&xq, &w, &wt, &xk, &xp, &ee0, &theta, &ig, &xsec, &nuphi);

  // Apply Jacobian to transform into the kPSFullNBody phase space
  double pLepFinal = fsLeptonP4.Vect().Mag();
  double ENf = fhad4.E();
  xsec *= 2. * ENf / ( kPi * pNf * pNf * pLepFinal );

  return xsec;
}
//____________________________________________________________________________
double FortranWrapperQELPXSec::Integral(const Interaction* in) const
{
  // We'll take care of this later. We need to use this function to get
  // total cross sections during spline generation.
  double integ = fXSecIntegrator->Integrate( this, in );
  return integ;
}
//____________________________________________________________________________
bool FortranWrapperQELPXSec::ValidProcess(const Interaction* interaction) const
{
  if ( interaction->TestBit(kISkipProcessChk) ) return true;

  const InitialState& init_state = interaction->InitState();
  const ProcessInfo&  proc_info  = interaction->ProcInfo();

  if ( !proc_info.IsQuasiElastic() ) return false;

  if ( !proc_info.IsEM() ) return false;

  return true;
}
//____________________________________________________________________________
void FortranWrapperQELPXSec::Configure(const Registry& config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void FortranWrapperQELPXSec::Configure(std::string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void FortranWrapperQELPXSec::LoadConfig(void)
{
  // Get access to the nuclear model for use in Pauli blocking, etc.
  fNuclModel = dynamic_cast< const NuclearModelI* >(
    this->SubAlg("NuclearModel") );
  assert( fNuclModel );

  // Load XSec Integrator
  fXSecIntegrator = dynamic_cast< const XSecIntegratorI* >(
    this->SubAlg("XSec-Integrator") );
  assert( fXSecIntegrator );
}
