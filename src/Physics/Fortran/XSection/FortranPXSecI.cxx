//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Steven Gardiner <gardiner \at fnal.gov>
         Fermi National Accelerator Laboratory

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"

#include "Physics/Fortran/XSection/FortranPXSecI.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::utils;

//____________________________________________________________________________
FortranPXSecI::FortranPXSecI()
  : XSecAlgorithmI("genie::FortranPXSecI")
{

}

//____________________________________________________________________________
FortranPXSecI::FortranPXSecI(const std::string& config)
  : XSecAlgorithmI("genie::FortranPXSecI", config)
{

}

//____________________________________________________________________________
FortranPXSecI::~FortranPXSecI()
{

}

//____________________________________________________________________________
double FortranPXSecI::XSec(const Interaction* interaction,
  KinePhaseSpace_t kps) const
{

  //if (! this -> ValidProcess    (interaction) ) return 0.;
  //if (! this -> ValidKinematics (interaction) ) return 0.;

  double xsec = 0.;

  std::cout << "C++ side: kps == " << kps << '\n';
  compute_my_diff_xsec(interaction, kps, &xsec);

  return xsec;
}
//____________________________________________________________________________
double FortranPXSecI::Integral(const Interaction* in) const
{
  // The hadron tensor is assumed to already be averaged over
  // the initial nucleon momenta, so just integrate
  // the differential cross section in the normal way.
  /// \todo Revisit this
  return fXSecIntegrator->Integrate(this, in);
}

//____________________________________________________________________________
bool FortranPXSecI::ValidProcess(const Interaction * interaction) const
{
  if ( interaction->TestBit(kISkipProcessChk) ) return true;

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();

  if ( !proc_info.IsQuasiElastic() ) return false;

  int hit_nucleon_pdg = init_state.Tgt().HitNucPdg();
  int probe_pdg = init_state.ProbePdg();

  bool is_P      = pdg::IsProton(hit_nucleon_pdg);
  bool is_N      = pdg::IsNeutron(hit_nucleon_pdg);
  bool is_nu     = pdg::IsNeutrino(probe_pdg);
  bool is_nubar  = pdg::IsAntiNeutrino(probe_pdg);
  bool is_e      = pdg::IsElectron(probe_pdg);

  if ( proc_info.IsWeakCC() && ((is_P && is_nubar) || (is_N && is_nu)) ) {
    // CCQE
    return true;
  }
  else if ( proc_info.IsWeakNC() && (is_nu || is_nubar) ) {
    // NCQE
    return true;
  }
  else if ( proc_info.IsEM() && is_e ) {
    // EMQE
    return true;
  }

  return false;
}
//____________________________________________________________________________
void FortranPXSecI::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void FortranPXSecI::Configure(std::string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void FortranPXSecI::LoadConfig(void)
{
   // load XSec Integrator
  fXSecIntegrator = dynamic_cast<const XSecIntegratorI*>(
    this->SubAlg("XSec-Integrator") );
  assert(fXSecIntegrator);
}
