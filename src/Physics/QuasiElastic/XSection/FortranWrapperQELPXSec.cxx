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

#include "Physics/QuasiElastic/XSection/FortranWrapperQELPXSec.h"
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

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;
using namespace genie::utils;

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
  // This is the main "missing piece" for interfacing with Noemi's code.
  // Given an Interaction object which will be pre-loaded with the correct
  // variables, we need to extract the inputs to Noemi's function, call that
  // function, and then return a differential cross section in natural units.
  // You can ignore the kps variable for now. Just work with the output
  // of Noemi's cc1 subroutine here.
  double xsec = 0.;
  return xsec;
}
//____________________________________________________________________________
double FortranWrapperQELPXSec::Integral(const Interaction* in) const
{
  // We'll take care of this later. We need to use this function to get
  // total cross sections during spline generation.
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
}
