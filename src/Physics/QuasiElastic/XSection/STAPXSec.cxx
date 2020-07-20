//_________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 For the class documentation see the corresponding header file.
*/
//_________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/HadronTensors/LabFrameHadronTensorI.h"
#include "Physics/HadronTensors/STAHadronTensorModel.h"
#include "Physics/QuasiElastic/XSection/STAPXSec.h"
#include "Physics/Multinucleon/XSection/MECUtils.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"

using namespace genie;

//_________________________________________________________________________
STAPXSec::STAPXSec() : XSecAlgorithmI("genie::STAPXSec")
{
}
//_________________________________________________________________________
STAPXSec::STAPXSec(string config)
  : XSecAlgorithmI("genie::STAPXSec", config)
{
}
//_________________________________________________________________________
STAPXSec::~STAPXSec()
{
}
//_________________________________________________________________________
double STAPXSec::XSec(const Interaction* interaction,
  KinePhaseSpace_t kps) const
{
  if ( !this->ValidProcess(interaction) ) return 0.;

  // Get the hadron tensor for the selected nuclide. Check the probe PDG code
  // to know whether to use the tensor for CC neutrino scattering or for
  // electron scattering
  int target_pdg = interaction->InitState().Tgt().Pdg();
  int probe_pdg = interaction->InitState().ProbePdg();

  HadronTensorType_t tensor_type = kHT_Undefined;
  if ( pdg::IsNeutrino(probe_pdg) || pdg::IsAntiNeutrino(probe_pdg) ) {
    tensor_type = kHT_QE_Full;
  }
  else {
    // If the probe is not a neutrino, assume that it's an electron
    tensor_type = kHT_QE_EM;
  }

  // The STA hadron tensors are defined using the same conventions as the
  // Valencia MEC (and SuSAv2-MEC) model, so we can use the same sort of tensor
  // object to describe them.
  const LabFrameHadronTensorI* tensor
    = dynamic_cast<const LabFrameHadronTensorI*>( fHadronTensorModel->GetTensor(target_pdg,
    tensor_type) );

  // If retrieving the tensor failed, complain and return zero
  if ( !tensor ) {
    LOG("STA", pWARN) << "Failed to load a hadronic tensor for the"
      " nuclide " << target_pdg;
    return 0.;
  }

  // Check that the input kinematical point is within the range in which hadron
  // tensors are known (for the chosen nuclear target)
  double Ev    = interaction->InitState().ProbeE( kRfLab );
  double Tl    = interaction->Kine().GetKV( kKVTl );
  double costl = interaction->Kine().GetKV( kKVctl );
  double ml    = interaction->FSPrimLepton()->Mass();
  double Q0    = 0.;
  double Q3    = 0.;

  genie::utils::mec::Getq0q3FromTlCostl( Tl, costl, Ev, ml, Q0, Q3 );

  double Q0min = tensor->q0Min();
  double Q0max = tensor->q0Max();
  double Q3min = tensor->qMagMin();
  double Q3max = tensor->qMagMax();
  if (Q0 < Q0min || Q0 > Q0max || Q3 < Q3min || Q3 > Q3max) {
    return 0.0;
  }

  // Saori's tensors for STA already have binding energy incorporated
  // in a way that requires no additional correction. Set the "Q-value"
  // to zero for this model.
  const double Q_value = 0.;

  // *** Enforce the global Q^2 cut (important for EM scattering) ***
  // Choose the appropriate minimum Q^2 value based on the interaction
  // mode (this is important for EM interactions since the differential
  // cross section blows up as Q^2 --> 0)
  double Q2min = genie::controls::kMinQ2Limit; // CC/NC limit
  if ( interaction->ProcInfo().IsEM() ) Q2min = genie::utils::kinematics
    ::electromagnetic::kMinQ2Limit; // EM limit

  double Q2 = Q3*Q3 - Q0*Q0;
  if ( Q2 < Q2min ) return 0.;

  // Compute the differential cross section using the hadron tensor
  // (dsigma / dTl / dctl in GeV^{-3) / atom)
  //double xsec = tensor->dSigma_dT_dCosTheta( interaction, Q_value );
  double xsec = tensor->dSigma_dT_dCosTheta_rosenbluth( interaction, Q_value );
  // TODO: revisit this!
  if ( std::isnan(xsec) ||  xsec < 0. ) xsec = 0.;
  //LOG("STA", pERROR) << "Tl = " << Tl << ", ctl = " << costl;
  //LOG("STA", pERROR) << "Q0 = " << Q0 << ", Q3 = " << Q3;
  //LOG("STA", pERROR) << "STA xsec = " << xsec;

  // Apply given overall scaling factor (for tuning)
  xsec *= fXSecScale;

  // Right now we don't have code ready for handling other phase spaces
  if ( kps != kPSTlctl ) {
    LOG("STA", pWARN)
      << "Doesn't support transformation from "
      << KinePhaseSpace::AsString( kPSTlctl ) << " to "
      << KinePhaseSpace::AsString( kps );
    xsec = 0.;
  }

  return xsec;
}
//_________________________________________________________________________
double STAPXSec::Integral(const Interaction* interaction) const
{
  double xsec = fXSecIntegrator->Integrate( this, interaction );
  return xsec;
}
//_________________________________________________________________________
bool STAPXSec::ValidProcess(const Interaction* interaction) const
{
  if ( interaction->TestBit(kISkipProcessChk) ) return true;

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();

  if ( !proc_info.IsQuasiElastic() ) return false;

  // The STA calculation is only appropriate for complex nuclear targets, not
  // free nucleons.
  if ( !init_state.Tgt().IsNucleus() ) return false;

  int  nuc = init_state.Tgt().HitNucPdg();
  int  nu  = init_state.ProbePdg();

  bool isP   = pdg::IsProton( nuc );
  bool isN   = pdg::IsNeutron( nuc );
  bool isnu  = pdg::IsNeutrino( nu );
  bool isnub = pdg::IsAntiNeutrino( nu );
  bool is_chgl = pdg::IsChargedLepton( nu );

  bool prcok = ( proc_info.IsWeakCC() && ( (isP && isnub) || (isN && isnu)) )
    || ( proc_info.IsEM() && is_chgl && (isP || isN) );
  if ( !prcok ) return false;

  return true;
}
//_________________________________________________________________________
void STAPXSec::Configure(const Registry& config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void STAPXSec::Configure(std::string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//_________________________________________________________________________
void STAPXSec::LoadConfig(void)
{
  // Cross section scaling factor
  GetParamDef( "XSecScale", fXSecScale, 1.0 ) ;

  // Load hadron tensor model
  fHadronTensorModel = dynamic_cast< const STAHadronTensorModel* >(
    this->SubAlg("HadronTensorAlg") );
  assert( fHadronTensorModel );

  // Load XSec Integrator
  fXSecIntegrator = dynamic_cast<const XSecIntegratorI *>(
    this->SubAlg("XSec-Integrator") );
  assert( fXSecIntegrator );
}
