//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Steven Gardiner <gardiner \at fnal.gov>
         Fermi National Accelerator Laboratory - Feb 26, 2019

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/RefFrame.h"
#include "Physics/QuasiElastic/XSection/NewQELXSec.h"

#include "Physics/NuclearState/NuclearModelI.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Numerical/GSLUtils.h"
#include "Physics/Common/VertexGenerator.h"
#include "Physics/NuclearState/PauliBlocker.h"
#include "Physics/NuclearState/NuclearModel.h"
#include "Physics/NuclearState/NuclearModelMap.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::utils::gsl;

//____________________________________________________________________________
NewQELXSec::NewQELXSec() : XSecIntegratorI("genie::NewQELXSec")
{

}
//____________________________________________________________________________
NewQELXSec::NewQELXSec(std::string config) : XSecIntegratorI("genie::NewQELXSec", config)
{

}
//____________________________________________________________________________
double NewQELXSec::Integrate(const XSecAlgorithmI* model, const Interaction* in) const
{
  LOG("NewQELXSec",pDEBUG) << "Beginning integrate";
  if ( !model->ValidProcess(in) ) return 0.;

  const KPhaseSpace& kps = in->PhaseSpace();
  if ( !kps.IsAboveThreshold() ) {
     LOG("NewQELXSec", pDEBUG)  << "*** Below energy threshold";
     return 0.;
  }

  Interaction* interaction = new Interaction( *in );
  interaction->SetBit( kISkipProcessChk );
  //interaction->SetBit( kISkipKinematicChk );

  // We're doing the integration over the nucleon momentum
  // distribution ourselves, so we don't want to apply the
  // nuclear suppression factor. To turn it off, we'll
  // set the "assume free nucleon" flag.
  interaction->SetBit( kIAssumeFreeNucleon );

  const NuclearModelI* nucl_model = dynamic_cast<const NuclearModelI*>(
    model->SubAlg("IntegralNuclearModel") );
  assert( nucl_model );

  AlgFactory* algf = AlgFactory::Instance();
  const VertexGenerator* vtx_gen = dynamic_cast<const VertexGenerator*>(
    algf->GetAlgorithm(fVertexGenID) );
  assert( vtx_gen );

  // Determine the appropriate integration mode and binding energy mode to use.
  // The defaults given here are for the case of a free nucleon.
  FullQELdXSecMode_t integration_mode = kFreeNucleonVars;
  QELEvGen_BindingMode_t bind_mode = kOnShell;
  const Target& tgt = interaction->InitState().Tgt();
  if ( tgt.IsNucleus() ) {
    // Use a different approach for a nuclear target
    integration_mode = this->IntegrationModeFromNuclearModel( *nucl_model, tgt );

    std::string bind_mode_str = model->GetConfig()
      .GetString( "IntegralNucleonBindingMode" );
    bind_mode = genie::utils::StringToQELBindingMode( bind_mode_str );
  }

  ROOT::Math::IBaseFunctionMultiDim* func = new utils::gsl::FullQELdXSec(model,
    interaction, integration_mode, fPauliBlock, fPauliBlockID, bind_mode, fVertexGenID);
  ROOT::Math::IntegrationMultiDim::Type ig_type =
    utils::gsl::IntegrationNDimTypeFromString( fGSLIntgType );

  double abstol = 1e-16; // We mostly care about relative tolerance
  ROOT::Math::IntegratorMultiDim ig(*func, ig_type, abstol, fGSLRelTol, fGSLMaxEval);

  // TODO: check min/max range dependence on hit nucleon radius (defaults to zero)
  Range1D_t pNi_lim( nucl_model->MinMomentum(tgt), nucl_model->MaxMomentum(tgt) );
  Range1D_t E_lim( nucl_model->MinRemovalEnergy(tgt), nucl_model->MaxRemovalEnergy(tgt) );
  Range1D_t r_lim( 0., vtx_gen->RMax(tgt.A()) );
  Range1D_t cos_theta_Ni_lim( -1., 1. );
  Range1D_t phi_Ni_lim( 0., 2.*kPi );
  Range1D_t cos_theta_0_lim( -1., 1. );
  Range1D_t phi_0_lim( 0., 2.*kPi );

  // This can be done more cleanly in C++11 with std::vector and initializer
  // lists
  double kine_min[7] = { cos_theta_0_lim.min, phi_0_lim.min, pNi_lim.min,
    cos_theta_Ni_lim.min, phi_Ni_lim.min, 0., 0. };
  double kine_max[7] = { cos_theta_0_lim.max, phi_0_lim.max, pNi_lim.max,
    cos_theta_Ni_lim.max, phi_Ni_lim.max, 0., 0. };

  // Initialize the extra elements if needed for a particular integration mode
  if ( integration_mode == kLocalFermiGasVars ) {
    kine_min[5] = r_lim.min;
    kine_max[5] = r_lim.max;
  }
  else if ( integration_mode == kSpectralFuncVars ) {
    kine_min[5] = E_lim.min;
    kine_max[5] = E_lim.max;
  }
  else if ( integration_mode == kAllQELVars ) {
    kine_min[5] = r_lim.min;
    kine_max[5] = r_lim.max;

    kine_min[6] = E_lim.min;
    kine_max[6] = E_lim.max;
  }

  double xsec = ig.Integral(kine_min, kine_max) * (1E-38 * units::cm2);

  delete func;
  delete interaction;

  return xsec;
}
//____________________________________________________________________________
void NewQELXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NewQELXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NewQELXSec::LoadConfig(void)
{
  // Get GSL integration type & relative tolerance
  GetParamDef( "gsl-integration-type", fGSLIntgType, std::string("adaptive") ) ;
  GetParamDef( "gsl-relative-tolerance", fGSLRelTol, 1e-2 ) ;
  int max;
  GetParamDef( "gsl-max-eval", max, 500000 ) ;
  fGSLMaxEval  = static_cast<unsigned int>( max );

  RgAlg vertexGenID;
  GetParamDef( "VertexGenAlg", vertexGenID, RgAlg("genie::VertexGenerator", "Default") );
  fVertexGenID = AlgId( vertexGenID );

  GetParamDef( "DoPauliBlocking", fPauliBlock, true );

  RgAlg pauliBlockID;
  GetParamDef( "PauliBlockerAlg", pauliBlockID, RgAlg("genie::PauliBlocker", "Default") );
  fPauliBlockID = AlgId( pauliBlockID );
}
//____________________________________________________________________________
FullQELdXSecMode_t NewQELXSec::IntegrationModeFromNuclearModel(
  const NuclearModelI& nucl_model, const Target& tgt) const
{
  std::string alg_name = nucl_model.Id().Name();
  if (alg_name == "genie::FGMBodekRitchie") return kFermiGasVars;
  if (alg_name == "genie::LocalFGM") return kLocalFermiGasVars;
  if (alg_name == "genie::SpectralFunc") return kSpectralFuncVars;
  else if (alg_name == "genie::NuclearModelMap") {
    const NuclearModelMap& nmm = dynamic_cast<const NuclearModelMap&>( nucl_model );
    NuclearModel_t type = nmm.ModelType( tgt );
    if ( type == kNucmFermiGas ) return kFermiGasVars;
    else if ( type == kNucmLocalFermiGas ) return kLocalFermiGasVars;
    else if ( type == kNucmSpectralFunc ) return kSpectralFuncVars;
    else {
      LOG("NewQELXSec", pFATAL) << "Unrecognized nuclear model type "
        << type << " returned by genie::NuclearModelMap within NewQELXSec"
        << "::IntegrationModeFromNuclearModel()";
      std::exit(1);
    }
  }
  else {
    LOG("NewQELXSec", pFATAL) << "Unrecognized nuclear model \""
      << alg_name << "\" passed to NewQELXSec::IntegrationModeFromNuclearModel()";
    std::exit(1);
  }

  // This should never execute
  return kFreeNucleonVars;
}


genie::utils::gsl::FullQELdXSec::FullQELdXSec(const XSecAlgorithmI* xsec_model,
  const Interaction* interaction, FullQELdXSecMode_t mode, bool do_Pauli_blocking,
  const AlgId& pauli_blocker_ID, QELEvGen_BindingMode_t binding_mode,
  const AlgId& vertex_gen_ID) : fXSecModel( xsec_model ),
  fInteraction( new Interaction(*interaction) ), fIntegrationMode( mode ),
  fDoPauliBlocking( do_Pauli_blocking ), fHitNucleonBindingMode( binding_mode )
{
  fNuclModel = dynamic_cast<const NuclearModelI*>( fXSecModel->SubAlg("IntegralNuclearModel") );
  assert( fNuclModel );

  AlgFactory* algf = AlgFactory::Instance();
  fPauliBlocker = dynamic_cast<const PauliBlocker*>( algf->GetAlgorithm(pauli_blocker_ID) );
  assert( fPauliBlocker );

  fVertexGenerator = dynamic_cast<const VertexGenerator*>( algf->GetAlgorithm(vertex_gen_ID) );
  assert( fVertexGenerator );
}

genie::utils::gsl::FullQELdXSec::~FullQELdXSec()
{
  delete fInteraction;
}

ROOT::Math::IBaseFunctionMultiDim* genie::utils::gsl::FullQELdXSec::Clone(void) const
{
  return new FullQELdXSec(fXSecModel, fInteraction, fIntegrationMode, fDoPauliBlocking,
    fPauliBlocker->Id(), fHitNucleonBindingMode, fVertexGenerator->Id());
}

unsigned int genie::utils::gsl::FullQELdXSec::NDim(void) const
{
  if ( fIntegrationMode == kFreeNucleonVars ) return 2;
  else if ( fIntegrationMode == kFermiGasVars ) return 5;
  else if ( fIntegrationMode == kLocalFermiGasVars ) return 6;
  else if ( fIntegrationMode == kSpectralFuncVars ) return 6;
  else if ( fIntegrationMode == kAllQELVars ) return 7;
  else {
    LOG("FullQELdXSec", pFATAL) << "Unrecognized FullQELdXSecMode_t value \""
      << fIntegrationMode << "\" encountered in FullQELdXSec::NDim()";
    std::exit(1);
  }
}

double genie::utils::gsl::FullQELdXSec::DoEval(const double* xin) const
{
  // *** Input for "FreeNucleonVars" mode ***
  //
  // element 0: "cos_theta0" = Cosine of theta0, the angle between the COM frame
  //                           3-momentum of the outgoing lepton and the COM frame velocity
  //                           as measured in the laboratory frame
  // element 1: "phi_theta0" = Azimuthal angle of the COM frame 3-momentum of the
  //                           outgoing lepton measured with respect to the COM frame
  //                           velocity as measured in the laboratory frame
  //
  // *** For all other modes, add ***
  //
  // element 2: "pNi" = Momentum of initial hit nucleon 3-momentum (GeV)
  // element 3: "cos_theta_Ni" = Cosine of the angle between the initial hit nucleon
  //                             3-momentum and the probe direction
  // element 4: "phi_Ni" = Azimuthal angle of the initial hit nucleon 3-momentum
  //                       with respect to probe direction
  //
  // *** For "LocalFermiGasVars" mode, also add ***
  //
  // element 5: "hitNucleonRadius" = Radial position of the initial struck nucleon
  //
  // *** For "SpectralFuncVars" mode, also add ***
  //
  // element 5: "E_remove" = Removal energy for the initial struck nucleon
  //
  // *** For "AllQELVars" mode, also add ***
  //
  // element 5: "hitNucleonRadius" = Radial position of the initial struck nucleon
  // element 6: "E_remove" = Removal energy for the initial struck nucleon

  double cos_theta0 = xin[0];
  double phi0 = xin[1];
  double pNi = xin[2];
  double cos_theta_Ni = xin[3];
  double sin_theta_Ni = std::sqrt(std::max(0., 1. - std::pow(cos_theta_Ni, 2)));
  double phi_Ni = xin[4];

  // Set the radius for the hit nucleon and the removal energy as needed
  const Target& tgt = fInteraction->InitState().Tgt();
  double hitNucleonRadius = 0.;
  double E_remove = 0.;
  if ( fIntegrationMode == kFermiGasVars ) {
   // In (global) Fermi gas mode, the removal energy is constant, so
   // just set it to its minimum value (equal to the maximum)
   E_remove = fNuclModel->MinRemovalEnergy(tgt);
  }
  else if ( fIntegrationMode == kLocalFermiGasVars ) {
    hitNucleonRadius = xin[5];
    E_remove = fNuclModel->MinRemovalEnergy(tgt, hitNucleonRadius);
  }
  else if ( fIntegrationMode == kSpectralFuncVars ) {
    E_remove = xin[5];
  }
  else if ( fIntegrationMode == kAllQELVars ) {
    hitNucleonRadius = xin[5];
    E_remove = xin[6];
  }
  else if ( fIntegrationMode != kFreeNucleonVars ) {
    LOG("FullQELdXSec", pFATAL) << "Unrecognized FullQELdXSecMode_t value \""
      << fIntegrationMode << "\" encountered in FullQELdXSec::DoEval()";
    std::exit(1);
  }

  // Tell the nuclear model about the current value of the removal energy. It needs
  // to be set in the nuclear model so that it will be handled correctly by
  // genie::utils::ComputeFullQELPXSec()
  fNuclModel->SetRemovalEnergy( E_remove );

  // Construct the initial hit nucleon's 3-momentum
  TVector3 p3Ni(0., 0., 0.);
  if ( fIntegrationMode != kFreeNucleonVars ) {
    p3Ni.SetX( pNi*sin_theta_Ni*std::cos(phi_Ni) );
    p3Ni.SetY( pNi*sin_theta_Ni*std::sin(phi_Ni) );
    p3Ni.SetZ( pNi*cos_theta_Ni );
  }

  // Set the 3-momentum of the nucleon to the specified value
  fNuclModel->SetMomentum3( p3Ni );

  // Dummy storage for the binding energy of the hit nucleon
  double dummy_Eb = 0.;

  // TODO: min_angle_EM!!!!
  double min_angle_EM = 0.;

  // Compute the full differential cross section
  double xsec = genie::utils::ComputeFullQELPXSec(fInteraction, fNuclModel,
    fXSecModel, cos_theta0, phi0, dummy_Eb, fHitNucleonBindingMode, min_angle_EM, true);

  // If we're integrating the cross section for a free nucleon, then just return
  // what we have so far (no additional factors needed)
  if ( fIntegrationMode == kFreeNucleonVars ) return xsec;

  // ComputeFullQELPXSec() sets the final nucleon's lab frame 4-momentum via
  // interaction->KinePtr()->SetHadSystP4(), so we can now check for Pauli
  // blocking

  // If the final nucleon would be Pauli blocked, then return zero immediately
  if ( fDoPauliBlocking ) {
    double kF = fPauliBlocker->GetFermiMomentum(tgt, fInteraction->RecoilNucleonPdg(),
      hitNucleonRadius);
    double pNf = fInteraction->Kine().HadSystP4().P();
    if ( pNf < kF ) return 0.;
  }

  // PDG code for the initial hit nucleon
  int pdg_Ni = tgt.HitNucPdg();
  // Number of active nucleons in the target
  int num_active_nucleons = genie::pdg::IsProton( pdg_Ni ) ? tgt.Z() : tgt.N();

  // Probability density for sampling the current initial nucleon momentum
  double prob_density_pNi_E = fNuclModel->ProbDensity(pNi, E_remove, tgt, hitNucleonRadius);

  xsec *= std::pow(pNi, 2) * num_active_nucleons * prob_density_pNi_E;

  // For radially-dependent nuclear models, also multiply by the probability density for
  // sampling the given hit nucleon radius
  if ( fIntegrationMode == kLocalFermiGasVars || fIntegrationMode == kAllQELVars ) {
    double prob_density_hNR = fVertexGenerator->ProbDensity( tgt.A(), hitNucleonRadius );
    xsec *= prob_density_hNR;
  }

  return xsec;
}
