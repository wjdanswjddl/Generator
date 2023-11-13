//____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Steven Gardiner <gardiner \at fnal.gov>
 Fermi National Accelerator Laboratory
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
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/GSLUtils.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/Range1.h"

#include "Physics/Resonance/XSection/SuSAv2InelXSec.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::utils::gsl;

//____________________________________________________________________________
SuSAv2InelXSec::SuSAv2InelXSec() : XSecIntegratorI( "genie::SuSAv2InelXSec" )
{

}

//____________________________________________________________________________
SuSAv2InelXSec::SuSAv2InelXSec( std::string config )
  : XSecIntegratorI( "genie::SuSAv2InelXSec", config )
{

}

//____________________________________________________________________________
double SuSAv2InelXSec::Integrate( const XSecAlgorithmI* model,
  const Interaction* in ) const
{
  LOG( "SuSAv2InelXSec",pDEBUG ) << "Beginning integrate";
  if ( !model->ValidProcess(in) ) return 0.;

  // Copy the input interaction so we can adjust the process check bit without
  // modifying the original
  Interaction* interaction = new Interaction( *in );
  interaction->SetBit( kISkipProcessChk );
  //interaction->SetBit( kISkipKinematicChk );

  // Create the helper object to manage the integration. This creates
  // a copy of our copy of the input interaction
  utils::gsl::SuSAv2IneldXSec* func = new utils::gsl::SuSAv2IneldXSec(model,
    interaction );

  ROOT::Math::IntegrationMultiDim::Type ig_type =
    utils::gsl::IntegrationNDimTypeFromString( fGSLIntgType );

  // Switch to using the copy of the interaction in the integrator rather than
  // the copy that we made in this function (which we delete since it is no
  // longer needed)
  delete interaction;
  interaction = func->GetInteractionPtr();

  double abstol = 1e-16; // We mostly care about relative tolerance
  ROOT::Math::IntegratorMultiDim ig( *func, ig_type, abstol, fGSLRelTol,
    fGSLMaxEval );

  // Get the lab-frame probe energy from the input interaction
  double Ev = in->InitState().ProbeE( kRfLab );

  // Get the final-state lepton mass from the input interaction
  double ml = in->FSPrimLepton()->Mass();

  // Set the integration ranges for the kinematic variables used in the
  // SuSAv2InelPXSec triple-differential cross section. The hadronic invariant
  // mass W is Lorentz-invariant, while the other variables are evaluated in
  // the laboratory frame. The integration ranges are set conservatively for
  // guaranteed full coverage, even if the kinematics are actually more
  // narrowly constrained in reality.
  Range1D_t W_lim( kProtonMass, Ev - kProtonMass );
  Range1D_t cos_theta_lim( -1., 1. );
  Range1D_t Tl_lim( 0., Ev - ml );

  double kine_min[3] = { W_lim.min, cos_theta_lim.min, Tl_lim.min };
  double kine_max[3] = { W_lim.max, cos_theta_lim.max, Tl_lim.max };

  // We're ready. Integrate over all three dimensions.
  double xsec_total = ig.Integral( kine_min, kine_max );
  delete func;

  // The triple-differential cross section also has uniform azimuthal lepton
  // angle dependence. Throw in a factor of two pi to account for this, then
  // return the result.
  xsec_total *= 2. * kPi;
  return xsec_total;

}

//____________________________________________________________________________
void SuSAv2InelXSec::Configure( const Registry & config )
{
  Algorithm::Configure( config );
  this->LoadConfig();
}

//____________________________________________________________________________
void SuSAv2InelXSec::Configure( string config )
{
  Algorithm::Configure( config );
  this->LoadConfig();
}

//____________________________________________________________________________
void SuSAv2InelXSec::LoadConfig()
{
  // Get GSL integration type & relative tolerance
  GetParamDef( "gsl-integration-type", fGSLIntgType, std::string("adaptive") ) ;
  GetParamDef( "gsl-relative-tolerance", fGSLRelTol, 1e-2 );

  int max;
  GetParamDef( "gsl-max-eval", max, 500000 );
  fGSLMaxEval = static_cast<unsigned int>( max );
}

//____________________________________________________________________________
genie::utils::gsl::SuSAv2IneldXSec::SuSAv2IneldXSec(
  const XSecAlgorithmI* xsec_model, const Interaction* interaction )
  : fXSecModel( xsec_model ), fInteraction( new Interaction(*interaction) )
{
}

//____________________________________________________________________________
genie::utils::gsl::SuSAv2IneldXSec::~SuSAv2IneldXSec()
{
  delete fInteraction;
}

//____________________________________________________________________________
Interaction* genie::utils::gsl::SuSAv2IneldXSec::GetInteractionPtr()
{
  return fInteraction;
}

//____________________________________________________________________________
const Interaction& genie::utils::gsl::SuSAv2IneldXSec::GetInteraction() const
{
  return *fInteraction;
}

//____________________________________________________________________________
ROOT::Math::IBaseFunctionMultiDim*
  genie::utils::gsl::SuSAv2IneldXSec::Clone() const
{
  return new SuSAv2IneldXSec( fXSecModel, fInteraction );
}

//____________________________________________________________________________
unsigned int genie::utils::gsl::SuSAv2IneldXSec::NDim() const
{
  return 3u;
}

//____________________________________________________________________________
double genie::utils::gsl::SuSAv2IneldXSec::DoEval( const double* xin ) const
{
  // Elements of "xin"
  //
  // element 0: "W" = Hadronic invariant mass (GeV)
  //
  // element 1: "ctl" = Scattering cosine of the final-state lepton as measured
  //                    in the laboratory frame
  //
  // element 2: "Tl" = Kinetic energy of the final-state lepton, as measured
  //                   in the laboratory frame

  // Retrieve the kinematic variables and store them in the owned interaction
  double W = xin[0];
  double ctl = xin[1];
  double Tl = xin[2];

  fInteraction->KinePtr()->SetW( W );
  fInteraction->KinePtr()->SetKV( genie::kKVctl, ctl );
  fInteraction->KinePtr()->SetKV( genie::kKVTl, Tl );

  // Compute the full differential cross section
  double xsec = fXSecModel->XSec( fInteraction, kPSTlctl );
  return xsec;
}
