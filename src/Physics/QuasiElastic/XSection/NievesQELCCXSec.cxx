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
#include "Physics/QuasiElastic/XSection/NievesQELCCXSec.h"

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
#include "Physics/NuclearState/NuclearModel.h"
#include "Physics/NuclearState/NuclearModelMap.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::utils::gsl;

//____________________________________________________________________________
NievesQELCCXSec::NievesQELCCXSec() : XSecIntegratorI("genie::NievesQELCCXSec")
{

}
//____________________________________________________________________________
NievesQELCCXSec::NievesQELCCXSec(std::string config) : XSecIntegratorI("genie::NievesQELCCXSec", config)
{

}
//____________________________________________________________________________
double NievesQELCCXSec::Integrate(const XSecAlgorithmI* model, const Interaction* in) const
{
  LOG("NievesQELCCXSec",pDEBUG) << "Beginning integrate";
  if ( !model->ValidProcess(in) ) return 0.;

  Interaction* interaction = new Interaction( *in );
  interaction->SetBit( kISkipProcessChk );
  interaction->SetBit( kISkipKinematicChk );

  AlgFactory* algf = AlgFactory::Instance();
  const VertexGenerator* vtx_gen = dynamic_cast<const VertexGenerator*>(
    algf->GetAlgorithm(fVertexGenID) );
  assert( vtx_gen );

  ROOT::Math::IBaseFunctionMultiDim* func = new utils::gsl::NievesQELdXSec(model,
    interaction);
  ROOT::Math::IntegrationMultiDim::Type ig_type =
    utils::gsl::IntegrationNDimTypeFromString( fGSLIntgType );

  double abstol = 1e-16; // We mostly care about relative tolerance
  ROOT::Math::IntegratorMultiDim ig(*func, ig_type, abstol, fGSLRelTol, fGSLMaxEval);

  double probeE = interaction->InitState().ProbeE( kRfLab );
  double ml = interaction->FSPrimLepton()->Mass();

  // TODO: This solution is fragile. If the formula used by VertexGenerator
  // changes, then this one will need to change too. Switch to using
  // a common function to get Rmax for both.
  double R0 = 1.4;
  //model->GetParamDef("NUCL-R0", R0, 1.4);
  double Rmax = 3. * R0 * std::pow(interaction->InitState().Tgt().A(), 1./3.);

  // Integration ranges for the lepton scattering cosine, the lepton kinetic
  // energy, and the hit nucleon radius
  Range1D_t ctl_lim( -1., 1. );
  Range1D_t Tl_lim( 0., probeE - ml);
  Range1D_t r_lim( 0., Rmax);

  double kine_min[3] = { ctl_lim.min, Tl_lim.min, r_lim.min };
  double kine_max[3] = { ctl_lim.max, Tl_lim.max, r_lim.max };

  //// Select a new position for the initial hit nucleon (needed for the local
  //// Fermi gas model, but other than slowing things down a bit, it doesn't
  //// hurt to do this for other models)
  //TVector3 vertex_pos = vtx_gen->GenerateVertex( interaction, tgt->A() );
  //double radius = vertex_pos.Mag();
  //tgt->SetHitNucPosition( radius );

  // The initial state variables have all been defined, so integrate over
  // the final lepton angles.
  double total_xsec = ig.Integral(kine_min, kine_max);

  delete func;
  delete interaction;

  return total_xsec;
}
//____________________________________________________________________________
void NievesQELCCXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NievesQELCCXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NievesQELCCXSec::LoadConfig(void)
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
}

genie::utils::gsl::NievesQELdXSec::NievesQELdXSec(const XSecAlgorithmI* xsec_model,
  const Interaction* interaction) : fXSecModel( xsec_model ),
  fInteraction( new Interaction(*interaction) )
{
  fInteraction->SetBit( kISkipProcessChk );
  fInteraction->SetBit( kISkipKinematicChk );
}

genie::utils::gsl::NievesQELdXSec::~NievesQELdXSec()
{
  delete fInteraction;
}

ROOT::Math::IBaseFunctionMultiDim* genie::utils::gsl::NievesQELdXSec::Clone(void) const
{
  return new NievesQELdXSec(fXSecModel, fInteraction);
}

unsigned int genie::utils::gsl::NievesQELdXSec::NDim(void) const
{
  return 3;
}

double genie::utils::gsl::NievesQELdXSec::DoEval(const double* xin) const
{
  // Elements of "xin"
  //
  // element 0: "ctl" = Scattering cosine of outgoing lepton
  //
  // element 1: "Tl"  = Kinetic energy of outgoing lepton
  //
  // element 2: "r"   = Radial position of initial struck nucleon

  double ctl = xin[0];
  double Tl = xin[1];
  double r = xin[2];

  fInteraction->InitState().TgtPtr()->SetHitNucPosition( r );
  fInteraction->KinePtr()->SetKV( kKVctl, ctl );
  fInteraction->KinePtr()->SetKV( kKVTl, Tl );

  double xsec = fXSecModel->XSec( fInteraction, kPSTlctl );

  return xsec;
}
