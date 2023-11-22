//____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/RunningThreadInfo.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/ParticleData/BaryonResonance.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/Resonance/EventGen/SuSAv2InelGenerator.h"

//___________________________________________________________________________
genie::SuSAv2InelGenerator::SuSAv2InelGenerator()
  : genie::KineGeneratorWithCache( "genie::SuSAv2InelGenerator" )
{
}

//___________________________________________________________________________
genie::SuSAv2InelGenerator::SuSAv2InelGenerator( std::string config )
  : genie::KineGeneratorWithCache( "genie::SuSAv2InelGenerator", config )
{
}

//___________________________________________________________________________
genie::SuSAv2InelGenerator::~SuSAv2InelGenerator()
{
}

//___________________________________________________________________________
void genie::SuSAv2InelGenerator::ProcessEventRecord( GHepRecord* evrec ) const
{
  // Access the cross section algorithm for running thread (presumably the
  // SuSAv2 inelastic one)
  genie::RunningThreadInfo* rt_info = RunningThreadInfo::Instance();
  const genie::EventGeneratorI* evg = rt_info->RunningThread();
  fXSecModel = evg->CrossSectionAlg();

  this->SelectLeptonKinematics( evrec );

  // TODO: construct the hadronic final state here
}

//___________________________________________________________________________
void genie::SuSAv2InelGenerator::SelectLeptonKinematics( GHepRecord* evrec )
  const
{
  // Get access to the interaction from within the event record
  genie::Interaction* interaction = evrec->Summary();
  interaction->SetBit( genie::kISkipProcessChk );

  // Check whether this is an EM or weak (CC only for now) process
  bool is_em = interaction->ProcInfo().IsEM();

  // Choose the appropriate minimum Q^2 value based on the interaction
  // mode (this is important for EM interactions since the differential
  // cross section blows up as Q^2 --> 0)
  double Q2min = genie::controls::kMinQ2Limit; // CC/NC limit
  if ( is_em ) Q2min = genie::utils::kinematics::electromagnetic::kMinQ2Limit;

  // Get access to the random number generators
  genie::RandomGen* rnd = RandomGen::Instance();

  // Retrieve the probe 4-momentum in the lab frame. Delete the cloned version
  // to avoid a memory leak.
  TLorentzVector* temp_p4v = interaction->InitState().GetProbeP4();
  TLorentzVector p4v( *temp_p4v );
  delete temp_p4v;

  // Get the lab-frame probe energy from the input interaction
  double Ev = p4v.E();

  // Get the final-state lepton mass from the input interaction
  double ml = interaction->FSPrimLepton()->Mass();

  // Set the kinematic limits for sampling
  // TODO: revisit these and assign better values if needed
  double W_min = genie::constants::kNeutronMass + genie::constants::kPionMass;
  double W_max = std::max( Ev - W_min , W_min );

  double cth_min = -1.;
  double cth_max =  1.;

  double Tl_min = 0.;
  double Tl_max = Ev - ml;

  // Calculate the maximum differential cross section or retrieve it
  double xsec_max = this->MaxXSec( evrec );

  double W, cth, Tl, Q2, xsec;
  TLorentzVector p4l;
  unsigned int iter = 0;
  bool accept = false;
  while ( true ) {
    ++iter;
    if( iter > genie::controls::kRjMaxIterations ) {
      LOG( "SuSAv2Inel", pWARN )
        << "*** Could not select a valid (W, cth, Tl) tuple after "
        << iter << " iterations";
      evrec->EventFlags()->SetBitNumber( genie::kKineGenErr, true );
      genie::exceptions::EVGThreadException exception;
      exception.SetReason( "Couldn't select kinematics" );
      exception.SwitchOnFastForward();
      throw exception;
    }

    // Currently sampled values of the kinematic variables of interest
    W = W_min + ( W_max - W_min ) * rnd->RndKine().Rndm();
    cth = cth_min + ( cth_max - cth_min ) * rnd->RndKine().Rndm();
    Tl = Tl_min + ( Tl_max - Tl_min ) * rnd->RndKine().Rndm();

    LOG( "SuSAv2Inel", pDEBUG ) << "Trying: W = " << W
      << ", cth = " << cth << ", Tl = " << Tl;

    // Update the interaction object with the sampled values
    interaction->KinePtr()->SetW( W );
    interaction->KinePtr()->SetKV( genie::kKVctl, cth );
    interaction->KinePtr()->SetKV( genie::kKVTl, Tl );

    // Compute the differential cross section for the current kinematics
    xsec = fXSecModel->XSec( interaction, genie::kPSTlctl );

    LOG( "SuSAv2Inel", pDEBUG ) << "xsec = " << xsec;

    // Do the usual check for rejection sampling of the kinematics
    double y = xsec_max * rnd->RndKine().Rndm();
    accept = (y < xsec);

    // Double-check that the precomputed maximum differential cross section
    // wasn't exceeded by the current one (otherwise we have a rejection
    // sampling problem)
    this->AssertXSecLimits( interaction, xsec, xsec_max );

    // If the generated kinematics are accepted, then double-check the Q^2
    // limit before breaking out of the rejection sampling loop. Skip the
    // check if we're not going to accept the variables thrown to save time.
    if ( accept ) {

      // Select an azimuthal angle for the outgoing lepton (uniformly
      // distributed)
      double phi = 2. * genie::constants::kPi * rnd->RndKine().Rndm();

      // Build the outgoing lepton 4-momentum in the lab frame
      double El = ml + Tl;
      double pl = std::sqrt( std::max(0., El*El - ml*ml) );
      double sth = std::sqrt( std::max(0., 1. - cth*cth) );

      double pxl = pl * sth * std::cos( phi );
      double pyl = pl * sth * std::sin( phi );
      double pzl = pl * cth;

      // Rotate lepton 3-momentum vector so that the prior z direction is along
      // the probe direction
      TVector3 unit_probe_dir = p4v.Vect().Unit();
      TVector3 p3l( pxl, pyl, pzl );
      p3l.RotateUz( unit_probe_dir );

      p4l = TLorentzVector( p3l, El );

      // Compute the 4-momentum transfer and Q^2
      TLorentzVector q4 = p4v - p4l;
      Q2 = -1. * q4.Mag2();

      // If we're below the minimum Q^2 value, then just go back to rejection
      // sampling
      if ( Q2 < Q2min ) continue;

      // Reset the 'trust' bits before exiting the loop
      interaction->ResetBit( kISkipProcessChk );
      interaction->ResetBit( kISkipKinematicChk );

      LOG( "SuSAv2Inel", pINFO ) << "Selected: W = " << W
        << ", cth = " << cth << ", Tl = " << Tl;

      break;
    }
  } // iterations of the rejection sampling loop

  // Store the differential cross section for the selected kinematics
  evrec->SetDiffXSec( xsec, kPSTlctl );

  // Compute x and y from W and Q^2
  // NOTE: The values are currently approximate because we assume an on-shell
  // initial nucleon at rest (the hadronic system has not been generated yet
  // since the SuSAv2 differential cross section is inclusive).
  // TODO: revisit and do a better job later
  double x = -1.;
  double y = -1.;
  double M = genie::constants::kNucleonMass;
  genie::utils::kinematics::WQ2toXY( Ev, M, W, Q2, x, y );

  // Store selected kinematics and clear running values
  interaction->KinePtr()->SetQ2( Q2, true );
  interaction->KinePtr()->SetW( W,  true );
  interaction->KinePtr()->Setx( x,  true );
  interaction->KinePtr()->Sety( y,  true );
  interaction->KinePtr()->SetFSLeptonP4( p4l );
  interaction->KinePtr()->ClearRunningValues();

  // Create a particle object for the outgoing lepton. Make it the daughter of
  // the probe and assign it the same 4-position in the event record.
  int pdgc = interaction->FSPrimLepton()->PdgCode();
  int mom_index = evrec->ProbePosition();
  TLorentzVector pos4v( *evrec->Probe()->X4() );
  evrec->AddParticle( pdgc, kIStStableFinalState, mom_index, -1, -1, -1,
    p4l, pos4v );
}

//___________________________________________________________________________
void genie::SuSAv2InelGenerator::Configure( const genie::Registry& config )
{
  genie::Algorithm::Configure( config );
  this->LoadConfig();
}

//____________________________________________________________________________
void genie::SuSAv2InelGenerator::Configure( std::string config )
{
  genie::Algorithm::Configure( config );
  this->LoadConfig();
}

//____________________________________________________________________________
void genie::SuSAv2InelGenerator::LoadConfig()
{
  // Safety factor for the maximum differential cross section
  this->GetParamDef( "MaxXSec-SafetyFactor", fSafetyFactor, 1.25 );

  // Minimum energy for which max xsec would be cached, forcing explicit
  // calculation for lower eneries
  this->GetParamDef( "Cache-MinEnergy", fEMin, 0.5 );

  // Load Wcut used in DIS/RES join scheme
  this->GetParam( "Wcut", fWcut );

  // Maximum allowed fractional cross section deviation from maxim cross
  // section used in rejection method
  this->GetParamDef( "MaxXSec-DiffTolerance", fMaxXSecDiffTolerance, 999999. );
  assert( fMaxXSecDiffTolerance >= 0. );
}

//____________________________________________________________________________
double genie::SuSAv2InelGenerator::ComputeMaxXSec(
  const genie::Interaction* interaction ) const
{
  // Computes the maximum differential cross section in the requested phase
  // space. This method overloads the KineGeneratorWithCache::ComputeMaxXSec()
  // method. The value is cached at a circular cache branch for retrieval
  // during subsequent event generation. The computed max differential cross
  // section does not need to be the exact maximum. The number used in the
  // rejection method will be scaled up by a safety factor. But this needs to
  // be fast -- do not use a very fine grid.

  double max_xsec = 0.;

  double Ev = interaction->InitState().ProbeE( genie::kRfLab );
  double ml = interaction->FSPrimLepton()->Mass();

  // Use a constant cos(theta) value for now
  const double cth = 0.9;
  interaction->KinePtr()->SetKV( genie::kKVctl, cth );

  // 2D scan in W, Tl
  int num_W = 20;
  double W_min = genie::constants::kNeutronMass + genie::constants::kPionMass;
  double W_max = std::max( Ev - W_min, 0. );

  if ( W_max < W_min ) return 0.;

  double dW = ( W_max - W_min ) / ( num_W - 1 );

  int num_Tl = 25;
  double Tl_min = 0.;
  double Tl_max = Ev - ml;

  double dTl = ( Tl_max - Tl_min ) / ( num_Tl - 1 );

  for ( int iw = 0; iw < num_W; ++iw ) {

    double W = W_min + iw*dW;
    interaction->KinePtr()->SetW( W );

    for ( int iTl = 0; iTl < num_Tl; ++iTl ) {

      double Tl = Tl_min + iTl*dTl;
      interaction->KinePtr()->SetKV( genie::kKVTl, Tl );

      double xsec = fXSecModel->XSec( interaction, kPSTlctl );
      LOG( "SuSAv2Inel", pDEBUG ) << "xsec(W= " << W << ", Tl= " << Tl
        << ", cth = " << cth << ") = " << xsec;
      max_xsec = std::max( xsec, max_xsec );
    } // Tl
  } // W

  // Apply safety factor
  max_xsec *= fSafetyFactor;

  return max_xsec;
}
