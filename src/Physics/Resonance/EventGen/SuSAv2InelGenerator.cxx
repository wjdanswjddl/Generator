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
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
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

    cth = cth_min + ( cth_max - cth_min ) * rnd->RndKine().Rndm();
    Tl = Tl_min + ( Tl_max - Tl_min ) * rnd->RndKine().Rndm();

    double W_max = genie::constants::kNeutronMass + Ev - Tl - ml;
    if ( W_max < W_min ) continue;

    W = W_min + ( W_max - W_min ) * rnd->RndKine().Rndm();

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

    LOG( "SuSAv2Inel", pDEBUG ) << "xsec_max = " << xsec_max << ", y = " << y;

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

  double Ev = interaction->InitState().ProbeE( genie::kRfLab );
  double ml = interaction->FSPrimLepton()->Mass();

  // Start with a guess of where the maximum xsec lies in phase space
  genie::PDGLibrary* pdg_library = genie::PDGLibrary::Instance();
  double mDelta = pdg_library->Find( genie::kPdgP33m1232_DeltaP )->Mass();

  double W = mDelta;
  double Tl = 0.9516*Ev - 0.3228;
  double cth = std::min( 1., 0.5395*std::exp(-1.984*Ev) + 0.9956 );

  interaction->KinePtr()->SetW( W );
  interaction->KinePtr()->SetKV( genie::kKVTl, Tl );
  interaction->KinePtr()->SetKV( genie::kKVctl, cth );

  double max_xsec = fXSecModel->XSec( interaction, kPSTlctl );
  LOG( "SuSAv2Inel", pDEBUG ) << "xsec(W= " << W << ", Tl= " << Tl
    << ", cth = " << cth << ") = " << max_xsec;

  // Step around the guess in all directions with a decreasing step size,
  // move on each iteration towards the maximum value
  double step_W = 0.5;
  double step_Tl = 0.5;
  double step_cth = 0.5;
  const int MAX_STEPS = 1000;
  const int MIN_STEP_W = 0.01;
  int step_iter = 0;

  while ( step_iter < MAX_STEPS && step_W > MIN_STEP_W ) {
    std::vector< double > Ws = { W - step_W, W, W + step_W };
    std::vector< double > Tls = { Tl - step_Tl, Tl, Tl + step_Tl };
    std::vector< double > cths = { cth - step_cth, cth, cth + step_cth };

    std::vector< genie::Interaction > interaction_vec;
    for ( const double& cur_W : Ws ) {
      for ( const double& cur_Tl : Tls ) {
        for ( const double& cur_cth : cths ) {
          interaction_vec.emplace_back( *interaction );
          auto& last_interaction = interaction_vec.back();
          last_interaction.KinePtr()->SetW( cur_W );
          last_interaction.KinePtr()->SetKV( genie::kKVTl, cur_Tl );
          last_interaction.KinePtr()->SetKV( genie::kKVctl, cur_cth );
        }
      }
    }

    const genie::Interaction* max_inter = interaction;
    for ( const auto& test_inter : interaction_vec ) {
      double xsec = fXSecModel->XSec( &test_inter, kPSTlctl );
      if ( xsec > max_xsec ) {
        LOG( "SuSAv2Inel", pDEBUG ) << "xsec(W= " << W << ", Tl= " << Tl
          << ", cth = " << cth << ") = " << xsec;
        max_inter = &test_inter;
        max_xsec = xsec;
        W = test_inter.Kine().W();
        Tl = test_inter.Kine().GetKV( genie::kKVTl );
        cth = test_inter.Kine().GetKV( genie::kKVctl );
      }
    }

    // If we found a phase-space point with a higher differential cross
    // section, then just continue the loop. If the current choice is still
    // best, decrease the step sizes for the next attempt.
    if ( max_inter == interaction ) {
      step_W /= 5.;
      step_Tl /= 5.;
      step_cth /= 5.;
    }
  }

  LOG( "SuSAv2Inel", pDEBUG ) << "max_xsec(W= " << W << ", Tl= " << Tl
    << ", cth = " << cth << ") = " << max_xsec;

  // Apply safety factor
  max_xsec *= fSafetyFactor;

  return max_xsec;
}
