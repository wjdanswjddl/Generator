// Standard library includes
#include <iostream>
#include <memory>

// ROOT includes
#include "TGenPhaseSpace.h"
#include "TRandom3.h"

// GENIE includes
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/Units.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/RunOpt.h"
#include "Physics/NuclearState/NuclearModelI.h"
#include "Physics/NuclearState/SpectralFunc.h"
#include "Physics/QuasiElastic/XSection/QELUtils.h"

#include "FortranWrapperXSecIntegrator.h"

// The number of Interaction objects to create in the event loop below
constexpr int NUM_EVENTS = 100;
// Upper limit on the number of phase space throws to use (avoids a potentially
// infinite loop)
constexpr int MAX_PS_ITERATIONS = 10000;

int main() {

  // Set up a ROOT random number generator
  TRandom3 myRandom;

  // Set up the default GENIE v3 model set
  auto* ro = genie::RunOpt::Instance();
  ro->SetTuneName( "G18_02a_00_000" );
  ro->BuildTune();

  // Configure the input interaction
  int probe_pdg = 11; // electron
  double Eprobe = 0.961; // GeV
  int target_pdg = 1000060120; // 12C
  int struck_nucleon_pdg = genie::kPdgProton; // genie::kPdgNeutron
  double probe_mass = genie::PDGLibrary::Instance()
    ->Find( probe_pdg )->Mass(); // GeV
  double p_probe = std::sqrt( Eprobe*Eprobe - probe_mass*probe_mass );

  // The probe is initially directed along +z (for simplicity)
  TLorentzVector p4Probe = TLorentzVector( 0., 0., p_probe, Eprobe );

  genie::Interaction* interaction = genie::Interaction::QELEM(
    target_pdg, struck_nucleon_pdg, probe_pdg, p4Probe );

  // Set the masses of the two-body final state for QE scattering. The
  // recoiling nucleus is not explicitly treated
  double ml = interaction->FSPrimLepton()->Mass();
  double mNf = interaction->RecoilNucleon()->Mass();
  constexpr unsigned NUM_FINAL = 2;
  double masses[NUM_FINAL] = { ml, mNf };

  genie::AlgFactory* factory = genie::AlgFactory::Instance();

  const auto* xsec_model = dynamic_cast< const genie::XSecAlgorithmI* >(
    factory->GetAlgorithm("genie::FortranWrapperQELPXSec", "Default") );

  const auto* nucl_model = dynamic_cast< const genie::NuclearModelI* >(
    factory->GetAlgorithm("genie::SpectralFunc", "Default") );

  // Compute the total cross section
  const genie::FortranWrapperXSecIntegrator* integrator
    = dynamic_cast< const genie::FortranWrapperXSecIntegrator* >(
    xsec_model->SubAlg( "XSec-Integrator" )
  );
  assert( integrator );

  double integ = integrator->Integrate( xsec_model, interaction );
  double max = integrator->GetMaxDiffXSec();
  std::cout << "TOTAL XSEC = " << integ << '\n';
  std::cout << "MAX XSEC = " << max << '\n';

  for ( int e = 0; e < NUM_EVENTS; ++e ) {

    // Draw a nucleon from the nuclear model and update its 4-momentum in the
    // interaction based on its 3-momentum and removal energy
    nucl_model->GenerateNucleon( interaction->InitState().Tgt() );
    double dummy = 0.;
    genie::utils::BindHitNucleon( *interaction, *nucl_model, dummy,
      genie::kUseNuclearModel );

    // Determine the total 4-momentum of the two-body system (probe, bound
    // nucleon) in the initial state
    TLorentzVector p4Ni = interaction->InitState().Tgt().HitNucP4();
    TLorentzVector p4tot = p4Probe + p4Ni;

    // Check that we're above the threshold for creating the final-state
    // particles. If we're not, then start over with a new nucleon
    double sqrt_s = interaction->InitState().CMEnergy();
    double diff = sqrt_s - ml - mNf;
    if ( diff < 0. ) {
      std::cerr << "This event is below threshold. I will try again.\n";
      continue;
    }

    // Set up the phase-space generator to randomly sample points in the
    // final-state phase space
    TGenPhaseSpace event;
    event.SetDecay( p4tot, NUM_FINAL, masses);

    // Maximum weight for sampling over the allowed phase space
    double max_weight_ps = event.GetWtMax();

    // Choose an unweighted phase space point uniformly using the generator
    bool accepted = false;
    int ps_counter = 0;
    do {

      // Generate a weighted event from within the allowed phase space
      Double_t weight_ps = event.Generate();

      // Use rejection sampling to decide whether to keep this event (thus
      // allowing us to treat it as unweighted if it is accepted)
      double y = myRandom.Uniform( 0., max_weight_ps );
      accepted = ( y <= weight_ps );

      ++ps_counter;

    } while ( !accepted && ps_counter < MAX_PS_ITERATIONS );

    // If the accepted flag is false here, we got stuck and iterated too long
    // trying to find a suitable phase space point
    if ( !accepted ) {
      std::cerr << "I got stuck sampling phase space points. Time to try"
        << " again with a new initial state.\n";
      continue;
    }

    // Retrieve the final particle 4-momenta from the phase-space generator
    TLorentzVector* p4Lep = event.GetDecay( 0 );
    TLorentzVector* p4Nf  = event.GetDecay( 1 );

    // Update the interaction object with the sampled values
    interaction->KinePtr()->SetFSLeptonP4( *p4Lep );
    interaction->KinePtr()->SetHadSystP4( *p4Nf );

    // Compute the differential cross section
    std::cout << "EVENT " << e << '\n';
    double xsec = xsec_model->XSec( interaction, genie::kPSFullNBody );
    std::cout << "XSEC = " << xsec << '\n';

  } // event loop

  std::cout << "DONE\n";
  return 0;
}
