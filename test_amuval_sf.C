// Run with the genie ROOT interpreter (instead of plain root) to load the
// needed class dictionaries

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Interaction/Kinematics.h"
#include "Framework/Utils/RunOpt.h"
#include "Physics/DeepInelastic/XSection/AMUValStrucFunc.h"

const int probe = 14; // numu
const int target = 1000260560; // 56Fe
const int hitnuc = 2112; // neutron

void test_amuval_sf() {

  // Set up the dummy tune to use for testing the AMU-Valencia structure
  // function calculation
  genie::RunOpt* ro = genie::RunOpt::Instance();
  ro->SetTuneName( "H20_02a_00_000" );
  ro->BuildTune();

  // Retrieve the algorithm that calculates the structure functions
  const genie::AMUValStrucFunc* sf = dynamic_cast<const genie::AMUValStrucFunc*>(
    genie::AlgFactory::Instance()->GetAlgorithm("genie::AMUValStrucFunc",
    "Default") );

  // Create a dummy DIS-CC Interaction object
  genie::Interaction* interaction = genie::Interaction::DISCC( target,
    hitnuc, probe );

  // Set the value of Q^2 and Bjorken x in the interaction
  double Q2 = 2.; // GeV^2
  double x = 0.5;

  interaction->KinePtr()->SetQ2( Q2 );
  interaction->KinePtr()->Setx( x );

  // Compute the structure functions
  sf->Calculate( interaction );

  // Output the result
  std::cout << "Q^2 = " << Q2 << " GeV\n";
  std::cout << "x = " << x << '\n';
  std::cout << "F1 = " << sf->F1() << '\n';
  std::cout << "F2 = " << sf->F2() << '\n';
  std::cout << "F3 = " << sf->F3() << '\n';
}
