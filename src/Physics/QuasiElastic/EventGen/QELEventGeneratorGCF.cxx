//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author:  Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

 For the class documentation see the corresponding header file.
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen//RunningThreadInfo.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Physics/QuasiElastic/EventGen/QELEventGeneratorGCF.h"


#include "Framework/Utils/Range1.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/Cache.h"
#include "Framework/Utils/CacheBranchFx.h"

using namespace genie;
using namespace genie::controls;
using namespace genie::utils;

//___________________________________________________________________________
QELEventGeneratorGCF::QELEventGeneratorGCF() :
KineGeneratorWithCache("genie::QELEventGeneratorGCF")
{
  fTree = NULL;
}
//___________________________________________________________________________
QELEventGeneratorGCF::QELEventGeneratorGCF(string config) :
KineGeneratorWithCache("genie::QELEventGeneratorGCF", config)
{
  fTree = NULL;
}
//___________________________________________________________________________
QELEventGeneratorGCF::~QELEventGeneratorGCF()
{
  if ( fTree ) delete fTree;
}
//___________________________________________________________________________
void QELEventGeneratorGCF::ProcessEventRecord(GHepRecord * evrec) const
{
  // Get the interaction and set the 'trust' bits
  Interaction * interaction = evrec->Summary();
  const InitialState & init_state = interaction -> InitState();
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  // Fetch the next event's information from the input ROOT file.
  // Use it to override information in the input event record.
  double Ev, weight;
  int recoil_nucleon_pdg;
  double p3mu[3], p3LeadN[3], p3RecN[3], p3Aminus2[3];
  fTree->SetBranchAddress( "Eneutrino", &Ev );
  fTree->SetBranchAddress( "rec_type", &recoil_nucleon_pdg );
  fTree->SetBranchAddress( "pk", &p3mu );
  fTree->SetBranchAddress( "pLead", &p3LeadN );
  fTree->SetBranchAddress( "pRec", &p3RecN );
  fTree->SetBranchAddress( "pAm2", &p3Aminus2 );
  fTree->SetBranchAddress( "weight", &weight );

  fTree->GetEntry( fTreeEntry );
  if ( fTreeEntry < fTree->GetEntries() - 1 ) ++fTreeEntry;


  // Set the event weight based on the input TTree contents
  evrec->SetWeight( weight );

  // Alter the probe, hit nucleon, and target in the GHepRecord and in the
  // interaction Summary. Ignore what was in there before.
  GHepParticle* probe = evrec->Probe();
  GHepParticle* hitnuc = evrec->HitNucleon();
  GHepParticle* target = evrec->TargetNucleus();

  // TODO: remove hard-coded stuff here
  probe->SetPdgCode( kPdgNuMu );
  probe->SetEnergy( Ev );
  probe->SetPx( 0. );
  probe->SetPy( 0. );
  probe->SetPz( Ev );

  const int ARGON_40 = 1000180400;
  target->SetPdgCode( ARGON_40 );
  interaction->InitStatePtr()->TgtPtr()->SetId( ARGON_40 );

  hitnuc->SetPdgCode( kPdgNeutron );

  interaction->InitStatePtr()->SetProbePdg( kPdgNuMu );
  // Also resets the 3-momentum components to match what we did above
  interaction->InitStatePtr()->SetProbeE( Ev );

  // Access the target from the interaction summary
  Target* tgt = init_state.TgtPtr();
  tgt->SetHitNucPdg( kPdgNeutron );

  double lepMass = interaction->FSPrimLepton()->Mass();
  double lepE = std::sqrt( p3mu[0]*p3mu[0] + p3mu[1]*p3mu[1] + p3mu[2]*p3mu[2]
    + lepMass*lepMass );

  TLorentzVector outLeptonMom( p3mu[0], p3mu[1], p3mu[2], lepE );

  double nuc1Mass = interaction->RecoilNucleon()->Mass();
  double nuc1E = std::sqrt( p3LeadN[0]*p3LeadN[0] + p3LeadN[1]*p3LeadN[1]
    + p3LeadN[2]*p3LeadN[2] + nuc1Mass*nuc1Mass );

  TLorentzVector outNucleonMom1( p3LeadN[0], p3LeadN[1], p3LeadN[2], nuc1E );

  double nuc2Mass = PDGLibrary::Instance()->Find( recoil_nucleon_pdg )->Mass();
  double nuc2E = std::sqrt( p3RecN[0]*p3RecN[0] + p3RecN[1]*p3RecN[1]
    + p3RecN[2]*p3RecN[2] + nuc2Mass*nuc2Mass );

  TLorentzVector outNucleonMom2( p3RecN[0], p3RecN[1], p3RecN[2], nuc2E );

  fTree->SetBranchAddress( "pAm2", &p3Aminus2 );

  TLorentzVector x4l( *(evrec->Probe())->X4() );

  evrec->AddParticle( interaction->FSPrimLeptonPdg(), kIStStableFinalState, evrec->ProbePosition(), -1, -1, -1, outLeptonMom, x4l);

  GHepStatus_t ist = tgt->IsNucleus() ? kIStHadronInTheNucleus : kIStStableFinalState;

  GHepParticle outNucleon1( interaction->RecoilNucleonPdg(), ist, evrec->HitNucleonPosition(), -1, -1, -1, outNucleonMom1, x4l);
  evrec->AddParticle( outNucleon1 );

  GHepParticle outNucleon2( recoil_nucleon_pdg, ist, evrec->HitNucleonPosition(), -1, -1, -1, outNucleonMom2, x4l);
  evrec->AddParticle( outNucleon2 );

  // Add a recoiled nucleus remnant. Skip if not a nuclear target.
  if ( tgt->IsNucleus() ) this->AddTargetNucleusRemnant(evrec);

  LOG("QELEvent", pINFO) << "Done generating QE event kinematics!";
}
//___________________________________________________________________________
void QELEventGeneratorGCF::AddTargetNucleusRemnant(GHepRecord * evrec) const
{
// add the remnant nuclear target at the GHEP record

  LOG("QELEvent", pINFO) << "Adding final state nucleus";

  double Px = 0;
  double Py = 0;
  double Pz = 0;
  double E  = 0;

  GHepParticle * nucleus = evrec->TargetNucleus();
  int A = nucleus->A();
  int Z = nucleus->Z();

  int fd = nucleus->FirstDaughter();
  int ld = nucleus->LastDaughter();

  for(int id = fd; id <= ld; id++) {

    // compute A,Z for final state nucleus & get its PDG code and its mass
    GHepParticle * particle = evrec->Particle(id);
    assert(particle);
    int  pdgc = particle->Pdg();
    bool is_p  = pdg::IsProton (pdgc);
    bool is_n  = pdg::IsNeutron(pdgc);

    if (is_p) Z--;
    if (is_p || is_n) A--;

    Px += particle->Px();
    Py += particle->Py();
    Pz += particle->Pz();
    E  += particle->E();

  }//daughters

  TParticlePDG * remn = 0;
  int ipdgc = pdg::IonPdgCode(A, Z);
  remn = PDGLibrary::Instance()->Find(ipdgc);
  if(!remn) {
    LOG("HadronicVtx", pFATAL)
          << "No particle with [A = " << A << ", Z = " << Z
                            << ", pdgc = " << ipdgc << "] in PDGLibrary!";
    assert(remn);
  }

  double Mi = nucleus->Mass();
  Px *= -1;
  Py *= -1;
  Pz *= -1;
  E = Mi-E;

  // Add the nucleus to the event record
  LOG("QELEvent", pINFO)
     << "Adding nucleus [A = " << A << ", Z = " << Z
     << ", pdgc = " << ipdgc << "]";

  int imom = evrec->TargetNucleusPosition();
  evrec->AddParticle(
       ipdgc,kIStStableFinalState, imom,-1,-1,-1, Px,Py,Pz,E, 0,0,0,0);

  LOG("QELEvent", pINFO) << "Done";
  LOG("QELEvent", pINFO) << *evrec;
}
//___________________________________________________________________________
void QELEventGeneratorGCF::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void QELEventGeneratorGCF::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void QELEventGeneratorGCF::LoadConfig(void)
{
  std::string events_file_name;
  GetParam( "EventsFile", events_file_name );

  std::string ttree_name;
  GetParamDef( "EventsTTree", ttree_name, std::string("genT") );

  if ( fTree ) delete fTree;

  fTree = new TChain( ttree_name.c_str() );
  fTree->Add( events_file_name.c_str() );

  fTreeEntry = 0;
}
//____________________________________________________________________________
double QELEventGeneratorGCF::ComputeMaxXSec(const Interaction*
  /*interaction*/) const
{
  return 1.;
}
