//____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Steven Gardiner <gardiner \at fnal.gov>
 Fermi National Accelerator Laboratory
*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/EventGen/InteractionList.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/Resonance/EventGen/SuSAv2InelInteractionListGenerator.h"

using namespace genie;

//___________________________________________________________________________
SuSAv2InelInteractionListGenerator::SuSAv2InelInteractionListGenerator() :
  InteractionListGeneratorI( "genie::SuSAv2InelInteractionListGenerator" )
{

}

//___________________________________________________________________________
SuSAv2InelInteractionListGenerator::SuSAv2InelInteractionListGenerator(
  std::string config ) : InteractionListGeneratorI(
  "genie::SuSAv2InelInteractionListGenerator", config )
{

}

//___________________________________________________________________________
SuSAv2InelInteractionListGenerator::~SuSAv2InelInteractionListGenerator()
{

}

//___________________________________________________________________________
InteractionList* SuSAv2InelInteractionListGenerator::CreateInteractionList(
  const InitialState& init_state ) const
{
  LOG( "IntLst", pINFO ) << "InitialState = " << init_state.AsString();

  // Specify the requested interaction type
  InteractionType_t inttype;
  if      (fIsCC) inttype = kIntWeakCC;
  else if (fIsNC) inttype = kIntWeakNC;
  else if (fIsEM) inttype = kIntEM;
  else {
    LOG( "IntLst", pWARN )
      << "Unknown InteractionType! Returning NULL InteractionList "
      << "for init-state: " << init_state.AsString();
    return 0;
  }

  // Create a process information object
  ProcessInfo proc_info( kScResonant, inttype );

  // Create an interaction list
  InteractionList* int_list = new InteractionList;

  // Learn whether the input nuclear or free target has protons and neutrons
  // available or not
  const Target& inp_target = init_state.Tgt();
  bool hasP = ( inp_target.Z() > 0 );
  bool hasN = ( inp_target.N() > 0 );

  // Possible hit nucleons
  const int hit_nucleon[2] = { kPdgProton, kPdgNeutron };

  // Loop over hit nucleons
  for ( int i = 0; i < 2; ++i ) {

    // Proceed only if the hit nucleon exists in the current initial state
    if ( hit_nucleon[i] == kPdgProton  && !hasP ) continue;
    if ( hit_nucleon[i] == kPdgNeutron && !hasN ) continue;

    // Create an interaction
    Interaction* interaction = new Interaction( init_state, proc_info );

    // Set the struck nucleon PDG code
    Target* target = interaction->InitStatePtr()->TgtPtr();
    target->SetHitNucPdg( hit_nucleon[i] );

    // Add the interaction to the interaction list
    int_list->push_back( interaction );
  }

  if ( int_list->size() == 0u ) {
    LOG( "IntLst", pERROR ) << "Returning NULL InteractionList for"
      " init-state: " << init_state.AsString();
    delete int_list;
    return 0;
  }

  return int_list;
}

//___________________________________________________________________________
void SuSAv2InelInteractionListGenerator::Configure( const Registry& config )
{
  Algorithm::Configure( config );
  this->LoadConfigData();
}

//____________________________________________________________________________
void SuSAv2InelInteractionListGenerator::Configure( std::string config )
{
  Algorithm::Configure( config );
  this->LoadConfigData();
}

//____________________________________________________________________________
void SuSAv2InelInteractionListGenerator::LoadConfigData()
{
  this->GetParamDef( "is-CC", fIsCC, false );
  this->GetParamDef( "is-NC", fIsNC, false );
  this->GetParamDef( "is-EM", fIsEM, false );
}
