//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Steve Dytman <dytman+@pitt.edu>, Pittsburgh Univ.
         Marc Vololoniaina <narymarc@yahoo.com>, Univ. Antananarivo, Madagascar/Pittsburgh Univ.

         May 23, 2017

 For the class documentation see the corresponding header file.

 Important revisions after version 2.12.8 :
*/
//____________________________________________________________________________
//---------------------
#include <cstdlib>
#include <sstream>

#include <iostream>
#include <iomanip>
#include <string>

#include <cassert>
#include <cstdlib>
#include <map>
#include <TMath.h>

#include "G4INCLCascade.hh"
#include "G4INCLConfigEnums.hh"
#include "DatafilePaths.hh" 
#include "G4INCLParticle.hh"

#include "Algorithm/AlgConfigPool.h"
#include "Algorithm/AlgFactory.h"
#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/Controls.h"
#include "EVGCore/EVGThreadException.h"
#include "GHEP/GHepFlags.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepParticle.h"
#include "HadronTransport/Intranuke.h"
#include "HadronTransport/HNIntranuke.h"
#include "HadronTransport/INukeException.h"
#include "HadronTransport/INukeHadroData.h"
#include "HadronTransport/INukeUtils.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Numerical/Spline.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGCodeList.h"
#include "PDG/PDGUtils.h"
#include "Utils/PrintUtils.h"
#include "Utils/NuclearUtils.h"
#include "G4INCLCascade.hh"
#include "G4INCLIAvatar.hh"
#include "G4INCLIPropagationModel.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLStandardPropagationModel.hh"
#include "G4INCLRandom.hh"
#include "G4INCLRanecu.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLKinematicsUtils.hh"

#include "HINCLCascade.h"
#include "INCLCascade.h"

// signal handler (for Linux and GCC)

#include "G4INCLConfig.hh"
#include "G4INCLVersion.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLUnorderedVector.hh"
#include "G4INCLParticle.hh"
#include "G4INCLEventInfo.hh"
#include "G4INCLStore.hh"
#include "G4INCLIAvatar.hh"


using namespace genie;
using namespace genie::utils;
using namespace genie::utils::intranuke;
using namespace genie::constants;
using namespace genie::controls;
using namespace G4INCL;
using std::ostringstream;

HINCLCascade::HINCLCascade() :
INCLCascade("genie::HINCLCascade")
{

}
//___________________________________________________________________________
HINCLCascade::HINCLCascade(string config) :
INCLCascade("genie::HINCLCascade",config)
{

}
//___________________________________________________________________________
HINCLCascade::~HINCLCascade()
{

}
void HINCLCascade::ProcessEventRecord(GHepRecord * evrec) const{
  LOG("HINCLCascade", pNOTICE)
     << "************ Running HINCLCascade MODE INTRANUKE ************";

 INCLCascade::ProcessEventRecord(evrec);
//INCLCascade::INCLcascade();
  LOG("HINCLCascade", pINFO) << "Done with this event";
  //this->TransportHadrons(evrec);
}


