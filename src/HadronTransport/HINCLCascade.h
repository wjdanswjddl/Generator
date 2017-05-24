//____________________________________________________________________________
/*!

\class    genie::HINCLCascade

\brief    Calls INCLCascade which links to INCL++.

          Current INTRANUKE development is led by Pittsburgh/Antananarivo group

\author   Steve Dytman <dytman+@pitt.edu>, Pittsburgh University
          Marc Vololoniaina <narymarc@yahoo.com>, Pittsburgh University/Antananarivo

\created  May 23, 2017

\cpright  Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _HINCL_test_H_
#define _HINCL_test_H_


#include <TGenPhaseSpace.h>

#include "Algorithm/AlgFactory.h"
#include "EVGCore/EventRecordVisitorI.h"
#include "HadronTransport/INukeMode.h"
#include "HadronTransport/INukeHadroFates.h"
#include "Interfaces/NuclearModelI.h"
#include "INCLCascade.h"

class TLorentzVector;
class TVector3;

namespace genie{
  class GHepParticle;
  class INukeHadroData;
  class PDGCodeList;
class HINCLCascade: public INCLCascade{
  friend class IntranukeTester;
public :
 HINCLCascade();
 HINCLCascade(string config);

 ~HINCLCascade();
//void HINCLCascadePr(void);
 void ProcessEventRecord(GHepRecord * event_rec) const;
};
}
#endif
