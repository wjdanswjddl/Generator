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
