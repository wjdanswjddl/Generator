#ifndef _INCL_test_H_
#define _INCL_test_H_
#include <string>
#include <TGenPhaseSpace.h>
#include "GHEP/GHepRecord.h"
#include "EVGCore/EventRecord.h"
#include "Algorithm/AlgFactory.h"
#include "EVGCore/EventRecordVisitorI.h"
#include "Conventions/GMode.h"
#include "HadronTransport/INukeMode.h"
#include "HadronTransport/INukeHadroFates.h"

class TLorentzVector;
class TVector3;
//string kDefOptEvFilePrefix="gntp.inuke";
namespace genie{

 class GHepParticle;
 class INukeHadroData;
 class PDGCodeList;
 class HINCLCascade;
 class INCLCascade: public EventRecordVisitorI{


  public :
  INCLCascade();
  INCLCascade(string name);
  INCLCascade(string name, string config);
  ~INCLCascade();
  int INCLcascade(int arg, char * test[],GHepRecord * event_rec) const;
//  int INCLtopdgcode(int A, int Z)const;
  int pdgcpiontoA(int pdgc)const;
  int pdgcpiontoZ(int pdgc)const;
    //void INCLprint() const;
    // implement the EventRecordVisitorI interface
  virtual void ProcessEventRecord(GHepRecord * event_rec) const;

    //virtual void ProcessEventRecordI(int arg, char *test[]) const;

    // override the Algorithm::Configure methods to load configuration
    // data to protected data members
protected:
 bool CanRescatter(const GHepParticle * p) const;
 bool IsInNucleus(const GHepParticle * p) const;
 void TransportHadrons(GHepRecord * evrec) const;
 bool NeedsRescattering(const GHepParticle * p) const;
   mutable int            fRemnA;         ///< remnant nucleus A
   mutable int            fRemnZ;         ///< remnant nucleus Z
   mutable double         fTrackingRadius;
   mutable TLorentzVector fRemnP4;        ///< P4 of remnant system
   mutable GEvGenMode_t   fGMode;
   double       fR0;           ///< effective nuclear size param
   double       fNR;
 };
}
#endif
