//____________________________________________________________________________
/*!

\class    genie::INCLCascade

\brief    Link to INCL++ intranuclear hadron transport MC.
          Is a concrete implementation of the EventRecordVisitorI interface.

          Current INTRANUKE development is led by Pittsburgh/Antananarivo group

\author   Steve Dytman <dytman+@pitt.edu>, Pittsburgh University
          Marc Vololoniaina <narymarc@yahoo.com>, Pittsburgh University/Antananarivo

\created  May 23, 2017

\cpright  Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

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
  int INCLcascade(int arg, char * test[]) const;
  int INCLtopdgcode(int A, int Z)const;
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
