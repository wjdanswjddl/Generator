//____________________________________________________________________________
/*!

\class    genie::QELEventGeneratorGCF

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Acclerator Laboratory

\created  Sep 11, 2020

\cpright  Copyright (c) 2003-2020, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________


#ifndef _QEL_EVENT_GENERATORGCF_H_
#define _QEL_EVENT_GENERATORGCF_H_

#include "TChain.h"

#include "Physics/Common/KineGeneratorWithCache.h"
#include "Framework/Conventions/Controls.h"


namespace genie {

class QELEventGeneratorGCF: public KineGeneratorWithCache {

public :
  QELEventGeneratorGCF();
  QELEventGeneratorGCF(string config);
 ~QELEventGeneratorGCF();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void   LoadConfig     (void);
  double ComputeMaxXSec(const Interaction* in) const;
  void AddTargetNucleusRemnant(GHepRecord* evrec) const;

  mutable TChain* fTree;
  mutable int fTreeEntry;

}; // class definition

} // genie namespace

#endif // _QEL_EVENT_GENERATORGCF_H_
