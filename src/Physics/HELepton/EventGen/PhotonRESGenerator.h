//____________________________________________________________________________
/*!

\class    genie::PhotonRESGenerator

\brief    Generator for trident production.

\author   Alfonso Garcia <aagarciasoto \at km3net.de>
          IFIC & Harvard University

\created  Dec 8, 2021

\cpright  Copyright (c) 2003-2024, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _PHOTON_RES_GENERATOR_H_
#define _PHOTON_RES_GENERATOR_H_

#include "Framework/EventGen/EventRecordVisitorI.h"

using namespace genie;

namespace genie {

class PhotonRESGenerator : public EventRecordVisitorI {

public :
  PhotonRESGenerator();
  PhotonRESGenerator(string config);
 ~PhotonRESGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord (GHepRecord * event) const;

  // Overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  virtual void Configure(const Registry & config);
  virtual void Configure(string config);

private:

  void LoadConfig         (void);

  const EventRecordVisitorI * fWDecayer; ///< PYTHIA W decayer

};

}      // genie namespace
#endif // _PHOTON_RES_GENERATOR_H_
