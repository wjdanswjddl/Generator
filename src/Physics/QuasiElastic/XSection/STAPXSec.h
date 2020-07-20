//____________________________________________________________________________
/*!

\class    genie::STAPXSec

\brief    Computes the inclusive differential cross section for lepton-nucleus
          scattering according to the short-time approximation (STA) approach.
          Uses precomputed hadron tensor tables.
          Is a concrete implementation of the XSecAlgorithmI interface.

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Acclerator Laboratory

\ref      \TODO: add reference(s)

\created  July 20, 2020

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _STA_PXSEC_H_
#define _STA_PXSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/HadronTensors/HadronTensorI.h"
#include "Physics/HadronTensors/HadronTensorModelI.h"

namespace genie {

class XSecIntegratorI;

class STAPXSec : public XSecAlgorithmI {

public:

  STAPXSec();
  STAPXSec(string config);
  virtual ~STAPXSec();

  // XSecAlgorithmI interface implementation
  double XSec(const Interaction* i, KinePhaseSpace_t k) const;
  double Integral(const Interaction* i) const;
  bool   ValidProcess(const Interaction* i) const;

  // override the Algorithm::Configure methods to load configuration
  // data to private data members
  void Configure (const Registry & config);
  void Configure (string config);

private:

  /// Load algorithm configuration
  void LoadConfig (void);

  /// External scaling factor for this cross section
  double fXSecScale;

  const HadronTensorModelI* fHadronTensorModel;

  /// GSL numerical integrator
  const XSecIntegratorI*  fXSecIntegrator;

};

} // genie namespace
#endif // _STA_PXSEC_H_
