//____________________________________________________________________________
/*!

\class    genie::FortranPXSecI

\brief    Interface that allows GENIE to compute a differential cross section
          using a model implemented in Fortran.
          Is an extension of the XSecAlgorithmI interface. \n

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  8 October 2018

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _FORTRAN_CROSS_SECTION_H_
#define _FORTRAN_CROSS_SECTION_H_

// standard library includes
#include <string>

// GENIE includes
#include "Framework/EventGen/XSecAlgorithmI.h"

// extern "C" wrappers for Fortran functions go here
extern "C" {

  // All function arguments that are not marked with the "value" attribute
  // should be passed as pointers. Fortran defaults to pass-by-reference
  // for all function arguments.
  void compute_my_diff_xsec(const genie::Interaction* interaction,
    genie::KinePhaseSpace_t kps, double* xsec);

}

namespace genie {

class XSecIntegratorI;

class FortranPXSecI : public XSecAlgorithmI {

public:

  FortranPXSecI();
  FortranPXSecI(const std::string& config);
  virtual ~FortranPXSecI();

  // XSecAlgorithmI interface implementation
  double XSec            (const Interaction* i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction* i) const;
  bool   ValidProcess    (const Interaction* i) const;

  // Override the Algorithm::Configure methods to load configuration
  // data to private data members
  void Configure (const Registry& config);
  void Configure (std::string param_set);

private:

  void LoadConfig (void);

  /// Algorithm to use to integrate this cross section
  const XSecIntegratorI* fXSecIntegrator;
};

} // genie namespace

#endif
