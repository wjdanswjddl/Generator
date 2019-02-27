//____________________________________________________________________________
/*!

\class    genie::SpectralFunc

\brief    A realistic spectral function - based nuclear model.
          Is a concrete implementation of the NuclearModelI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 07, 2004

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE

*/
//____________________________________________________________________________

#ifndef _SPECTRAL_FUNCTION_H_
#define _SPECTRAL_FUNCTION_H_

#include "Physics/NuclearState/NuclearModelI.h"

class TNtupleD;
class TGraph2D;

namespace genie {

class SpectralFunc : public NuclearModelI {

public:
  SpectralFunc();
  SpectralFunc(string config);
  virtual ~SpectralFunc();

  using NuclearModelI::GenerateNucleon;  // inherit versions not overridden here
  using NuclearModelI::Prob;

  //-- implement the NuclearModelI interface
  bool           GenerateNucleon (const Target & t) const;

  // TODO: Right now, SpectralFunc::Prob() appears to return the
  // probability *density*, not the probability. Think about changing the
  // behavior if doing so makes sense.
  double         Prob            (double p, double w, const Target & t) const;
  NuclearModel_t ModelType       (const Target &) const
  {
    return kNucmSpectralFunc;
  }

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string config);

  // Add "override" specifiers when GENIE moves to C++11
  // Note that, for this nuclear model, the removal energy is constant
  virtual double MinRemovalEnergy(const Target& t, double) const /*override*/;
  virtual double MaxRemovalEnergy(const Target& t, double) const /*override*/;

  virtual double MinMomentum(const Target&, double) const /*override*/;
  virtual double MaxMomentum(const Target&, double) const /*override*/;

  virtual double ProbDensity(double p, double w, const Target & t, double r = 0.)
    const /*override*/;

private:
  void       LoadConfig             (void);
  TGraph2D * Convert2Graph          (TNtupleD & data) const;
  TGraph2D * SelectSpectralFunction (const Target & target) const;

  TGraph2D * fSfFe56;   ///< Benhar's Fe56 SF
  TGraph2D * fSfC12;    ///< Benhar's C12 SF
};

}      // genie namespace
#endif // _SPECTRAL_FUNCTION_H_
