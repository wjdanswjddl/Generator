//____________________________________________________________________________
/*!

\class    genie::SpectralFunc1d

\brief    Simpler approach to using spectral functions.
	  A beta version.
          Implements the NuclearModelI interface.

\ref

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  October 09, 2004

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _SPECTRAL_FUNCTION_1D_H_
#define _SPECTRAL_FUNCTION_1D_H_

#include <map>

#include "Physics/NuclearState/NuclearModelI.h"

using std::map;

namespace genie {

class Spline;
class SpectralFunc1d : public NuclearModelI {

public:
  SpectralFunc1d();
  SpectralFunc1d(string config);
  virtual ~SpectralFunc1d();

  using NuclearModelI::GenerateNucleon;  // inherit versions not overridden here
  using NuclearModelI::Prob;

  //-- implement the NuclearModelI interface
  bool           GenerateNucleon (const Target & t) const;
  double         Prob            (double p, double w, const Target & t) const;
  NuclearModel_t ModelType       (const Target &) const
  {
    return kNucmFermiGas; /// is not really a spectral func model, just a FG model with different momentum distribution
  }

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

  // TODO: IMPLEMENT THESE!
  // Add "override" specifiers when GENIE moves to C++11
  // Note that, for this nuclear model, the removal energy is constant
  inline virtual double MinRemovalEnergy(const Target&, double) const /*override*/
    { return 0.; }
  inline virtual double MaxRemovalEnergy(const Target&, double) const /*override*/
    { return 0.; }

  inline virtual double MinMomentum(const Target&, double) const /*override*/
    { return 0.; }
  inline virtual double MaxMomentum(const Target&, double) const /*override*/
    { return 0.; }

  inline virtual double ProbDensity(double /*p*/, double /*w*/, const Target& /*t*/,
    double /*r*/) const /*override*/ { return 0.; }

private:
  void LoadConfig (void);
  void CleanUp    (void);

  // Spectral function data
  // Hopefully, analytical expressions for spectral functions will become available soon.
  //
  bool   fUseRFGRemovalE;
  bool   fUseRFGMomentumCutoff ;
  double fPCutOff;
  map<int, Spline *> fSFk;     ///< All available spectral funcs integrated over removal energy
  map<int, Spline *> fSFw;     ///< Average nucleon removal as a function of pF - computed from the spectral function
  map<int, double>   fNucRmvE; ///< Removal energies as used in FG model
  map<int, double>   fMaxProb; ///< Max SF(k) probability used in rejection method
};

}         // genie namespace
#endif    // _SPECTRAL_FUNCTION_1D_H_

