//____________________________________________________________________________
/*!

\class    genie::FGMBodekRitchie

\brief    The Bodek Richie Fermi Gass model. Implements the NuclearModelI
          interface.

\ref

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  October 09, 2004

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _FGM_BODEK_RITCHIE_H_
#define _FGM_BODEK_RITCHIE_H_

#include <map>

#include <TH1D.h>

#include "Physics/NuclearState/NuclearModelI.h"

using std::map;

namespace genie {

class FGMBodekRitchie : public NuclearModelI {

public:
  FGMBodekRitchie();
  FGMBodekRitchie(string config);
  virtual ~FGMBodekRitchie();

  using NuclearModelI::GenerateNucleon;  // inherit versions not overridden here
  using NuclearModelI::Prob;

  //-- implement the NuclearModelI interface
  bool           GenerateNucleon (const Target & t) const;
  double         Prob            (double mom, double w, const Target & t) const;
  NuclearModel_t ModelType       (const Target &) const
  {
    return kNucmFermiGas;
  }

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

  // Add "override" specifiers when GENIE moves to C++11
  // Note that, for this nuclear model, the removal energy is constant
  inline virtual double MinRemovalEnergy(const Target& t, double) const /*override*/
    { return GetRemovalEnergy(t); }
  inline virtual double MaxRemovalEnergy(const Target& t, double) const /*override*/
    { return GetRemovalEnergy(t); }

  inline virtual double MinMomentum(const Target&, double) const /*override*/
    { return 0.; }
  inline virtual double MaxMomentum(const Target&, double) const /*override*/
    { return fPMax; }

  double ProbDensity(double mom, double w, const Target& t, double r = 0.) const /*override*/;

private:
  void   LoadConfig (void);
  double GetRemovalEnergy (const Target& t) const;

  /// Helper private version of ProbDensity that also retrieves the width of
  /// the relevant bin
  double ProbDensity(double mom, double w, const Target& target,
    double r, double& bin_width) const;

  TH1D * ProbDistro (const Target & t) const;

  mutable map<string, TH1D *> fProbDistroMap;

  map<int, double> fNucRmvE;

  double fPMax;
  double fPCutOff;
  string fKFTable;
  bool fUseParametrization;
};

}         // genie namespace
#endif    // _FGM_BODEK_RITCHIE_H_

