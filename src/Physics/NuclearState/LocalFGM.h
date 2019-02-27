//____________________________________________________________________________
/*!

\class    genie::LocalFGM

\brief    local Fermi gas model. Implements the NuclearModelI
          interface.

\ref

\author   Joe Johnston, Steven Dytman

\created  December 2015

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _LOCAL_FGM_H_
#define _LOCAL_FGM_H_

#include <map>

#include <TH1D.h>

#include "Physics/NuclearState/NuclearModelI.h"

using std::map;

namespace genie {

class LocalFGM : public NuclearModelI {

public:
  LocalFGM();
  LocalFGM(string config);
  virtual ~LocalFGM();

  using NuclearModelI::GenerateNucleon;  // inherit versions not overridden here
  using NuclearModelI::Prob;

  //-- allow methods to be called with a radius
  bool   GenerateNucleon (const Target & t, double hitNucleonRadius) const;
  double Prob            (double p, double w, const Target & t,
			  double hitNucleonRadius) const;

  //-- implement the NuclearModelI interface
  bool GenerateNucleon (const Target & t) const {
    return GenerateNucleon(t,0.0);
  }
  double Prob (double p, double w, const Target & t) const {
    return Prob(p,w,t,0.0);
  }
  NuclearModel_t ModelType       (const Target &) const
  {
    return kNucmLocalFermiGas;
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
  TH1D * ProbDistro (const Target & t, double r) const;

  /// Helper private version of ProbDensity that also retrieves the width of
  /// the relevant bin
  double ProbDensity(double mom, double w, const Target& target,
    double r, double& bin_width) const;

  map<int, double> fNucRmvE;

  double fPMax;
};

}         // genie namespace
#endif    // _LOCAL_FGM_H_

