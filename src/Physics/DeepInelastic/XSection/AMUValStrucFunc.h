//____________________________________________________________________________
/*!

\class    genie::AMUValStrucFunc

\brief    Implements the treatment of DIS structure functions developed
          by researchers at Aliargh Muslim University and the University
          of Valencia.
          Concrete implementation of the DISStructureFuncModelI interface.

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

          Huma Haider <huma.haider8 \at gmail.com>
          Aliargh Muslim University

\created  Feb 20, 2020

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _AMUVAL_DIS_SF_H_
#define _AMUVAL_DIS_SF_H_

#include "Physics/DeepInelastic/XSection/DISStructureFuncModelI.h"
#include "Framework/Interaction/Interaction.h"
#include "Physics/PartonDistributions/PDF.h"

// Declare external functions from Huma's code
// TODO: organize this in a better way (make a namespace, perhaps?)

// Hide these function declarations from rootcint / rootcling. We won't need them
// when generating ROOT dictionaries, and they can cause trouble.
#ifndef __MAKECINT__
#include "Physics/PartonDistributions/LHAPDF6.h"
void start_spectral_function();

// Helper functions for F1
void spectralF1_Iso(double& x, double& q2, const std::string& setname, int& imem, int& iset,
  const LHAPDF::PDF* pdf, double& resaf);
void spectralF1_proton(double& x, double& q2, const std::string& setname, int& imem, int& iset,
  const LHAPDF::PDF* pdf, double& resafp);
void spectralF1_neutron(double& x, double& q2, const std::string& setname, int& imem, int& iset,
  const LHAPDF::PDF* pdf, double& resafn);
double dshadowingF1(double x, double q2,std::string setname, int imem,int iset,
  const LHAPDF::PDF* pdf);
double F1ApionEM(double x, double q2);

// Helper functions for F2
void spectralF2_Iso(double& x, double& q2, const std::string& setname, int& imem, int& iset,
  const LHAPDF::PDF* pdf, double& resaf);
void spectralF2_proton(double& x, double& q2, const std::string& setname, int& imem, int& iset,
  const LHAPDF::PDF* pdf, double& resafp);
void spectralF2_neutron(double& x, double& q2, const std::string& setname, int& imem, int& iset,
  const LHAPDF::PDF* pdf, double& resafn);
double dshadowingF2(double x, double q2, std::string setname, int imem, int iset,
  const LHAPDF::PDF* pdf);
double F2ApionEM(double x, double q2);

// Helper functions for F3
void spectralF3_Iso(double& x, double& q2, const std::string& setname, int& imem, int& iset,
  const LHAPDF::PDF* pdf, double& resaf);
void spectralF3_proton(double& x, double& q2, const std::string& setname, int& imem, int& iset,
  const LHAPDF::PDF* pdf, double& resafp);
void spectralF3_neutron(double& x, double& q2, const std::string& setname, int& imem, int& iset,
  const LHAPDF::PDF* pdf, double& resafn);
double dshadowingf3(double x, double q2, std::string setname, int imem, int iset,
  const LHAPDF::PDF* pdf);
#endif

//------------------------------------------------------------------------------
namespace genie {

class AMUValStrucFunc : public DISStructureFuncModelI {

public:

  AMUValStrucFunc();
  AMUValStrucFunc(string config);

  virtual ~AMUValStrucFunc();

  // Common code for all DISFormFactorsModelI interface implementations
  virtual double F1 (void) const /*override*/ { return fF1; }
  virtual double F2 (void) const /*override*/ { return fF2; }
  virtual double F3 (void) const /*override*/ { return fF3; }
  virtual double F4 (void) const /*override*/ { return fF4; }
  virtual double F5 (void) const /*override*/ { return fF5; }
  virtual double F6 (void) const /*override*/ { return fF6; }

  virtual void Calculate(const Interaction* interaction) const;

  // Overload Algorithm's Configure() to set data members
  // from the configuration registry
  void Configure(const Registry& config);
  void Configure(std::string param_set);

protected:

  virtual void LoadConfig (void);

  /// Helper function that retrieves (or calculates) Q^2
  double Q2(const Interaction* interaction) const;

  /// Helper function that retrieves (or calculates) Bjorken x
  double ScalingVar(const Interaction* interaction) const;

  // Structure function values (set by Calculate())
  mutable double fF1;
  mutable double fF2;
  mutable double fF3;
  mutable double fF4;
  mutable double fF5;
  mutable double fF6;

  const LHAPDF::PDF* fPDF;
  std::string fPDFSetName;
};

}      // genie namespace
#endif // _AMUVAL_DIS_SF_H_
