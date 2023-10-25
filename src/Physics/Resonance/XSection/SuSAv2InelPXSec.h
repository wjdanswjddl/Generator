//____________________________________________________________________________
/*!

\class    genie::SuSAv2InelPXSec

\brief    Computes the differential cross section for inelastic (RES + DIS)
          electron and CC neutrino scattering according to the SuSAv2 treatment.

          Is a concrete implementation of the XSecAlgorithmI interface.

\ref      J. González-Rosa et al., Phys. Rev. D 105, 093009 (2022)
          https://doi.org/10.1103/PhysRevD.105.093009
          https://arxiv.org/abs/2203.12308

\author   Jesús González-Rosa <jgrosa \at us.es>
          University of Seville

          Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  10 October 2023

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _SUSAV2_INEL_PXSEC_H_
#define _SUSAV2_INEL_PXSEC_H_

// Standard library includes
#include <string>

// GENIE includes
#include "Framework/EventGen/XSecAlgorithmI.h"
//#include "Framework/ParticleData/BaryonResonance.h"

namespace genie {

class SuSAv2InelPXSec : public XSecAlgorithmI {

public:

  SuSAv2InelPXSec();
  SuSAv2InelPXSec( std::string config );
  virtual ~SuSAv2InelPXSec();

  // Implement the XSecAlgorithmI interface
  double XSec         ( const Interaction* i, KinePhaseSpace_t k ) const override;
  double Integral     ( const Interaction* i ) const override;
  bool   ValidProcess ( const Interaction* i ) const override;

  // Overload the Algorithm::Configure() methods to load private data members
  // from configuration options
  void Configure( const Registry& config ) override;
  void Configure( std::string config ) override;

private:

  /// Helper function for the Configure() methods
  void LoadConfig();

  /// Loads pre-computed inelastic structure functions from a file
  void LoadStructureFunctions( const std::string& input_file_name );

  // Storage for tabulated inelastic structure functions
  // TODO: refactor to use GENIE's bilinear interpolation code
  std::vector< double > Q2vec;
  std::vector< double > Wvec;
  std::vector< double > w1pvec;
  std::vector< double > w1nvec;
  std::vector< double > w2pvec;
  std::vector< double > w2nvec;
  std::vector< double > w3pvec;
  std::vector< double > w3nvec;

  double pF,shifte,esep,xW;
   double a1, a2, a3, a4;
  double b1, b2, b3, b4, b5, b6;
  double qi0, qi1, qi00, qi11, w0;
};

}       // genie namespace

#endif  // _SUSAV2_INEL_PXSEC_H_
