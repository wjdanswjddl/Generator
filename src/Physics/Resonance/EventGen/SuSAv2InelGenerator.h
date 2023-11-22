//____________________________________________________________________________
/*!

\class    genie::SuSAv2InelGenerator

\brief    Generates inelastic (RES+DIS) events according to the SuSAv2 model.
          Is a concrete implementation of the EventRecordVisitorI interface.

\ref      J. González-Rosa et al., Phys. Rev. D 105, 093009 (2022)
          https://doi.org/10.1103/PhysRevD.105.093009
          https://arxiv.org/abs/2203.12308

\author   Jesús González-Rosa <jgrosa \at us.es>
          University of Seville

          Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  21 November 2023

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _SUSAV2_INEL_GENERATOR_H_
#define _SUSAV2_INEL_GENERATOR_H_

#include "Framework/Utils/Range1.h"
#include "Physics/Common/KineGeneratorWithCache.h"

namespace genie {

class SuSAv2InelGenerator : public KineGeneratorWithCache {

public :

  SuSAv2InelGenerator();
  SuSAv2InelGenerator( std::string config );
 ~SuSAv2InelGenerator();

  // Implement the EventRecordVisitorI interface
  void ProcessEventRecord( GHepRecord* evrec ) const;

  // Overload the Algorithm::Configure() methods to load private data members
  // from the configuration options
  void Configure( const Registry& config );
  void Configure( std::string config );

private:

  void LoadConfig();

  void SelectLeptonKinematics( GHepRecord* evrec ) const;

  double ComputeMaxXSec( const Interaction* interaction ) const;

  double fWcut; ///< Wcut parameter in DIS/RES joining scheme
  double fSafetyFactor;
};

}      // genie namespace
#endif // _SUSAV2_INEL_GENERATOR_H_
