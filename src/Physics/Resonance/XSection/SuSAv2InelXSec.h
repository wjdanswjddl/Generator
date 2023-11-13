//____________________________________________________________________________
/*!

\class    genie::SuSAv2InelXSec

\brief    Computes the total inelastic (RES+DIS) cross section for
          the SuSAv2 treatment (SuSAv2InelPXSec). \n
          Is a concrete implementation of the XSecIntegratorI interface. \n

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  November 13, 2023

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _SUSAV2_INEL_XSEC_H_
#define _SUSAV2_INEL_XSEC_H_

#include "Physics/XSectionIntegration/XSecIntegratorI.h"

#include "TMath.h"
#include "Math/IFunction.h"
#include "Math/Integrator.h"

namespace genie {

namespace utils {

  namespace gsl   {

    class SuSAv2IneldXSec : public ROOT::Math::IBaseFunctionMultiDim {

     public:

       SuSAv2IneldXSec( const XSecAlgorithmI* xsec_model,
         const Interaction* interaction );

       virtual ~SuSAv2IneldXSec();

       // ROOT::Math::IBaseFunctionMultiDim interface
       unsigned int NDim() const;
       double DoEval( const double* xin ) const;
       ROOT::Math::IBaseFunctionMultiDim* Clone() const;

       Interaction* GetInteractionPtr();
       const Interaction& GetInteraction() const;

     private:

       const XSecAlgorithmI* fXSecModel;
       Interaction* fInteraction;

    };

  } // gsl   namespace

} // utils namespace

class SuSAv2InelXSec : public XSecIntegratorI {

public:

  SuSAv2InelXSec();
  SuSAv2InelXSec( std::string config );

  /// XSecIntegratorI interface implementation
  double Integrate( const XSecAlgorithmI* model, const Interaction* i ) const;

  /// Overload the Algorithm::Configure() methods to load private data
  /// members from configuration options
  void Configure( const Registry& config );
  void Configure( std::string config );

private:

  void LoadConfig();

  // XML configuration parameters
  std::string fGSLIntgType;
  double fGSLRelTol;
  unsigned int fGSLMaxEval;

};


} // genie namespace

#endif  // _SUSAV2_INEL_XSEC_H_
