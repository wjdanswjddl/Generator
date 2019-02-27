//____________________________________________________________________________
/*!

\class    genie::NewQELXSec

\brief    Computes the Quasi Elastic (QEL) total cross section. \n
          Is a concrete implementation of the XSecIntegratorI interface. \n

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  February 26, 2019

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NEW_QEL_XSEC_H_
#define _NEW_QEL_XSEC_H_

#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/QuasiElastic/XSection/QELUtils.h"

#include "TMath.h"
#include "Math/IFunction.h"
#include "Math/Integrator.h"

namespace genie {

class NuclearModelI;
class PauliBlocker;
class VertexGenerator;

namespace utils {
  namespace gsl   {

    typedef enum EFullQELdXSecMode {
      kFreeNucleonVars,
      kFermiGasVars,
      kLocalFermiGasVars,
      kSpectralFuncVars,
      kAllQELVars
    } FullQELdXSecMode_t;

    class FullQELdXSec : public ROOT::Math::IBaseFunctionMultiDim
    {
     public:
       FullQELdXSec(const XSecAlgorithmI* xsec_model, const Interaction* interaction,
         FullQELdXSecMode_t mode, bool do_Pauli_blocking, const AlgId& pauli_blocker_ID,
         QELEvGen_BindingMode_t binding_mode, const AlgId& vertex_gen_ID);
       virtual ~FullQELdXSec();

       // ROOT::Math::IBaseFunctionMultiDim interface
       unsigned int NDim(void) const;
       double DoEval(const double* xin) const;
       ROOT::Math::IBaseFunctionMultiDim* Clone(void) const;

     private:
       const XSecAlgorithmI* fXSecModel;
       const NuclearModelI* fNuclModel;
       const PauliBlocker* fPauliBlocker;
       const VertexGenerator* fVertexGenerator;
       Interaction* fInteraction;
       FullQELdXSecMode_t fIntegrationMode;
       bool fDoPauliBlocking;
       QELEvGen_BindingMode_t fHitNucleonBindingMode;
    };

  } // gsl   namespace
} // utils namespace

class NewQELXSec : public XSecIntegratorI {

public:

  NewQELXSec(void);
  NewQELXSec(std::string config);

  /// XSecIntegratorI interface implementation
  double Integrate(const XSecAlgorithmI* model, const Interaction* i) const;

  /// Overload the Algorithm::Configure() methods to load private data
  /// members from configuration options
  void Configure(const Registry& config);
  void Configure(std::string config);


private:

  void LoadConfig (void);

  utils::gsl::FullQELdXSecMode_t
    IntegrationModeFromNuclearModel(const NuclearModelI&, const Target& tgt) const;

  // Configuration obtained from cross section model
  utils::gsl::FullQELdXSecMode_t fMode;
  QELEvGen_BindingMode_t fBindingMode;

  // XML configuration parameters
  std::string fGSLIntgType;
  double fGSLRelTol;
  unsigned int fGSLMaxEval;
  bool fPauliBlock;
  AlgId fPauliBlockID;
  AlgId fVertexGenID;

};


} // genie namespace

#endif  // _NEW_QEL_XSEC_H_
