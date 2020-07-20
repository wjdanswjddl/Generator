//____________________________________________________________________________
/*!

\class    genie::STAHadronTensorModel

\brief    Creates hadron tensor objects for calculations of inclusive
          cross sections using the short-time approximation (STA) approach

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  July 20, 2020

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _STA_HADRON_TENSOR_MODEL_H_
#define _STA_HADRON_TENSOR_MODEL_H_

// GENIE includes
#include "Physics/HadronTensors/TabulatedHadronTensorModelI.h"

namespace genie {

class STAHadronTensorModel : public TabulatedHadronTensorModelI {

public:

  STAHadronTensorModel();
  STAHadronTensorModel(std::string config);

  virtual ~STAHadronTensorModel();

protected:

  // Implementation of TabulatedHadronTensorModelI interface
  virtual HadronTensorI* ParseTensorFile( const std::string& full_file_name ) const;

};

} // namespace genie

#endif // _STA_HADRON_TENSOR_MODEL_H_
