//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Steven Gardiner <gardiner \at fnal.gov>
         Fermi National Accelerator Laboratory

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

// GENIE includes
#include "Physics/HadronTensors/STAHadronTensorModel.h"
#include "Physics/HadronTensors/TabulatedHadronTensorModelI.h"
#include "Physics/HadronTensors/TabulatedLabFrameHadronTensor.h"

//____________________________________________________________________________
genie::STAHadronTensorModel::STAHadronTensorModel()
  : genie::TabulatedHadronTensorModelI("genie::STAHadronTensorModel")
{

}

//____________________________________________________________________________
genie::STAHadronTensorModel::STAHadronTensorModel(std::string config)
  : genie::TabulatedHadronTensorModelI("genie::STAHadronTensorModel", config)
{

}

//____________________________________________________________________________
genie::STAHadronTensorModel::~STAHadronTensorModel()
{

}

//____________________________________________________________________________
genie::HadronTensorI* genie::STAHadronTensorModel::ParseTensorFile(
  const std::string& full_file_name) const
{
  return new TabulatedLabFrameHadronTensor( full_file_name );
}
