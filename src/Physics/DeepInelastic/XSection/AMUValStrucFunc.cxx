//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Steven Gardiner <gardiner \at fnal.gov>
 Fermi National Accelerator Laboratory

 Huma Haider <huma.haider8 \at gmail.com>
 Aliargh Muslim University
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/DeepInelastic/XSection/AMUValStrucFunc.h"
#include "Physics/PartonDistributions/PDFModelI.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Framework/Utils/PhysUtils.h"

#include "LHAPDF/LHAPDF.h"

using namespace genie;

//____________________________________________________________________________
AMUValStrucFunc::AMUValStrucFunc() :
DISStructureFuncModelI("genie::AMUValStrucFunc")
{

}
//____________________________________________________________________________
AMUValStrucFunc::AMUValStrucFunc(string config) :
DISStructureFuncModelI("genie::AMUValStrucFunc", config)
{

}
//____________________________________________________________________________
AMUValStrucFunc::~AMUValStrucFunc()
{

}
//____________________________________________________________________________
void AMUValStrucFunc::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void AMUValStrucFunc::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void AMUValStrucFunc::LoadConfig(void)
{
  LOG("DISSF", pDEBUG) << "Loading configuration...";

  // Get a pointer to the owned PDF algorithm
  const PDFModelI* pdf_model = dynamic_cast< const PDFModelI* >(
    this->SubAlg("PDF-Set") );

  // Store the PDF set name from the configuration
  fPDFSetName = pdf_model->Id().Config();

  // Right now, this structure function calculation is only compatible
  // with LHAPDF6. If we're using a different PDF model, complain and exit.
  // TODO: Refactor Huma's code to use genie::PDFModelI instead of LHAPDF::PDF
  const LHAPDF6* lhapdf6 = dynamic_cast< const LHAPDF6* >( pdf_model );
  if ( !lhapdf6 ) {
    LOG("DISSF", pERROR) << "The AMUValDISStrucFunc algorithm is currently"
      " only compatible with LHAPDF6.";
    std::exit(1);
  }

  // Get a pointer to the right LHAPDF::PDF to use
  fPDF = lhapdf6->GetPDF();

  // Initialize the spectral function
  // TODO: understand how this works. Are there global variables/static members?
  start_spectral_function();

  LOG("DISSF", pDEBUG) << "Done loading configuration";
}
//____________________________________________________________________________
void AMUValStrucFunc::Calculate(const Interaction * interaction) const
{
  // Reset mutable members
  fF1 = 0.;
  fF2 = 0.;
  fF3 = 0.;
  fF4 = 0.;
  fF5 = 0.;
  fF6 = 0.;

  // Get process info & perform various checks
  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  const InitialState & init_state = interaction->InitState();
  const Target & tgt = init_state.Tgt();

  int  nuc_pdgc    = tgt.HitNucPdg();
  int  probe_pdgc  = init_state.ProbePdg();
  bool is_p        = pdg::IsProton       ( nuc_pdgc    );
  bool is_n        = pdg::IsNeutron      ( nuc_pdgc    );
  bool is_nu       = pdg::IsNeutrino     ( probe_pdgc  );
  bool is_nubar    = pdg::IsAntiNeutrino ( probe_pdgc  );
  bool is_lepton   = pdg::IsLepton       ( probe_pdgc  );
  bool is_dm       = pdg::IsDarkMatter   ( probe_pdgc  );
  bool is_CC       = proc_info.IsWeakCC();

  if ( !is_lepton && !is_dm ) return;
  if ( !is_p && !is_n       ) return;
  if ( tgt.N() == 0 && is_n ) return;
  if ( tgt.Z() == 0 && is_p ) return;

  if ( tgt.HitQrkIsSet() ) {

     int  qpdg = tgt.HitQrkPdg();

     bool is_u    = pdg::IsUQuark     (qpdg);
     bool is_ubar = pdg::IsAntiUQuark (qpdg);
     bool is_d    = pdg::IsDQuark     (qpdg);
     bool is_dbar = pdg::IsAntiDQuark (qpdg);
     bool is_s    = pdg::IsSQuark     (qpdg);
     bool is_sbar = pdg::IsAntiSQuark (qpdg);
     bool is_c    = pdg::IsCQuark     (qpdg);
     bool is_cbar = pdg::IsAntiCQuark (qpdg);

    // Make sure user inputs make sense
    if(is_nu    && is_CC && is_u   ) return;
    if(is_nu    && is_CC && is_c   ) return;
    if(is_nu    && is_CC && is_dbar) return;
    if(is_nu    && is_CC && is_sbar) return;
    if(is_nubar && is_CC && is_ubar) return;
    if(is_nubar && is_CC && is_cbar) return;
    if(is_nubar && is_CC && is_d   ) return;
    if(is_nubar && is_CC && is_s   ) return;
  }

  // Q^2 and Bjorken x
  double Q2val = this->Q2(interaction);
  double x = this->ScalingVar(interaction);

  //// Proton and neutron numbers of the target nucleus
  //int Z = interaction->InitState().Tgt().Z();
  //int N = interaction->InitState().Tgt().N();

  // Check whether we're working with an isoscalar target
  // TODO: revisit this
  bool isIso = true; //(Z == N);

//////////////// NEW

  // TODO: revisit this. What are imem and iset?
  int imem = 1; // lexical_cast<int>(smem);
  int iset = 2;

  //const genie::PDFModelI* pdf_model = dynamic_cast<const genie::PDFModelI*>(
  //  genie::AlgFactory::Instance()->GetAlgorithm("genie::LHAPDF6", "CT10"));

  // Numbers from param.h
  // TODO: remove these
  const double GeVtfm = 0.197330;
  const double Normal_Iso = 57.74176;
  const double Normal_NonIso = 59.07; // 27.40+31.671;

  double Q2_in_fm = Q2val / std::pow(GeVtfm, 2); // from GeV^2 to fm^2

  double f1nuc = 0;
  double f2nuc = 0;
  double f3nuc = 0;

  if ( Q2val >= 1. ) {
    if ( isIso ) {

      double resaf = 0.;

      spectralF1_Iso(x, Q2_in_fm, fPDFSetName, imem, iset, fPDF,resaf);
      f1nuc = resaf / Normal_Iso;

      spectralF2_Iso(x, Q2_in_fm, fPDFSetName, imem, iset, fPDF,resaf);
      f2nuc = resaf / Normal_Iso;

      spectralF3_Iso(x, Q2_in_fm, fPDFSetName, imem, iset, fPDF,resaf);
      f3nuc = resaf * x / Normal_Iso;
    }
    else {
      // Non-isoscalar case

      double resafp = 0.;
      double resafn = 0.;

      spectralF1_proton(x, Q2_in_fm, fPDFSetName, imem, iset, fPDF, resafp);
      spectralF1_neutron(x, Q2_in_fm, fPDFSetName, imem, iset, fPDF, resafn);
      f1nuc = (resafp + resafn) / Normal_NonIso;

      spectralF2_proton(x, Q2_in_fm, fPDFSetName, imem, iset, fPDF, resafp);
      spectralF2_neutron(x, Q2_in_fm, fPDFSetName, imem, iset, fPDF, resafn);
      f2nuc = (resafp + resafn) / Normal_NonIso;

      spectralF3_proton(x, Q2_in_fm, fPDFSetName, imem, iset, fPDF, resafp);
      spectralF3_neutron(x, Q2_in_fm, fPDFSetName, imem, iset, fPDF, resafn);
      f3nuc = (resafp + resafn) * x / Normal_NonIso; // if you want xf3
    }

    double f1shad = dshadowingF1(x, Q2_in_fm, fPDFSetName, imem, iset, fPDF);
    double f1pion = F1ApionEM(x, Q2val);

    fF1 = f1nuc + f1shad + f1pion;

    double f2shad = dshadowingF2(x, Q2_in_fm, fPDFSetName, imem, iset, fPDF);
    double f2pion = F2ApionEM(x, Q2val);

    fF2 = f2nuc + f2pion + f2shad;

    double f3shad = dshadowingf3(x, Q2_in_fm, fPDFSetName, imem, iset, fPDF)*x;

    fF3 = f3nuc + f3shad;
  }

//////////////// END NEW

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DISSF", pDEBUG)
     << "F1-F5 = "
     << fF1 << ", " << fF2 << ", " << fF3 << ", " << fF4 << ", " << fF5;
#endif

}
//____________________________________________________________________________
double AMUValStrucFunc::Q2(const Interaction * interaction) const
{
  // Return Q2 from the kinematics or, if not set, compute it from x,y
  const Kinematics& kinematics = interaction->Kine();

  // if Q2 (or q2) is set then prefer this value
  if (kinematics.KVSet(kKVQ2) || kinematics.KVSet(kKVq2)) {
    double Q2val = kinematics.Q2();
    return Q2val;
  }
  // If Q2 was not set, then compute it from x,y,Ev,Mnucleon
  if (kinematics.KVSet(kKVy)) {
    const InitialState & init_state = interaction->InitState();
    double Mn = init_state.Tgt().HitNucP4Ptr()->M(); // could be off-shell
    double x     = kinematics.x();
    double y     = kinematics.y();
    double Ev    = init_state.ProbeE(kRfHitNucRest);
    double Q2val = 2.*Mn*Ev*x*y;
    return Q2val;
  }
  LOG("DISSF", pERROR) << "Could not compute Q2!";
  return 0;
}
//____________________________________________________________________________
double AMUValStrucFunc::ScalingVar(const Interaction* interaction) const
{
  // The scaling variable is set to the normal Bjorken x.
  return interaction->Kine().x();
}
