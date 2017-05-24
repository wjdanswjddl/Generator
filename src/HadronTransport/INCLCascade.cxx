//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Steve Dytman <dytman+@pitt.edu>, Pittsburgh Univ.
         Marc Vololoniaina <narymarc@yahoo.com>, Univ. Antananarivo, Madagascar/Pittsburgh Univ.

         May 23, 2017

 For the class documentation see the corresponding header file.

 Important revisions after version 2.12.8 :
*/
//____________________________________________________________________________
//---------------------
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <ostream>
#include <sstream>
#include <cassert>
#include <cstdlib>
#include <map>
#include <cstdlib>
#include <sstream>

#include <TMath.h>

//#include "G4INCLNuclearMassTable.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLGlobals.hh" 


#include "G4INCLCascade.hh"
#include "G4INCLConfigEnums.hh"
#include "G4INCLParticle.hh"
#include "G4INCLIAvatar.hh"
#include "G4INCLIPropagationModel.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLStandardPropagationModel.hh"
#include "G4INCLRandom.hh"
#include "G4INCLRanecu.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLKinematicsUtils.hh"
#include "HINCLCascade.h"
#include "INCLCascade.h"
#include "G4INCLVersion.hh"
#include "G4INCLUnorderedVector.hh"
#include "G4INCLStore.hh"
#include "G4INCLIAvatar.hh"
#include "G4INCLEventInfo.hh"
#include "G4INCLIPropagationModel.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLStandardPropagationModel.hh"
#include "G4INCLRandom.hh"
#include "G4INCLRanecu.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLKinematicsUtils.hh"



// signal handler (for Linux and GCC)

#include "G4INCLConfig.hh"
#include "G4INCLVersion.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLUnorderedVector.hh"
#include "G4INCLParticle.hh"
#include "G4INCLEventInfo.hh"
#include "G4INCLStore.hh"
#include "G4INCLIAvatar.hh"

// signal handler (for Linux and GCC)
#include "G4INCLSignalHandling.hh"


// For I/O
#include "IWriter.hh"
#include "ASCIIWriter.hh"
#include "ProtobufWriter.hh"
#include "INCLTree.hh"
#include "ROOTWriter.hh"
#include "HDF5Writer.hh"

// For configuration
#include "G4INCLConfig.hh"
#include "ConfigParser.hh"

// For logging
#include "G4INCLLogger.hh"

// Generic de-excitation interface
#include "G4INCLIDeExcitation.hh"

// ABLA v3p de-excitation
#ifdef INCL_DEEXCITATION_ABLAXX
#include "G4INCLAblaInterface.hh"
#endif

// ABLA07 de-excitation
#ifdef INCL_DEEXCITATION_ABLA07
#include "G4INCLAbla07Interface.hh"
#endif

// SMM de-excitation
#ifdef INCL_DEEXCITATION_SMM
#include "G4INCLSMMInterface.hh"
#endif

// GEMINIXX de-excitation
#ifdef INCL_DEEXCITATION_GEMINIXX
#include "G4INCLGEMINIXXInterface.hh"
#endif

#ifdef HAS_BOOST_DATE_TIME
#include <boost/date_time/posix_time/posix_time.hpp>
namespace bpt = boost::posix_time;
#endif

#ifdef HAS_BOOST_TIMER
#include <boost/timer/timer.hpp>
namespace bt = boost::timer;
#endif

#include "INCLCascade.h"
#include <TMath.h>
#include <TRootIOCtor.h>
#include <TLorentzVector.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>

// --------------------------------------Include for GENIE---------------------
#include "Algorithm/AlgConfigPool.h"
#include "Algorithm/AlgFactory.h"
#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/Controls.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepParticle.h"
#include "HadronTransport/INCLCascade.h"
#include "HadronTransport/INukeHadroData.h"
#include "HadronTransport/INukeHadroFates.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Numerical/Spline.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGCodeList.h"
#include "Utils/PrintUtils.h"
#include "Utils/NuclearUtils.h"

#include "EVGCore/EventRecord.h"
#include "EVGDrivers/GMCJMonitor.h"
#include "EVGCore/EventRecordVisitorI.h"
#include "Ntuple/NtpWriter.h"
#include "Ntuple/NtpMCFormat.h"
#include "Utils/AppInit.h"
#include "Utils/StringUtils.h"
#include "Utils/RunOpt.h"
#include "PDG/PDGUtils.h"
#include "Apps/GenieItoINCL.h"
#include "HadronTransport/INukeHadroFates.h"
#include "HadronTransport/INukeUtils.h"

#include "INCLparticletype.hh"


using namespace genie;
using namespace genie::utils;
using namespace genie::utils::intranuke;
using namespace genie::constants;
using namespace genie::controls;
using namespace G4INCL;
using std::ostringstream;
using namespace std;

INCLCascade::INCLCascade() :
EventRecordVisitorI()
{

}
//___________________________________________________________________________
INCLCascade::INCLCascade(string name) :
EventRecordVisitorI(name)
{

}
//___________________________________________________________________________
INCLCascade::INCLCascade(string name, string config) :
EventRecordVisitorI(name, config)
{

}
//___________________________________________________________________________
INCLCascade::~INCLCascade()
{

}
// this is the call to INCL++ for hadron interactions
int INCLCascade::INCLcascade(int argc1, char *test[]) const {
  ConfigParser theParser;
  G4INCL::Config *theConfig = theParser.parse(argc1,test);
  if(!theConfig)
    return 0;

#ifdef INCL_SIGNAL_HANDLING
  enableSignalHandling();
#endif
  G4INCL::INCL *theINCLModel = new G4INCL::INCL(theConfig);

  int nShots = theConfig->getNumberOfShots();

  G4INCL::IDeExcitation *theDeExcitation = 0;
  switch(theConfig->getDeExcitationType()) {
#ifdef INCL_DEEXCITATION_ABLAXX
    case G4INCL::DeExcitationABLAv3p:
    theDeExcitation = new G4INCLAblaInterface(theConfig);
    break;
#endif
#ifdef INCL_DEEXCITATION_ABLA07
    case G4INCL::DeExcitationABLA07:
    theDeExcitation = new ABLA07CXX::Abla07Interface(theConfig);
    break;
#endif
#ifdef INCL_DEEXCITATION_SMM
    case G4INCL::DeExcitationSMM:
    theDeExcitation = new SMMCXX::SMMInterface(theConfig);
    break;
#endif
#ifdef INCL_DEEXCITATION_GEMINIXX
    case G4INCL::DeExcitationGEMINIXX:
    theDeExcitation = new G4INCLGEMINIXXInterface(theConfig);
    break;
#endif
    default:
    std::stringstream ss;
    ss << "########################################################\n"
    << "###              WARNING WARNING WARNING             ###\n"
    << "###                                                  ###\n"
    << "### You are running the code without any coupling to ###\n"
    << "###              a de-excitation model!              ###\n"
    << "###    Results will be INCOMPLETE and UNPHYSICAL!    ###\n"
    << "###    Are you sure this is what you want to do?     ###\n"
    << "########################################################\n";
    INCL_INFO(ss.str());
    std::cout << ss.str();
    break;
  }


  G4INCL::Random::SeedVector const theInitialSeeds = G4INCL::Random::getSeeds();
  /*
   *
   ***************GENIE**************************
   *
   */
  GHepRecord::SetPrintLevel(RunOpt::Instance()->EventRecordPrintLevel());
   //Initialize an NTuple writer to save GHEP record into a ROOT tree.
     //ntpw.CustomizeFilenamePrefix(kDefOptEvFilePrefix); //To do later for user specified options
  NtpWriter ntpwINCL(kNFGHEP, 1);
  ntpwINCL.Initialize();
  GMCJMonitor mcjmonitor(1);
   //
  int Countevent = 0;
  int Countevent_tr =0;
  GHepStatus_t    ist   = kIStInitialState;
  GHepStatus_t    ist1  = kIStStableFinalState;
  GHepStatus_t    ist2  = kIStFinalStateNuclearRemnant;

  int pdg_codeProbe = 0, pdg_codeP(0);
  double m_probe(0), m_pnP(0),E_pnP(0), EKinP(0);
  pdg_codeProbe =  INCLpartycleSpecietoPDGCODE(theConfig);
  for(int i = 0; i < nShots; i++) {

    G4INCL::Random::saveSeeds();

    INCL_DEBUG("Starting event " << i << std::endl);
    INCL_DEBUG("Random seeds " << G4INCL::Random::getSeeds() << std::endl);

    G4INCL::EventInfo result;
    result = theINCLModel->processEvent();
    result.event = i;




    if(theConfig->getProjectileSpecies().theType != Composite){
     m_probe = ParticleTable::getRealMass(theConfig->getProjectileSpecies().theType);}
     else {
       m_probe = ParticleTable::getRealMass(theConfig->getProjectileSpecies().theA,theConfig->getProjectileSpecies().theZ);}
       double m_target = ParticleTable::getTableMass(result.At, result.Zt);
       double E_probe = m_probe + theConfig->getProjectileKineticEnergy();
       double pz_p = TMath::Sqrt(TMath::Power(E_probe,2)-TMath::Power(m_probe,2));
       EventRecord * evrec = new EventRecord();
       Interaction * interaction = new Interaction;
       evrec->AttachSummary(interaction);
       TLorentzVector p4h   (0.,0.,pz_p/1000,E_probe/1000);
       TLorentzVector x4null(0.,0.,0.,0.);
       TLorentzVector p4tgt (0.,0.,0.,m_target/1000);
       int pdg_codeTarget= genie::pdg::IonPdgCode(theConfig->getTargetA(), theConfig->getTargetZ());

      //Add The Probe and the target to EventRecord 
       evrec->AddParticle(pdg_codeProbe, ist, -1,-1,-1,-1, p4h,x4null);
       evrec->AddParticle(pdg_codeTarget,ist,-1,-1,result.nParticles,-1,p4tgt,x4null);
       if (result.transparent==1) {
        Countevent_tr++; 
        evrec->AddParticle(pdg_codeProbe, ist1, 0,-1,-1,-1, p4h,x4null);
        evrec->AddParticle(pdg_codeTarget,ist1,1,-1,-1,-1,p4tgt,x4null);
        LOG("gevgen_hadron", pNOTICE ) << *evrec;
        ntpwINCL.AddEventRecord(i,evrec);
        delete evrec;
      }
        //Before INCL debug
      INCL_DEBUG("End of event " << i << std::endl);
      if( result.transparent ) {
        INCL_DEBUG("Transparent event" << std::endl);
      } else {
        INCL_DEBUG("Number of produced particles: " << result.nParticles << std::endl);
        if(theDeExcitation != 0) {
          theDeExcitation->deExcite(&result);
          INCL_DEBUG("Number of produced particles after de-excitation: " << result.nParticles << std::endl);
        }
      // Fill the inverse kinematics variables, if appropriate
        if(theConfig->getInverseKinematics()) {

          double gamma = G4INCL::KinematicsUtils::gammaFromKineticEnergy(theConfig->getProjectileSpecies(), theConfig->getProjectileKineticEnergy());
          result.fillInverseKinematics(gamma);
        }
        for (int nP = 0; nP < result.nParticles; nP++){
          pdg_codeP = INCLtopdgcode(result.A[nP],result.Z[nP]);
          EKinP = result.EKin[nP];
                   // double E_pnP = TMath::Sqrt((result.px[nP]*result.px[nP]+result.py[nP]*result.py[nP]+result.pz[nP]*result.pz[nP] +
           // TMath::Power(m_pnP,2))*TMath::Power(10,-6)) ;
      //      double E_pnP = TMath::Sqrt((result.px[nP]/1000)*(result.px[nP]/1000)+(result.py[nP]/1000)*(result.py[nP]/1000)+(result.pz[nP]/1000)*(result.pz[nP]/1000) +TMath::Power((m_pnP/1000),2)) ;          // Rotate lepton momentum vector from the reference frame (x'y'z') where z'is the Projectile direction, z'x':(theta plane) to the LAB
          TVector3 p3M(result.px[nP]/1000,result.py[nP]/1000,result.pz[nP]/1000);

           // Take a unit vector along the direction of the projectile.

          if (result.A[nP]>5) {
           m_pnP = ParticleTable::getTableMass(result.A[nP], result.Z[nP]);
           E_pnP= EKinP +m_pnP;
           TLorentzVector p4tgtf(p3M,E_pnP/1000);
           evrec->AddParticle(pdg_codeP,ist2,1,-1,-1,-1,p4tgtf,x4null);
         }else{
          m_pnP = 0.5*((result.px[nP])*(result.px[nP]) + (result.py[nP])*(result.py[nP]) + (result.pz[nP])*(result.pz[nP]) + EKinP*EKinP)/(EKinP);
          std::cout<<" Mass m_pnP = "<< m_pnP << std::endl;
          if (m_pnP<10)
          {
            std::cout<<" pdg_gamma = "<< m_pnP <<std::endl;
            pdg_codeP = 22;
            E_pnP = TMath::Sqrt((result.px[nP])*(result.px[nP])+(result.py[nP])*(result.py[nP])+(result.pz[nP])*(result.pz[nP])) ;  
          }else{
            double Mass_prodPar = PDGLibrary::Instance()->Find(pdg_codeP)->Mass();
            E_pnP = EKinP + Mass_prodPar*1000;
          }
          TLorentzVector p4tgtf(p3M,E_pnP/1000);
          evrec->AddParticle(pdg_codeP,ist1,0,-1,-1,-1,p4tgtf,x4null);
                      //evrec->AddParticle(pdg_codeP,kIStIntermediateState,1,-1,-1,-1,p4tgtf,x4null);
        }
      }
      LOG("gevgen_hadron", pNOTICE ) << *evrec;
      ntpwINCL.AddEventRecord(i,evrec);
         // mcjmonitor.Update(i,evrec);
      delete evrec;

    }


  }
// ____________________________GENIE SAVE NTUPLE 
  ntpwINCL.Save();

  if(theDeExcitation != 0) {
    delete theDeExcitation;
  }
  delete theINCLModel;
  return 0;
}

void INCLCascade::ProcessEventRecord(GHepRecord * evrec) const{
 ifstream inputincl;
 int arg =2;
 inputincl.open("inputincl.in");

 char* incl[] ={" ","inputincl.in"};

 fGMode = evrec->EventGenerationMode();
 if(fGMode == kGMdHadronNucleus ||
   fGMode == kGMdPhotonNucleus)
 {
   INCLCascade::INCLcascade(arg,incl);  
 }else if(fGMode == kGMdLeptonNucleus){
  INCLCascade::TransportHadrons(evrec);}
  inputincl.close();
}





bool INCLCascade::CanRescatter(const GHepParticle * p) const
{
// checks whether a particle that needs to be rescattered, can in fact be
// rescattered by this cascade MC

  assert(p);
  return  ( p->Pdg() == kPdgPiP     ||
    p->Pdg() == kPdgPiM     ||
    p->Pdg() == kPdgPi0     ||
    p->Pdg() == kPdgProton  ||
    p->Pdg() == kPdgNeutron 
    );
}

void INCLCascade::TransportHadrons(GHepRecord * evrec) const{
  int inucl = -1;
  if(fGMode == kGMdHadronNucleus ||
   fGMode == kGMdPhotonNucleus)
  {
   inucl = evrec->TargetNucleusPosition();
 }
 else
  if(fGMode == kGMdLeptonNucleus ||
   fGMode == kGMdNucleonDecay  ||
   fGMode == kGMdNeutronOsc)
  {
   inucl = evrec->RemnantNucleusPosition();
 }

 LOG("Intranuke", pNOTICE)
 << "Propagating hadrons within nucleus found in position = " << inucl;
 GHepParticle * nucl = evrec->Particle(inucl);
 if(!nucl) {
  LOG("Intranuke", pERROR)
  << "No nucleus found in position = " << inucl;
  LOG("Intranuke", pERROR)
  << *evrec;
  return;
}

fRemnA = nucl->A();
fRemnZ = nucl->Z();


LOG("Intranuke", pNOTICE)
<< "Nucleus (A,Z) = (" << fRemnA << ", " << fRemnZ << ")";

const TLorentzVector & p4nucl = *(nucl->P4());
fRemnP4 = p4nucl;


TObjArrayIter piter(evrec);
GHepParticle * p = 0;
G4INCL::Particle *p_incl = 0; //INCL Particle

int icurr = -1;

char * inclinit[] = {"NULL","-pp","-tFe56","-N1","-E10","-dabla07"};
int argcc = 6;
ConfigParser theParser;

G4INCL::Config *theConfig = theParser.parse(argcc,inclinit);
  //theConfig->setINCLXXDataFilePath("/home/nari/Documents/INCL_test/inclxx-v5.2.9.5-6aca7e6/data");
      #ifdef INCL_SIGNAL_HANDLING
enableSignalHandling();
      #endif
long lastAutoSaveEvent = -1;

  // Create the INCL- Model at the first Use.
if(!theINCLModel){
  theINCLModel = new G4INCL::INCL(theConfig);}

  int flags(1);
  bool is_DIS = evrec->Summary()->ProcInfo().IsDeepInelastic();
  bool is_RES = evrec->Summary()->ProcInfo().IsResonant();
  bool flagsRmncorr = false;
  bool DRRmnCorr = false;
  TLorentzVector * p_4 = nucl->P4();
 // momentum of the remnant nucleus.
  double pxRemn = p_4->Px();
  double pyRemn = p_4->Py();
  double pzRemn = p_4->Pz();
  int pdg_codeTargetRemn= genie::pdg::IonPdgCode(nucl->A(),nucl->Z());
  TLorentzVector p4tgf(p_4->Px(),p_4->Py(),p_4->Pz(),0.0);

 // Loop over GHEP and run intranuclear rescattering on handled particles
std::vector<G4INCL::Particle> pFsiListe;
std::vector<GHepParticle> copieListe;
std::vector<G4INCL::ParticleSpecies> SpecieListe; // Specie of the Particle needed by INCL++
std::vector<double> EKinList;
std::vector<double> ExciSum;
std::vector<G4INCL::ParticleSpecies> Speciep; // Specie of the Particle needed by INCL++
std::vector<int> pnumber;
  while( (p = (GHepParticle *) piter.Next()) )
  {
    icurr++;

    // Check whether the particle needs rescattering, otherwise skip it
    if( ! this->NeedsRescattering(p) ) continue;
    GHepParticle * sp = new GHepParticle(*p);

    // Set clone's mom to be the hadron that was cloned
    sp->SetFirstMother(icurr);

    // Check whether the particle can be rescattered
    if(!this->CanRescatter(sp)) {

       // if I can't rescatter it, I will just take it out of the nucleus
     LOG("Intranuke", pNOTICE)
     << "... Current version can't rescatter a " << sp->Name();
     sp->SetFirstMother(icurr);
     sp->SetStatus(kIStStableFinalState);
     evrec->AddParticle(*sp);
     delete sp;
     continue; // <-- skip to next GHEP entry
     }
      
      TLorentzVector *v4= sp->GetX4();

      ThreeVector thePosition(0.,0.,0.);
      ThreeVector momentum (0.,0.,0.);
      thePosition.setX(v4->X());
      thePosition.setY(v4->Y());
      thePosition.setZ(v4->Z());
      TLorentzVector * p4 = sp->P4();

      momentum.setX(p4->Px()*1000);
      momentum.setY(p4->Py()*1000);
      momentum.setZ(p4->Pz()*1000);

      int pdgc = sp->Pdg();

      const ParticleType theType = toINCLparticletype(pdgc); 

      double  E    = (p4->Energy())*1000;
      double massp = G4INCL::ParticleTable::getRealMass(theType);
      double M4 = momentum.getX()*momentum.getX() + momentum.getY()*momentum.getY() + momentum.getZ()*momentum.getZ();
      double EKin = E - massp;

      //************************idx of the particle 
      pnumber.push_back(icurr);

      EKinList.push_back(EKin);
        //ParticleType theType = G4INCL::Neutron;

      G4INCL::ParticleSpecies  theSpecies;
      theSpecies.theType=theType;
      theSpecies.theA=pdgcpiontoA(sp->Pdg());
      theSpecies.theZ=pdgcpiontoZ(sp->Pdg());
      
      
      Speciep.push_back(theSpecies);
      
      G4INCL::Particle *pincl=  new G4INCL::Particle( theType , E , momentum, thePosition);
      pFsiListe.push_back(*pincl);
      copieListe.push_back(*sp);
  
  //}  
      //delete sp;
  } //Ghep-entry
// Looop Over List of Hadron In the Nucleus and transport them using INCL++
  std::vector<G4INCL::EventInfo> ListeOfINCLresult;
  Double_t ExcitationE(0);
   for (int it=0; it<pFsiListe.size();it++){
      //*******************************Desexitation************************************
p_incl =new Particle(pFsiListe.at(it));
GHepParticle *sp = new GHepParticle(copieListe.at(it));
 // Model->Process(Specie, pincl, EKin, A,Z);
      G4INCL::IDeExcitation *theDeExcitation = 0;
    switch(theConfig->getDeExcitationType()) {
#ifdef INCL_DEEXCITATION_ABLAXX
        case G4INCL::DeExcitationABLAv3p:
        theDeExcitation = new G4INCLAblaInterface(theConfig);
        break;
#endif
#ifdef INCL_DEEXCITATION_ABLA07
        case G4INCL::DeExcitationABLA07:
        theDeExcitation = new ABLA07CXX::Abla07Interface(theConfig);
        break;
#endif
#ifdef INCL_DEEXCITATION_SMM
        case G4INCL::DeExcitationSMM:
        theDeExcitation = new SMMCXX::SMMInterface(theConfig);
        break;
#endif
#ifdef INCL_DEEXCITATION_GEMINIXX
        case G4INCL::DeExcitationGEMINIXX:
        theDeExcitation = new G4INCLGEMINIXXInterface(theConfig);
        break;
#endif
        default:
        std::stringstream ss;
        ss << "########################################################\n"
        << "###              WARNING WARNING WARNING             ###\n"
        << "###                                                  ###\n"
        << "### You are running the code without any coupling to ###\n"
        << "###              a de-excitation model!              ###\n"
        << "###    Results will be INCOMPLETE and UNPHYSICAL!    ###\n"
        << "###    Are you sure this is what you want to do?     ###\n"
        << "########################################################\n";
        INCL_INFO(ss.str());
        std::cout << ss.str();
        break;
      }
      //*******************************************************************************
      
      std::stringstream progressSS;
     #ifdef HAS_BOOST_DATE_TIME
      bpt::ptime lapTime, endTime, etaTime;
      bpt::ptime startTime = bpt::microsec_clock::local_time();
      INCL_INFO("Event loop started at local time: " << startTime << std::endl);
      #endif
        #ifdef HAS_BOOST_TIMER
      std::stringstream occupationSS;
      bt::cpu_timer theTimer;
      #endif

          G4INCL::Random::SeedVector const theInitialSeeds = G4INCL::Random::getSeeds();
      int i =0;
            G4INCL::Random::saveSeeds();
      INCL_DEBUG("Random seeds " << G4INCL::Random::getSeeds() << std::endl);

      G4INCL::EventInfo result;
      result=theINCLModel->processEvent(Speciep[it],p_incl,EKinList[it],nucl->A(),nucl->Z());
      ListeOfINCLresult.push_back(result); 
      result.event = i;

      double kinE(0), m_P(0),E_pnP(0);
      int pdg_codeP(0);
      TLorentzVector x4null(0.,0.,0.,0.);
        //TLorentzVector p4tgtf(0.,0.,0.,0.);
      

      INCL_DEBUG("End of event " << i << std::endl);
      //INCL++ can generate transparent event
      std::cout<<" nParticle before INCL_debug "<< result.nParticles <<std::endl;
      if(result.transparent){     
        sp->SetStatus(kIStStableFinalState);
       evrec->AddParticle(*sp);
       evrec->Particle(sp->FirstMother())->SetRescatterCode(1);
      if(it==(pFsiListe.size()-1)&&(is_DIS==true||is_RES==true)){
        if(ExcitationE==0){
          GHepParticle remnantRES(pdg_codeTargetRemn, kIStFinalStateNuclearRemnant, inucl,-1,-1,-1, fRemnP4, x4null);
          evrec->AddParticle(remnantRES);
        }else{
        ListeOfINCLresult.at(it-1).EStarRem[i]=ExcitationE;
        theDeExcitation->deExcite(&ListeOfINCLresult.at(it-1));
        for(int nP = 0; nP < ListeOfINCLresult.at(it-1).nParticles; ++nP){
                        kinE = ListeOfINCLresult.at(it-1).EKin[nP];

          pdg_codeP = INCLtopdgcode(ListeOfINCLresult.at(it-1).A[nP],ListeOfINCLresult.at(it-1).Z[nP]);
          if(nP == ListeOfINCLresult.at(it-1).nParticles-1){ 
           m_P = ParticleTable::getTableMass(ListeOfINCLresult.at(it-1).A[nP],ListeOfINCLresult.at(it-1).Z[nP] );
           E_pnP = TMath::Sqrt((ListeOfINCLresult.at(it-1).px[nP]*ListeOfINCLresult.at(it-1).px[nP]+ListeOfINCLresult.at(it-1).py[nP]*ListeOfINCLresult.at(it-1).py[nP]+ListeOfINCLresult.at(it-1).pz[nP]*ListeOfINCLresult.at(it-1).pz[nP]+TMath::Power(m_P,2) )*TMath::Power(10,-6));
           TLorentzVector p4REM(ListeOfINCLresult.at(it-1).px[nP]/1000 + pxRemn,ListeOfINCLresult.at(it-1).py[nP]/1000+pyRemn,ListeOfINCLresult.at(it-1).pz[nP]/1000+pzRemn,E_pnP);
           evrec->AddParticle(pdg_codeP, kIStFinalStateNuclearRemnant,inucl,-1,-1,-1,p4REM,x4null);
          }
          else if(nP<ListeOfINCLresult.at(it-1).nParticles-1) {
           m_P = 0.5*((ListeOfINCLresult.at(it-1).px[nP])*(ListeOfINCLresult.at(it-1).px[nP]) + (ListeOfINCLresult.at(it-1).py[nP])*(ListeOfINCLresult.at(it-1).py[nP]) + (ListeOfINCLresult.at(it-1).pz[nP])*(ListeOfINCLresult.at(it-1).pz[nP]) + kinE*kinE)/(kinE);
           if (m_P<8){
            pdg_codeP = 22;
            E_pnP = TMath::Sqrt((ListeOfINCLresult.at(it-1).px[nP])*(ListeOfINCLresult.at(it-1).px[nP])+(ListeOfINCLresult.at(it-1).py[nP])*(ListeOfINCLresult.at(it-1).py[nP])+(ListeOfINCLresult.at(it-1).pz[nP])*(ListeOfINCLresult.at(it-1).pz[nP])) ;
          }else{
            double Mass_prodPar = PDGLibrary::Instance()->Find(pdg_codeP)->Mass();
            E_pnP = kinE + Mass_prodPar*1000;
          }
          TLorentzVector p4Pproduced(ListeOfINCLresult.at(it-1).px[nP]/1000,ListeOfINCLresult.at(it-1).py[nP]/1000,ListeOfINCLresult.at(it-1).pz[nP]/1000,E_pnP/1000);                                    
          evrec->AddParticle(pdg_codeP, kIStStableFinalState, pnumber.at(it-1), inucl,-1,-1,p4Pproduced,x4null);
        }
          
          }

       }
      //pdg_codeP=INCLtopdgcode(p_incl->getA(),p_incl->getA());
      //TLorentzVector p4transparent(p_incl->getMomentum().getX()/1000,p_incl->getMomentum().getY()/1000,p_incl->getMomentum().getZ()/1000,p_incl->getEnergy()/1000);
      //evrec->AddParticle(pdg_codeP, kIStStableFinalState,inucl,-1,-1,-1,p4transparent,x4null);     
      }  // Else tsy misy Excitation 
      }else {
         INCL_DEBUG("Number of produced particles: " << result.nParticles << std::endl);

  if(is_DIS==true)  {
     ExcitationE += result.EStarRem[i];
    result.EStarRem[i]=ExcitationE;
  }else if(is_RES==true) { 
     ExcitationE += result.EStarRem[i];
    result.EStarRem[i]=ExcitationE;
}else{
          theDeExcitation->deExcite(&result);
          INCL_DEBUG("Number of produced particles after de-excitation: " << result.nParticles << std::endl);
}
if (result.nParticles==1&&it!=(pFsiListe.size()-1))
{
  TLorentzVector gammad(0.,0.,0.,0.);
  evrec->AddParticle(22, kIStStableFinalState, icurr, pnumber.at(it),-1,-1,gammad,x4null);
}
if(theDeExcitation != 0 && it==(pFsiListe.size()-1)&&(is_DIS==true || is_RES==true)) {
          theDeExcitation->deExcite(&result);
          INCL_DEBUG("Number of produced particles after de-excitation: " << result.nParticles << std::endl);
}

std::cout<<"  icurr = "<< pnumber.at(it) << " number of produced particle = "<<result.nParticles<<std::endl;
        for(int nP = 0; nP < result.nParticles; ++nP){

          kinE = result.EKin[nP];

          pdg_codeP = INCLtopdgcode(result.A[nP],result.Z[nP]);
          if(nP==result.nParticles-1&&it==(pFsiListe.size()-1)){ 
           m_P = ParticleTable::getTableMass(result.A[nP],result.Z[nP] );
           E_pnP = TMath::Sqrt((result.px[nP]*result.px[nP]+result.py[nP]*result.py[nP]+result.pz[nP]*result.pz[nP]+TMath::Power(m_P,2) )*TMath::Power(10,-6));
           TLorentzVector p4tgtf(result.px[nP]/1000 + pxRemn,result.py[nP]/1000+pyRemn,result.pz[nP]/1000+pzRemn,E_pnP);
           evrec->AddParticle(pdg_codeP, kIStFinalStateNuclearRemnant,inucl,-1,-1,-1,p4tgtf,x4null);
           flagsRmncorr =true;
         }
         else if(nP<result.nParticles-1) {
           m_P = 0.5*((result.px[nP])*(result.px[nP]) + (result.py[nP])*(result.py[nP]) + (result.pz[nP])*(result.pz[nP]) + kinE*kinE)/(kinE);
           if (m_P<8){
            pdg_codeP = 22;
            E_pnP = TMath::Sqrt((result.px[nP])*(result.px[nP])+(result.py[nP])*(result.py[nP])+(result.pz[nP])*(result.pz[nP])) ;
          }else{
            double Mass_prodPar = PDGLibrary::Instance()->Find(pdg_codeP)->Mass();
            E_pnP = kinE + Mass_prodPar*1000;
          }
          TLorentzVector p4Pproduced(result.px[nP]/1000,result.py[nP]/1000,result.pz[nP]/1000,E_pnP/1000);                                    
          evrec->AddParticle(pdg_codeP, kIStStableFinalState, pnumber.at(it), inucl,-1,-1,p4Pproduced,x4null);
        }

  }

} 

       GHepParticle remnant(pdg_codeTargetRemn, kIStFinalStateNuclearRemnant, inucl,-1,-1,-1, fRemnP4, x4null);
    GHepParticle remnant_nucleus(kPdgHadronicBlob, kIStFinalStateNuclearRemnant, inucl,-1,-1,-1, p4tgf, x4null);

    if (flagsRmncorr==false&&is_RES==false&&is_DIS==false){         
        evrec->AddParticle(remnant);
    }//else if((is_RES==true || is_DIS==true)&&result.nParticles==0){

       // if(DRRmnCorr==false) {
         //           DRRmnCorr = true;
        //}else {
          //evrec->AddParticle(remnant);
          //flagsRmncorr=true;          
        ///}
      //}
  
}
evrec->Particle(inucl)->SetStatus(kIStIntermediateState);  

} 

int INCLCascade::INCLtopdgcode(int A, int Z)const {
  int pdg_codeP(0);
  if (A==1 && Z==1)       return pdg_codeP = 2212;
  else if(A==1 && Z==0)   return pdg_codeP = 2112;
  else if ( A==0 && Z==0) return pdg_codeP=111;
  else if (A==0 && Z==1)  return pdg_codeP = 211;
  else if (A==0 && Z==-1) return pdg_codeP= -211;
  else return pdg_codeP = genie::pdg::IonPdgCode( A , Z );
}
  //____________________________________________________________________________________
int INCLCascade::pdgcpiontoA(int pdgc)const{
  int A(0);
  if(pdgc == 2212 || pdgc==2112) return A=1;
  else if(pdgc == 211|| pdgc == -211 || pdgc == 111) return A=0;
}
  //_____________________________________________________________________________________
int INCLCascade::pdgcpiontoZ(int pdgc)const{
  int Z(0);
  if(pdgc == 2212 || pdgc==211) return Z=1;
  else if (pdgc==2112|| pdgc == 111) return Z=0;
  else if (pdgc==-211) return Z=-1;
}
  //_____________________________________________________________________________________
bool INCLCascade::NeedsRescattering(const GHepParticle * p) const
{
    // checks whether the particle should be rescattered
  assert(p);
    // attempt to rescatter anything marked as 'hadron in the nucleus'
  return (p->Status() == kIStHadronInTheNucleus);
}
