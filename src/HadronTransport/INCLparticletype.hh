#ifndef INCLparticletype_hh
#define INCLparticletype_hh 1

#include "G4INCLCascade.hh"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include <list>
#include <sstream>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
using namespace genie;
namespace G4INCL{

//TTree *tree;
//TFile *outfile;	
G4INCL::INCL *theINCLModel;	
	int INCLpartycleSpecietoPDGCODE(G4INCL::Config *theConfig){
		int pdg_codeProbe(0);
		if(theConfig->getProjectileSpecies().theType != Composite){
			if(ParticleTable::getName(theConfig->getProjectileSpecies().theType)=="pi0")   return pdg_codeProbe =111;
			else if(ParticleTable::getName(theConfig->getProjectileSpecies().theType)=="pi+") return  pdg_codeProbe =211;
			else if(ParticleTable::getName(theConfig->getProjectileSpecies().theType)=="pi-") return pdg_codeProbe =-211;         
		}
		if(theConfig->getProjectileSpecies().theA==1&&theConfig->getProjectileSpecies().theZ==1) return pdg_codeProbe = 2212;  
		else if(theConfig->getProjectileSpecies().theA ==1 && theConfig->getProjectileSpecies().theZ==0) return pdg_codeProbe = 2112;
	}
	int INCLtopdgcode(int A, int Z){
  int pdg_codeP(0);
  if (A==1 && Z==1)       return pdg_codeP = 2212;
  else if(A==1 && Z==0)   return pdg_codeP = 2112;
  else if ( A==0 && Z==0) return pdg_codeP=111;
  else if (A==0 && Z==1)  return pdg_codeP = 211;
  else if (A==0 && Z==-1) return pdg_codeP= -211;
  else return pdg_codeP = genie::pdg::IonPdgCode( A , Z );
}
	G4INCL::ParticleType toINCLparticletype(int pdgc){

		if       (pdgc == 2212)         return     G4INCL::Proton;
		else if(pdgc == 2112)         return     G4INCL::Neutron;
		else if(pdgc == 211)          return     G4INCL::PiPlus;
		else if(pdgc == -211)         return     G4INCL::PiMinus;
		else if(pdgc == 111)          return     G4INCL::PiZero;
		else if(pdgc == 1000010020)     return   G4INCL::Composite;
		else if(pdgc == 1000010030)      return  G4INCL::Composite;
		else if(pdgc == 1000020030)      return  G4INCL::Composite;
		else if(pdgc == 1000020040)      return  G4INCL::Composite;
		else                             return  G4INCL::UnknownParticle;
	}
GHepParticle *INCLtoGenieParticle(G4INCL::EventInfo result, int nP, GHepStatus_t    ist, int mom1, int mom2) {
	int pdg_codeP(0);
	double E_pnP(0), EKinP = result.EKin[nP];
	TVector3 p3M(result.px[nP]/1000,result.py[nP]/1000,result.pz[nP]/1000);
	TLorentzVector x4null(0.,0.,0.,0.);
double 	m_pnP = 0.5*((result.px[nP])*(result.px[nP]) + (result.py[nP])*(result.py[nP]) + (result.pz[nP])*(result.pz[nP]) + EKinP*EKinP)/(EKinP);
if (m_pnP<10&&result.A[nP]==0&&result.Z[nP]==0)
          {
            pdg_codeP = 22;
            E_pnP = TMath::Sqrt((result.px[nP])*(result.px[nP])+(result.py[nP])*(result.py[nP])+(result.pz[nP])*(result.pz[nP])) ;  
          }else{
          	pdg_codeP  = INCLtopdgcode(result.A[nP],result.Z[nP]);
            double Mass_prodPar = PDGLibrary::Instance()->Find(pdg_codeP)->Mass();
            E_pnP = EKinP + Mass_prodPar*1000;
          }
          TLorentzVector p4tgtf(p3M,E_pnP/1000);
          return new GHepParticle(pdg_codeP,ist,mom1,mom2,-1,-1,p4tgtf,x4null); 
}   
  /*G4INCL::Particle *GenieparticletoINCLparticle(GHepParticle * p, G4INCL::ParticleType type) {
    TLorentzVector *v4= p->GetX4();;
    ThreeVector thePosition(0,0,0.);
    ThreeVector momentum (0.,0.,0.);
     thePosition.setX(v4->X());
     thePosition.setY(v4->Y());
     thePosition.setZ(v4->Z());
    // G4INCL::ParticleType type = GenieINCLInterface::GenieparticletoINCLParticleType(p->Pdg());
      // Creat a particle for INCL
      TLorentzVector * p4 = p->P4();
      momentum.setX(p4->Px()*1000);
      momentum.setY(p4->Py()*1000);
      momentum.setZ(p4->Pz()*1000);
      double  E    = p4->Energy()*1000;

      //G4INCL::Particle *pincl = new Particle( type , E , momentum, thePosition);

 return  new G4INCL::Particle(type, E, momentum, thePosition);
 //return pincl;  
  }
*/

}
#endif 
