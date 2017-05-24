#ifndef INCLparticletype_hh
#define INCLparticletype_hh 1

#include "G4INCLCascade.hh"
#include "GHEP/GHepParticle.h"
#include <list>
#include <sstream>
#include <TFile.h>
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
