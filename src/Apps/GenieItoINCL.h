#ifndef _Genie_toINCL_H_
#define _Genie_toINCL_H_
#include <string>
#include "PDG/PDGUtils.h"

string InputProbeINCL(int pdg_probe){
string probe;
if (pdg_probe==2212){
  probe="-pp";
}else if(pdg_probe==2112){
  probe="-pn";
}
else if(pdg_probe==111){
  probe="-ppizero";
}
else if(pdg_probe==211){
  probe="-ppi+";
}
else if(pdg_probe==-211){
  probe="-ppi-";
}
return probe;
}

string InputtargetINCL(int pdg_target){
string target;
int At=genie::pdg::IonPdgCodeToA(pdg_target);
int Zt=genie::pdg::IonPdgCodeToZ(pdg_target);
if(Zt==26){
  target="-tFe56";
}else if (Zt==6) {
  target="-tC12";
}else if (Zt==8) {
  target="-tO16";
}else if (Zt==4){
target="-tBe9";
}else if (Zt==3){
target="-tLi7";
}else if(Zt==7){
target="-tN14";
}else if(Zt==9){
target="-tF19";
}else if(Zt==10){
target="-tNe20";
}
return target;
}
#endif
