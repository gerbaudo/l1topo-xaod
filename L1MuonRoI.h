#ifndef __L1MuonRoI__
#define __L1MuonRoI__
#include "TLorentzVector.h"
#include <stdint.h>

class L1MuonRoI {
public :

  TLorentzVector p4;
  uint32_t roiWord;
  float  thrValue;
  int isVetoed;

  //L1MuonRoI() : TLorentzVector(), roiWord(0), thrValue(0){}
  L1MuonRoI() : p4(), roiWord(0), thrValue(0),isVetoed(0) {}

enum RoISource {   Barrel,  Endcap, Forward };
enum Hemisphere {  Positive, Negative };
enum Charge {  Neg = 0, Pos = 1, Undef = 100  };
int getSectorAddress() const { return ( ( roiWord >> 14 ) & 0xff ); }
bool isFirstCandidate() const { return ( ( roiWord >> 22 ) & 0x1 ); }
bool isMoreCandInRoI() const {  return ( ( roiWord >> 1 ) & 0x1 ); }
bool isMoreCandInSector() const {return ( roiWord & 0x1 );  }
RoISource getSource() const {
      if( this->getSectorAddress() & 0x80 ) {
         return Endcap;
      } else if( this->getSectorAddress() & 0x40 ) {
         return Forward;
      } else {
         return Barrel;
      }
   }
Hemisphere getHemisphere() const {
      if( this->getSectorAddress() & 0x1 )return Positive;
      else return Negative;
   }
Charge getCharge() const {
      if( getSource() == Barrel ) return Undef;
      if( roiWord & 0x8000000 ) return Pos;
      else   return Neg;
   }
  bool isVetoedRoI() const { return ( roiWord & 0x10000000 ); }

ClassDef (L1MuonRoI,1);

};

#endif 
