// Dear emacs, this is -*- c++ -*-
/********************************************************
 * Class  : MIOCTPhase0TopoRoI                          *
 *                                                      *
 * Description: This class holds info about a RoI.      *
 * Designed for internal use in the package, but could  *
 * perhaps be useful also outside (e.g. in L1Topo sim)  *
 *                                                      *
 * Authors (alphabetical):                              *
 *      Christian Ohm (CERN)                            *
 ********************************************************/

#ifndef L1TOPO_MIOCTPHASE0TOPOROI_H
#define L1TOPO_MIOCTPHASE0TOPOROI_H

#include "L1TopoRDO/L1TopoRoI.h"

namespace L1Topo {

  class MIOCTPhase0TopoRoI : public L1Topo::L1TopoRoI {

   public:

    MIOCTPhase0TopoRoI() : L1TopoRoI(0, 0, 0,1), m_thrNumber(0), m_charge(100),
                           m_isFwd(false), m_sectorAddress(255) {};
    MIOCTPhase0TopoRoI(float eta, float phi, int thrNumber, int charge, bool isFwd, int sectorAddress)
    : L1TopoRoI(eta, phi, thrNumber, 1), m_thrNumber(thrNumber), m_charge(charge), m_isFwd(isFwd),
      m_sectorAddress(sectorAddress) {};
    ~MIOCTPhase0TopoRoI() {};

    // member functions
    int   thrNumber() const {return m_thrNumber;};
    int   charge() const {return m_charge;};
    bool  isFwd() const {return m_isFwd;};
    int   sectorAddress() const {return m_sectorAddress;};
    void  setThrNumber(int thrNumber) {m_thrNumber = thrNumber;};
    void  setCharge(int charge) {m_charge = charge;};
    void  setIsFwd(bool isFwd) {m_isFwd = isFwd;};
    void  setSectorAddress(int sa) { m_sectorAddress = sa; };
    void  print() {
      std::cout << "  Muon RoI: eta = " << m_eta << ", phi = " << m_phi << ", pT = " << m_pT << ", thrNumber = " << m_thrNumber << ", charge = " << m_charge << ", sectorAddress = " << m_sectorAddress <<
          ", isFwd = " << m_isFwd << std::endl;
    }

    // for sorting (in case there are more than we can send in a MIOCT)
    bool operator < (const MIOCTPhase0TopoRoI& roi) const {
      return (m_thrNumber < roi.thrNumber());
    }

   private:

    int   m_thrNumber;
    int   m_charge;
    bool  m_isFwd;
    int   m_sectorAddress;

  }; // class MIOCTPhase0TopoRoI

} // namespace L1Topo

#endif
