// Dear emacs, this is -*- c++ -*-
/********************************************************
 * Class  : L1TopoRoI                                   *
 *                                                      *
 * Description: This class holds info about a L1 RoI.   *
 * Derived classes can be implemented for L1 muon and   *
 * L1Calo RoIs. The purpose is to provide a common base *
 * class for RoIs that L1Topo algs can operate on.      *
 *                                                      *
 * Authors (alphabetical):                              *
 *      Christian Ohm (CERN)                            *
 ********************************************************/

#ifndef L1TOPO_L1TOPOROI_H
#define L1TOPO_L1TOPOROI_H

#include <iostream>

namespace L1Topo {
  
  // define generic L1Topo RoI base class
  class L1TopoRoI {
    
   public:
    
    L1TopoRoI() : m_eta(0), m_phi(0), m_pT(0), m_type(0) {};
    L1TopoRoI(float eta, float phi, float pT, int type)
    : m_eta(eta), m_phi(phi), m_pT(pT), m_type(type) {};
    virtual ~L1TopoRoI() {};
    
    // member functions
    float eta() const {return m_eta;};
    float phi() const {return m_phi;};
    int   pT() const {return m_pT;};
    int type() const {return m_type;};
    void  setEta(float eta) {m_eta = eta;};
    void  setPhi(float phi) {m_phi = phi;};
    void  setpT(float pT) {m_pT = pT;};
    void  setType(int type) {m_type = type;};
    virtual void  print() {
      std::cout << "  RoI: eta=" << m_eta << ", phi=" << m_phi
                << ", pT=" << m_pT << ", type=" << m_type << std::endl;
    };
    
    // for sorting (in case there are more than we can send in a MIOCT)
    bool operator < (const L1TopoRoI& roi) const {
      return (m_pT < roi.pT());
    }
    
   protected:
    
    float m_eta;
    float m_phi;
    int   m_pT;
    int m_type;
    
  }; // class L1TopoRoI
  
} // namespace L1Topo

#endif
