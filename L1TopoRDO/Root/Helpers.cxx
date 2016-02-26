#include "L1TopoRDO/Helpers.h"

#include <vector>
#include <bitset>
#include <algorithm> 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdint>
//#include <boost/io/ios_state.hpp>
// #include "L1TopoRDO/L1TopoRDO.h"
// #include "L1TopoRDO/L1TopoRDOCollection.h"
#include "L1TopoRDO/Helpers.h"
#include "L1TopoRDO/Header.h"
#include "L1TopoRDO/Fibre.h"
#include "L1TopoRDO/Status.h"
#include "L1TopoRDO/L1TopoTOB.h"
#include "L1TopoRDO/BlockTypes.h"


namespace L1Topo{
  
  std::string formatHex8(uint32_t word){
    std::ostringstream s;
    s << std::showbase << std::hex << std::internal << std::setfill ('0') << std::setw(10) << word << std::dec << std::noshowbase;
    return s.str();
  }

  std::string formatHex4(uint32_t word){
    std::ostringstream s;
    s << std::showbase << std::hex << std::internal << std::setfill ('0') << std::setw( 6 ) << word << std::dec << std::noshowbase;
    return s.str();
  }

  unsigned int triggerBitIndex(uint32_t moduleId, L1Topo::L1TopoTOB c){
    uint32_t module = (moduleId >>4) & 0x1;
    uint32_t index = 64*module + 32*c.fpga() + 16*c.clock() + 8*c.index();
    return index;
  }

  // convenience wrapper for use with L1TopoResult - creates dependency loop with TrigT1Result
  /*
  std::pair< std::bitset<128>,std::bitset<128> > getDecisionAndOverflowBits(const std::vector<L1TopoResult>& res){
    L1TopoRDOCollection col;
    for (auto & r : res){
      col.push_back(new L1TopoRDO(r.rdo()));
    }
    return getDecisionAndOverflowBits(col);
  }
  */


} // namespace L1Topo
