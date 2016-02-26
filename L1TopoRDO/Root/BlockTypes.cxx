#include "L1TopoRDO/BlockTypes.h"

namespace L1Topo{

  L1Topo::BlockTypes blockType(const uint32_t word){
    return static_cast<BlockTypes>(word>>28 & 0x0f);
  }

} // namespace L1Topo
