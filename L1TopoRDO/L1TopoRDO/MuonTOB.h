#ifndef L1TOPORDO_MUONTOB_H
#define L1TOPORDO_MUONTOB_H

#include <iostream>
#include <cstdint>

namespace L1Topo {

  /** @short Represents the Muon TOB word of the L1Topo ROI data, with decode and encoder
   *
   * Used to decode/encode the 32-bit (16?) word that contains the L1
   * muon trigger input to the L1Topo processor.
   * This class is a thin wrapper around MIOCTPhase0TopoDataDecoder
   * that provides the same behaviour as the existing decoders by
   * Simon.
   */
  class MuonTOB {
  public:
    //! Construct with word and decode contents
    MuonTOB(const int32_t word);
    /* //! \todo provide access to decoded minirois */
    /* uint32_t minirois() const; */
   protected:
    //! method used by constructor to decode word
    void decode();
  private:
    uint32_t m_word;
  };

} // namespace L1Topo

//! Helper for printing
std::ostream& operator<<(std::ostream&, const L1Topo::MuonTOB&);

// these might be needed?
/* //! Comparison operators, based on word */
/* bool operator==(const L1Topo::MuonTOB&, const L1Topo::MuonTOB&); */
/* //! Comparison operators, based on word */
/* bool operator!=(const L1Topo::MuonTOB&, const L1Topo::MuonTOB&); */
/* //! Comparison operators, based on word */
/* bool operator< (const L1Topo::MuonTOB&, const L1Topo::MuonTOB&); */
/* //! Comparison operators, based on word */
/* bool operator> (const L1Topo::MuonTOB&, const L1Topo::MuonTOB&); */

#endif // L1TOPORDO_MUONTOB_H

