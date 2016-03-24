#ifndef L1TOPORDO_MUONTOB_H
#define L1TOPORDO_MUONTOB_H

#include <iostream>
#include <cstdint>

namespace L1Topo {

  /** @short Represents the Muon TOB word of the L1Topo ROI data, with decode and encoder
   *
   * Used to decode/encode the 32-bit word that contains the L1
   * muon trigger input to the L1Topo processor.
   * The 32 bits are specified in
   * https://twiki.cern.ch/twiki/pub/Atlas/L1CaloUpgrade/ROD_data_format_v1.0.9.xlsx
   * C1 and C2 are the two candidates
   * | 4b     | 1b   | 3b   | 1b      | 7b       | 3b     | 3b     | 2b    | 3b     | 3b     | 2b    |
   * | 4b     | 1b   | 3b   | 1b      | 7b       |          8b             |          8b             |
   * | header | side | part | delayed | reserved | eta C2 | phi C2 | pt C2 | eta C1 | phi C1 | pt C1 |
   *  31--------------------------------------bit numbers-------------------------------------------00
   * | 3322     2      222    2         2221111    111      111      00      000      000       00
   * | 1098     7      654    3         2109876    543      210      98      765      432       10

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
    //! whether this octant has more than two candidates
    bool overflow() const;
    //! whether these are delayed muons
    bool delayed() const;
    //! side: A=0, C=1
    bool side() const;
   protected:
    //! method used by constructor to decode word
    void decode();
  private:
    /// \todo DG-2016-03-23 do we need to store on the side also the other 16 bits?
    uint32_t m_candidate1_bits;
    uint32_t m_candidate2_bits;
    uint32_t m_word;
    const uint32_t m_bits_per_candidate = 8;
    const uint32_t m_side_bit_number = 27;
    const uint32_t m_delayed_bit_number = 23;
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

