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
   *  31--------------------------------------bit numbers----------------------------------------00
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
    //! ATLAS convention: A side --> z>0 or eta>0
    uint32_t signed_side() const;
    //! octant: 0--7
    uint32_t octant() const;
    uint32_t word() const;
    uint32_t candidate1() const;
    uint32_t candidate2() const;
    bool has_candidate1() const;
    bool has_candidate2() const;
  public:
    static const uint32_t m_side_bit_number = 27;
    static const uint32_t m_delayed_bit_number = 23;
    static const uint32_t m_octant_lower_bit_number = 24;
    static const uint32_t m_number_of_bits_per_candidate = 8;
    static const uint32_t m_number_of_eta_bits = 3;
    static const uint32_t m_number_of_phi_bits = 3;
    static const uint32_t m_number_of_pt_bits = 2;
    //! this eta value is reserved to indicate no candidate for these 8 bits
    static const uint32_t m_eta_means_no_candidates = 0b111;
    //! decode eta from lowest 8 bits
    static uint32_t decode_eta_from_8_bits(const uint32_t &bits);
    //! decode phi from lowest 8 bits
    static uint32_t decode_phi_from_8_bits(const uint32_t &bits);
    //! decode pt from lowest 8 bits
    static uint32_t decode_pt_from_8_bits(const uint32_t &bits);
   protected:
    //! method used by constructor to decode word
    void decode();
  private:
    /// \todo DG-2016-03-23 do we need to store on the side also the other 16 bits?
    uint32_t m_candidate1_bits;
    uint32_t m_candidate2_bits;
    uint32_t m_word;
  };

  /** @short One of the two possible muon candidates from the muon TOB word
   *
   * Helper class, used to convert to geometrical values from int to double and print.
   */
  class MuonCandidate {
  public:
    MuonCandidate(const int32_t full_word, const int32_t candidate_word);
    double m_eta;
    double m_phi;
    int m_pt;
  };

} // namespace L1Topo

//! Helper for printing
std::ostream& operator<<(std::ostream&, const L1Topo::MuonTOB&);

std::ostream& operator<<(std::ostream&, const L1Topo::MuonCandidate&);

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

