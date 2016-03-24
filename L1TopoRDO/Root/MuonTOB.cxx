#include "L1TopoRDO/MuonTOB.h"
#include "L1TopoRDO/Helpers.h"
#include "L1TopoRDO/BlockTypes.h"
#include "L1TopoRDO/MIOCTPhase0TopoRoI.h"
#include "L1TopoRDO/MuCTPIConstants.h"

#include "TMath.h"
#include "TVector2.h"

#include <iostream>
#include <iomanip>
#include <bitset>
#include <vector>
#include <cmath>

namespace L1Topo {

using std::vector;
// DG-2016-03-23
// this is a temporary declaration of the functions from Christian's encoder/decoder
std::vector<MIOCTPhase0TopoRoI> decode(const std::vector<std::vector<bool> >& dataWords);
std::vector<MIOCTPhase0TopoRoI> decode_miniroi2cands16bitFinal(const std::vector<bool>& dataWord,
                                                               unsigned int mioct, int nFwdEtaBins);
MIOCTPhase0TopoRoI decodeRoIFromBits(const std::vector<bool>& bits, const int nEtaBits,
                                     const int nPhiBits, const int nPtBits, const bool decodeCharge,
                                     const int mioct, const int nFwdEtaBins);

int decodeIntFromBits(const std::vector<bool>& bits);
std::vector<bool> uint32_t2bits(const uint32_t &v);
//  end of temporary declarations

MuonTOB::MuonTOB(const int32_t word)
    : m_word(word){

    this->decode();
}

void MuonTOB::decode()
{
    const uint32_t eight_1 = 0b11111111;
    m_candidate1_bits = m_word & eight_1;
    m_candidate2_bits = (m_word >> m_number_of_bits_per_candidate) & eight_1;

    const int nEtaBits = m_number_of_eta_bits;
    const int nPhiBits = m_number_of_phi_bits;
    const int nPtBits = m_number_of_pt_bits;
    const bool decodeCharge = false;
    const int mioct = octant(); // DG-2016-03-24 check with Olya
    const int nFwdEtaBins = 1; // from 'miniroi2cands16bitFinal1' encoding
    vector<bool> bitsc1 = uint32_t2bits(m_candidate1_bits); bitsc1.resize(8);
    vector<bool> bitsc2 = uint32_t2bits(m_candidate2_bits); bitsc2.resize(8);
    MIOCTPhase0TopoRoI candidate1 = decodeRoIFromBits(bitsc1, nEtaBits, nPhiBits, nPtBits,
                                                      decodeCharge, mioct, nFwdEtaBins);
    MIOCTPhase0TopoRoI candidate2 = decodeRoIFromBits(bitsc2, nEtaBits, nPhiBits, nPtBits,
                                                      decodeCharge, mioct, nFwdEtaBins);
    return;
}

bool MuonTOB::overflow() const
{
    // both pt bits of C2 set to 1
    // DG-2016-03-24 ask Olya whether it's C1 or C2; the order seems to be swapped between the xml and the wiki
    const uint32_t two_1 = 0b11;
    return m_candidate2_bits & two_1;
}
bool MuonTOB::delayed() const
{
    return (m_word >> m_delayed_bit_number) & 1;
}
bool MuonTOB::side() const
{
    return (m_word >> m_side_bit_number) & 1;
}

uint32_t MuonTOB::signed_side() const
{
    return side() ? -1 : +1;
}

uint32_t MuonTOB::octant() const
{
    const uint32_t three_1 = 0b111;
    return (m_word >> m_octant_lower_bit_number) & three_1;
}

uint32_t MuonTOB::word() const
{
    return m_word;
}

uint32_t MuonTOB::candidate1() const
{
    return m_candidate1_bits;
}

uint32_t MuonTOB::candidate2() const
{
    return m_candidate2_bits;
}

bool MuonTOB::has_candidate1() const
{
    return decode_eta_from_8_bits(m_candidate1_bits) != m_eta_means_no_candidates;
}

bool MuonTOB::has_candidate2() const
{
    return decode_eta_from_8_bits(m_candidate2_bits) != m_eta_means_no_candidates;
}

uint32_t MuonTOB::decode_eta_from_8_bits(const uint32_t &bits)
{
    const uint32_t three_1 = 0b111;
    return (bits >> 5 ) & three_1;
}

uint32_t MuonTOB::decode_phi_from_8_bits(const uint32_t &bits)
{
    const uint32_t three_1 = 0b111;
    return (bits >> 2 ) & three_1;
}

uint32_t MuonTOB::decode_pt_from_8_bits(const uint32_t &bits)
{
    const uint32_t two_1 = 0b11;
    return bits & two_1;
}

MuonCandidate::MuonCandidate(const int32_t full_word, const int32_t candidate_word):
    m_eta(0.0),
    m_phi(0.0),
    m_pt(0)
{
    MuonTOB tob(full_word);
    int side = tob.signed_side();
    int octant = tob.octant();
    int eta = MuonTOB::decode_eta_from_8_bits(candidate_word);
    int phi = MuonTOB::decode_phi_from_8_bits(candidate_word);
    // convert the bin values to coordinates, corresponding to the center values of the bins
    bool isForward = eta==6;
    const double phi_octant_width = TMath::Pi()/4.0;
    const double phi_bin_width = mioctDeltaPhiMax / 8.0; // 8 = pow(2, nPhiBits=3)
    m_eta = (isForward ?
             (mioctMinEtaForward+mioctMaxEtaForward)/2 : // fw candidate in bin 6
             mioctMinEtaBarrel+(eta+0.5)*(mioctMaxEtaEndcap-mioctMinEtaBarrel)/6); // barrel or endcap
    m_eta *= side;
    // DG-2016-03-24 I dont understand the line below; ask Olya or Christian.
    // eta *= ((octant > 7) ? 1 : -1); // set sign of eta based on implicit hemisphere from MIOCT number
    m_phi = TVector2::Phi_mpi_pi(mioct0MinPhi + phi_octant_width*octant + (phi+0.5)*phi_bin_width);
    if(isForward) // if this is a forward candidate, compensate for difference in MIOCT phi coverage
        m_phi += 0.14;
    m_pt = MuonTOB::decode_pt_from_8_bits(candidate_word);
}

/**
 * Does the decoding of one MIOCT for the miniroi2cands16bitFinal[1,2] schemes
 *
 * @param dataword The bits containing the L1 muon RoI variables
 * @param mioct The number of the MIOCT module the candidate was found in
 * @param nFwdEtaBins The number of eta bins used for the forward region
 */
std::vector<MIOCTPhase0TopoRoI> decode_miniroi2cands16bitFinal(const std::vector<bool>& dataWord,
                                                               unsigned int mioct, int nFwdEtaBins)
{
    bool verbose = false;
    std::vector<MIOCTPhase0TopoRoI> rois;
    if (dataWord.size() != 16) {
        std::cout << "In decode_miniroi2cands16bitFinal(): ERROR - asked to decode MIOCT data word with "
                  << dataWord.size() << " bits but expected 16 => will exit." << std::endl;
        exit(1);
    }
    const int number_of_candidates_per_word = 2;
    const int number_of_bits_per_candidate = 8;
    const int number_of_eta_bits_per_candidate = 3;
    for (int roi = 0; roi < number_of_candidates_per_word; ++roi) {
        std::vector<bool>::const_iterator first = dataWord.begin()+roi*number_of_bits_per_candidate;
        std::vector<bool>::const_iterator last = dataWord.begin()+(roi+1)*number_of_bits_per_candidate;
        std::vector<bool> bits(first, last);
        if (verbose) {
            std::cout << "Will now decode " << bits.size() << " bits representing an RoI" << std::endl;
        }
        const int eta_means_no_candidate = 7; // the eta value 0b111 = 7 is reserved to signal no candidate
        if (decodeIntFromBits(std::vector<bool>(first, first + number_of_eta_bits_per_candidate)) == eta_means_no_candidate) {
            continue;
        }
        rois.push_back(decodeRoIFromBits(bits, 3, 3, 2, false, mioct, nFwdEtaBins));
        if (verbose) {
            (rois[rois.size()-1]).print();
        }
        // check for overflow in the octant (11 in second candidate's pT bits)
        if (roi == 1) {
            /*
              std::cout << "Will now check for overflow when decoding second RoI - last two bits: "
              << get_binrep_boolvec(std::vector<bool>(last-2, last)) << " => "
              << decodeIntFromBits(std::vector<bool>(last-2, last)) << std::endl;
            */
            if (decodeIntFromBits(std::vector<bool>(last-2, last)) == 3) {
                (rois[rois.size()-1]).setThrNumber(0); // change to mock pT value since unknown
                (rois[rois.size()-1]).setpT(0); // change to mock pT value since unknown
                //rois.push_back(new MIOCTPhase0TopoRoI(0, 0, 0, 0, false, 0)); // for debugging: dummy RoI to show overflow
                /*std::cout << "Added dummy RoI since there was overflow: pT = "
                  << get_binrep_boolvec(std::vector<bool>(first+6, first+8)) << std::endl;*/
            }
        }
    }
    // if (verbose) {
    //     std::cout << "Decoded dataword " << get_binrep_boolvec(dataWord)
    //               << ", found " << rois.size() << " RoI(s)" << std::endl;
    // }
    return rois;
}


MIOCTPhase0TopoRoI decodeRoIFromBits(const std::vector<bool>& bits, const int nEtaBits,
                                                const int nPhiBits, const int nPtBits, const bool decodeCharge,
                                                const int mioct, const int nFwdEtaBins)
{
    // extract the bin values from the bits
    int offset = 0;
    std::vector<bool>::const_iterator first = bits.begin();
    int eta = decodeIntFromBits(std::vector<bool>(first, first + nEtaBits));
    offset += nEtaBits;
    int phi = decodeIntFromBits(std::vector<bool>(first + offset, first + offset + nPhiBits));
    offset += nPhiBits;
    int pT = 0;
    if (nPtBits > 0) {
      pT = decodeIntFromBits(std::vector<bool>(first + offset, first + offset + nPtBits))+1;
      offset += nPtBits;
    }
    int charge = 100; // default value, should be interpreted as N/A
    if (decodeCharge) { // WARNING: we can't know if a candidate is from RPC or TGC endcap! TODO: fix?
      charge = ((bits[(offset += 1)]) ? 1 : -1);
    }
    bool isFwd = false;
    // convert the bin values to coordinates, corresponding to the center values of the bins
    float phiValue = mioct0MinPhi+TMath::Pi()/4*mioct+(phi+0.5)*(mioctDeltaPhiMax)/pow(2, nPhiBits);
    // convert back to ]-pi, pi] interval
    while (phiValue > TMath::Pi())
      phiValue -= 2*TMath::Pi();
    float etaValue = 0.0;
    if (nEtaBits == 3) {
      if (nFwdEtaBins < 0) {
        if (eta > 5) { // forward candidate, two bins (6 and 7)
          etaValue = mioctMinEtaForward+(eta-6+0.5)*(mioctMaxEtaForward-mioctMinEtaForward)/2;
          isFwd = true;
        }
        else { // barrel or endcap
          etaValue = mioctMinEtaBarrel+(eta+0.5)*(mioctMaxEtaEndcap-mioctMinEtaBarrel)/6;
        }
      }
      // below two cases valid for miniroi2cands16bitFinal[1,2] encoding schemes
      else if (nFwdEtaBins == 1) {
        if (eta == 6) { // forward candidate in bin 6
          etaValue = (mioctMinEtaForward+mioctMaxEtaForward)/2;
          isFwd = true;
        }
        else { // barrel or endcap
          etaValue = mioctMinEtaBarrel+(eta+0.5)*(mioctMaxEtaEndcap-mioctMinEtaBarrel)/6;
        }
      }
      else if (nFwdEtaBins == 2) {
        if (eta > 4) { // forward candidate in bins 5 and 6
          etaValue = mioctMinEtaForward+(eta-5+0.5)*(mioctMaxEtaForward-mioctMinEtaForward)/2;
          isFwd = true;
        }
        else { // barrel or endcap in bins 0-4
          etaValue = mioctMinEtaBarrel+(eta+0.5)*(mioctMaxEtaEndcap-mioctMinEtaBarrel)/5;
        }
      }
    } // end if nEtaBits == 3
    else if (nEtaBits == 2) {
      if (eta > 2) { // forward candidate, one bin (3)
        etaValue = (mioctMinEtaForward+mioctMaxEtaForward)/2;
        isFwd = true;
      }
      else { // barrel or endcap
        etaValue = mioctMinEtaBarrel+(eta+0.5)*(mioctMaxEtaEndcap-mioctMinEtaBarrel)/3;
      }
    }
    else {
      std::cout << "In decodeRoIFromBits(): ERROR - unsupported number of eta bits ("
                << nEtaBits << "), cannot proceed" << std::endl;
      exit(1);
    }
    etaValue *= ((mioct > 7) ? 1 : -1); // set sign of eta based on implicit hemisphere from MIOCT number
    if (isFwd) // if this is a forward candidate, compensate for difference in MIOCT phi coverage
      phiValue += 0.14;
    return MIOCTPhase0TopoRoI(etaValue, phiValue, pT, charge, isFwd, -1); // dummy sector address -1
}

/**
 * Decodes a vector of bools into an integer
 *
 * @param bits The vector of bools containing the encoded integer
 */
int decodeIntFromBits(const std::vector<bool>& bits) {
    int x = 0;
    for (unsigned int bit = 0; bit < bits.size(); ++bit) {
        x = (x << 1) | bits[bit];
    }
    return x;
}

std::vector<bool> uint32_t2bits(const uint32_t &v) {
    std::bitset<32> bits(v);
    vector<bool> vbits(bits.size());
    for(size_t i=0; i<bits.size(); ++i)
        vbits[i] = bits[i];
    return vbits;
}

} // namespace L1Topo

std::ostream& operator<<(std::ostream& os, const L1Topo::MuonTOB& c)
{
    os << "     MuonTOB: "
       << " side" << c.side()
       << " octant "<< c.octant()
       << " delayed "<< c.delayed()
       << " candidate 1 "<< std::bitset<8>(c.candidate1());
    if(c.has_candidate1())
        os <<"("<<L1Topo::MuonCandidate(c.word(), c.candidate1())<<")";
    else
        os <<"(--)";
    os << " candidate 2 "<< std::bitset<8>(c.candidate2());
    if(c.has_candidate2())
        os <<"("<<L1Topo::MuonCandidate(c.word(), c.candidate2())<<")";
    else
        os <<"(--)";
    os << std::dec << "\n";
    return os;
}

std::ostream& operator<<(std::ostream& os, const L1Topo::MuonCandidate& c)
{
  os << "     MuonCandidate: "
     << " pt "<<c.m_pt
     << " eta "<<c.m_eta
     << " phi "<<c.m_phi;
  return os;
}


