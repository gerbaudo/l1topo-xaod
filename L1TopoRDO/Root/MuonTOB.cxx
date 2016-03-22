#include <iostream>
#include <iomanip>
#include <bitset>
#include <vector>
#include "L1TopoRDO/L1TopoTOB.h"
#include "L1TopoRDO/Helpers.h"
#include "L1TopoRDO/BlockTypes.h"
#include "L1TopoRDO/MIOCTPhase0TopoRoI.h"
#include "L1TopoRDO/MuCTPIConstants.h"

namespace L1Topo {


std::vector<MIOCTPhase0TopoRoI> decode(const std::vector<std::vector<bool> >& dataWords);
std::vector<MIOCTPhase0TopoRoI> decode_miniroi2cands16bitFinal(const std::vector<bool>& dataWord,
                                                               unsigned int mioct, int nFwdEtaBins);
MIOCTPhase0TopoRoI decodeRoIFromBits(const std::vector<bool>& bits, const int nEtaBits,
                                     const int nPhiBits, const int nPtBits, const bool decodeCharge,
                                     const int mioct, const int nFwdEtaBins)

L1TopoTOB::L1TopoTOB(const int32_t word)
    : m_word(word){

    this->decode();
}

void L1TopoTOB::decode()
{
    /// \todo how do we go from 32-bit to 16-bit? are we just talking about 2 datawords?
    // anyway for now pretend we know how to do this and convert to vec<bool>
    std::vector<std::vector<bool> > dataWords;
    std::vector<MIOCTPhase0TopoRoI> minirois = decode(dataWords);
}
 


// this is where the decoding happens
std::vector<MIOCTPhase0TopoRoI> decode(const std::vector<std::vector<bool> >& dataWords)
{
    
    std::vector<MIOCTPhase0TopoRoI> rois;
    for (unsigned int mioct = 0; mioct < dataWords.size(); ++mioct) {
        int nFwdEtaBins = 1;
        std::vector<MIOCTPhase0TopoRoI> roisInMioct = decode_miniroi2cands16bitFinal(dataWords[mioct], mioct, nFwdEtaBins);
        if (m_outputLevel >= INFO) {
            std::cout << "  Decoded " << get_binrep_boolvec(dataWords[mioct]) << " => found "
                      << roisInMioct.size() << " RoIs" << std::endl;
        }
        rois.insert(rois.end(), roisInMioct.begin(), roisInMioct.end());
    }
    return rois;
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
        if (m_outputLevel >= DEBUG) {
            std::cout << "Will now decode " << bits.size() << " bits representing an RoI" << std::endl;
        }
        const int eta_means_no_candidate = 7; // the eta value 0b111 = 7 is reserved to signal no candidate
        if (decodeIntFromBits(std::vector<bool>(first, first + number_of_eta_bits_per_candidate)) == eta_means_no_candidate) {
            continue;
        }
        rois.push_back(decodeRoIFromBits(bits, 3, 3, 2, false, mioct, nFwdEtaBins));
        if (m_outputLevel >= DEBUG) {
            (*rois[rois.size()-1]).print();
        }
        // check for overflow in the octant (11 in second candidate's pT bits)
        if (roi == 1) {
            /*
              std::cout << "Will now check for overflow when decoding second RoI - last two bits: "
              << get_binrep_boolvec(std::vector<bool>(last-2, last)) << " => "
              << decodeIntFromBits(std::vector<bool>(last-2, last)) << std::endl;
            */
            if (decodeIntFromBits(std::vector<bool>(last-2, last)) == 3) {
                (*rois[rois.size()-1]).setThrNumber(0); // change to mock pT value since unknown
                (*rois[rois.size()-1]).setpT(0); // change to mock pT value since unknown
                //rois.push_back(new MIOCTPhase0TopoRoI(0, 0, 0, 0, false, 0)); // for debugging: dummy RoI to show overflow
                /*std::cout << "Added dummy RoI since there was overflow: pT = "
                  << get_binrep_boolvec(std::vector<bool>(first+6, first+8)) << std::endl;*/
            }
        }
    }
    if (m_outputLevel >= DEBUG) {
        std::cout << "Decoded dataword " << get_binrep_boolvec(dataWord)
                  << ", found " << rois.size() << " RoI(s)" << std::endl;
    }
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
    //std::cout << "  decoded bits " << get_binrep_boolvec(bits) << " => " << x << std::endl;
    return x;
}

} // namespace L1Topo
