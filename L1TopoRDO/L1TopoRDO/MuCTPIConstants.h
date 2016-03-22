// Dear emacs, this is -*- c++ -*-
/*********************************************************
 * MuCTPIConstants.h                                     *
 *                                                       *
 * Description: This header file defines constants       *
 * related to the segmentation of the MuCTPI system      *
 * and an enum defining the various proposed encoding    *
 * schemes supported by the encoder and decoder classes. *
 *                                                       *
 * Authors (alphabetical):                               *
 *      Christian Ohm (CERN)                             *
 *********************************************************/

#ifndef L1TOPO_MUCTPICONSTANTS_H
#define L1TOPO_MUCTPICONSTANTS_H

namespace L1Topo {

  // enum specifying chosen encoding scheme
  enum MioctEncodingScheme {
    one4x2hitmappT16bit,    // One 2x4 hitmap, 2 bits per bin for pT
    one4x2hitmapmult16bit,  // One 2x4 hitmap, 2 bits per bin for multiplicity
    //two4x2hitmappT16bit,    // Two 2x4 hitmaps, 1 bit filled/empty per bin (worse than one2x4hitmapmult16bit!)
    two2x2hitmappT16bit,  // Two 2x2 hitmaps, 2 bits per bin for pT
    one2x2hitmappT8bit,     // One 2x2 hitmap, 2 bits per bin for pT
    one2x2hitmapmult8bit,   // One 2x2 hitmap, 2 bits per bin for multiplicity
    
    miniroi1cand16bit,      // 6 bits eta, 5 bits phi, 3 bits pT, 1 bit charge (alternative: RoI number instead of eta/phi bits?)
    miniroi2cands16bit,     // (3 bits eta, 2 bits phi, 2 bits pT)/cand, 2 last bits:  00: no cand seen
    //                                                         01: only 1st cand seen
    //                                                         10: two cands seen
    //                                                         11: Overflow - more than 2 cands
    miniroi2cands16bit2,     // (3 bits eta, 3 bits phi, 1 bit pT)/cand, 2 last bits:  00: no cand seen
    miniroi2cands16bitFinal1, // (3 bits eta, 3 bits phi, 2 bits pT)/cand, eta bits = 111 means no candidate, pT bits = 11 in second candidate signals overflow (pT bits of first still meaningful and gives highest pT threshold passed)
    miniroi2cands16bitFinal2, // like miniroi2cands16bitFinal1 but with 2 eta bins for the forward region and 5 for barrel+endcap, rather than 1 and 6.
    //                                                         01: only 1st cand seen
    //                                                         10: two cands seen
    //                                                         11: Overflow - more than 2 cands
    miniroi1cand8bit,       // 3 bits eta, 2 bits phi, 2 bits pT, last bit - 0: no cand seen, 1: at least one cand seen
    miniroi2cands8bit,      // (2 bits eta, 1 bit phi)/cand, 2 last bits:              00: no cand seen
    //                                                         01: only 1st cand seen
    //                                                         10: two cands seen
    //                                                         11: Overflow - more than 2 cands
    miniroi3cands16bit      // (2 bits eta, 2 bits phi)/cand, (2 bits???) 2 last bits: 00: no cand seen
    //                                                         01: only 1st cand seen
    //                                                         10: 1st and 2nd cands seen
    //                                                         11: all three cands seen
  };
  
  // enum for output levels
  enum OutputLevel {
    ERROR = 0,
    WARNING = 1,
    INFO = 2,
    DEBUG = 3
  };
  
  // constants
  const float mioctMinEtaBarrel = 0.0;
  const float mioctMaxEtaBarrel = 1.06;
  const float mioctMinEtaEndcap = 1.06;
  const float mioctMaxEtaEndcap = 1.965;
  const float mioctMinEtaForward = 1.94;
  const float mioctMaxEtaForward = 2.46;
  const float mioctDeltaPhiBarrel = 0.7;
  const float mioctDeltaPhiEndcap = 0.77;
  const float mioctDeltaPhiForward = 0.74;
  const float mioctDeltaPhiMax = mioctDeltaPhiEndcap;
  const float mioct0MinPhiBarrel = -0.35;
  const float mioct0MinPhiEndcap = -0.38;
  const float mioct0MinPhiForward = -0.24;
  const float mioct0MinPhi = mioct0MinPhiEndcap;
}

#endif
