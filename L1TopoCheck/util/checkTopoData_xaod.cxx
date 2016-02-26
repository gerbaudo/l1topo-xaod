#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/tools/ReturnCheck.h"

#include "xAODRootAccess/TEvent.h"
#include "xAODEventInfo/EventInfo.h"

#include "xAODTrigL1Calo/L1TopoRawData.h"
#include "xAODTrigL1Calo/L1TopoRawDataContainer.h"
#include "L1TopoRDO/BlockTypes.h"
#include "L1TopoRDO/Fibre.h"
#include "L1TopoRDO/Header.h"
#include "L1TopoRDO/Helpers.h"
#include "L1TopoRDO/L1TopoTOB.h"

#include "TChain.h"
#include "TError.h"
#include "TFile.h"
#include "TH1F.h"

#include <bitset>
#include <iostream>
#include <ios>

using namespace std;

int main(int argc, char* argv[])
{
    // gROOT->SetBatch(1);
    // gROOT->Macro("$ROOTCOREDIR/scripts/load_packages.C");
    if(not xAOD::Init().isSuccess())
        cout<<"Failed xAOD.Init()"<<endl;

    string default_input_file = "/afs/cern.ch/user/g/gerbaudo/public/tmp/for_simon/xAOD.L1Calo.00287924_lb0100.pool.root";
    string input_filename = default_input_file;

    string treeName = "CollectionTree";
    // TFile *input_file = TFile::Open(input_filename.c_str());

    TChain chain(treeName.c_str());
    chain.Add(input_filename.c_str());
    // kClassAccess (default) will cause a crash because of type changes in some branches (int -> char)

    // TTree *tree = xAOD::MakeTransientTree(&ch, xAOD::TEvent::kBranchAccess);

    xAOD::TEvent event(xAOD::TEvent::kBranchAccess);
    event.readFrom(&chain); //, treeName.c_str());

    size_t nEntries = chain.GetEntries();
    nEntries = 10;
    for (size_t iEntry = 0; iEntry < nEntries; ++iEntry) {
        event.getEntry(iEntry);
        cout<<" entry "<<iEntry<<endl;

        /* event info not available in the test xaod

        const xAOD::EventInfo *info = nullptr;
        event.retrieve(info, "EventInfo");
        printf("\n>>> counter, event global_id, LB, bc_id, lvl1_id, l1tt:"
               "%zu %d %d %d %d %d",
               iEntry,
               0, //info.global_id(),
               info->lumiBlock(),
               info->bcid(),
               info->extendedLevel1ID(),
               info->level1TriggerType());
        */

        const xAOD::L1TopoRawDataContainer *l1toporawdatas = nullptr;
        event.retrieve(l1toporawdatas, "L1TopoRawData");
        // initialise header: beware, this can make a valid-looking
        // header and be misinterpreted; set version 15, BCN -7, which
        // is unlikely:
        L1Topo::Header header(0xf,0,0,0,0,1,0x7);
        static const unsigned int nTopoCTPOutputs = 128; //! Number of CTP outputs
        std::bitset<nTopoCTPOutputs> triggerBits; //! trigger bits sent to CTP
        std::bitset<nTopoCTPOutputs> overflowBits; //! overflow bits corresponding to CTP output
        // first loop: just count the words
        unsigned int nl1topowords = 0;
        for(auto &l1topo : *l1toporawdatas) {
            for(auto &word : l1topo->dataWords()){
                if(L1Topo::BlockTypes::L1TOPO_TOB==L1Topo::blockType(word))
                    nl1topowords++;
            }
        }
        // second loop: just print out the words
        cout<<"L1Topo data words: "<<nl1topowords<<endl;
        for(auto &l1topo : *l1toporawdatas) {
            for(auto &word : l1topo->dataWords()){
                if(L1Topo::BlockTypes::L1TOPO_TOB==L1Topo::blockType(word))
                    cout<<L1Topo::formatHex8(word)<<endl;
            }
        }
        for(auto &l1topo : *l1toporawdatas) {
            cout<<"l1topo.sourceID "<<std::hex<<l1topo->sourceID() <<std::dec<<endl;
            for(auto word : l1topo->dataWords()){
                if (L1Topo::BlockTypes::HEADER==L1Topo::blockType(word)){
                    header = L1Topo::Header(word);
                    cout<<header<<endl;
                } else if(L1Topo::BlockTypes::L1TOPO_TOB==L1Topo::blockType(word)) {
                    auto tob = L1Topo::L1TopoTOB(word);
                    cout<<"dataword 0x"<<std::hex<<word<<std::dec<<endl;
                    cout<<tob<<endl;
                    // collect trigger and overflow bits in bitsets
                    for (unsigned int i=0; i<8; ++i){  // 8 bits/word?
                        unsigned int index = L1Topo::triggerBitIndexNew(l1topo->sourceID(), tob, i);
                        triggerBits[index]  = (tob.trigger_bits()>>i)&1;
                        overflowBits[index] = (tob.overflow_bits()>>i)&1;
                    }
                } else if(L1Topo::BlockTypes::FIBRE==L1Topo::blockType(word)) {
                    // check fibers errors (see fibers in monitoring) -- status flag
                    auto fibreBlock = L1Topo::Fibre(word);
                    auto statuses = fibreBlock.status();
                    for (unsigned int i=0; i<statuses.size(); ++i){
                        if (statuses.at(i)!=0){
                            cout<<" Warning: Fibre status set for fibre "
                                <<i<<" of ROB "<<L1Topo::formatHex8(l1topo->sourceID())<<" header "<<header<<endl;
                        }
                    }
                }
            } // for(word)
        } // for(l1topo)
        cout<<"trigger  bits from RoI Cnv: " <<triggerBits <<endl;
        cout<<"overflow bits from RoI Cnv: " <<overflowBits<<endl;
    } // for iEntry
}
