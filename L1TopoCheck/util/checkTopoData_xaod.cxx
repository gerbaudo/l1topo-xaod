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
#include "L1TopoRDO/Fibre.h"
#include "L1TopoRDO/L1TopoTOB.h"
#include "L1TopoRDO/Status.h"

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
    default_input_file = "/afs/cern.ch/user/b/beallen/public/ForDavide/Data16_cos.00291671.physics_L1Calo.part1.root";
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
    nEntries = 20;
    for (size_t iEntry = 0; iEntry < nEntries; ++iEntry) {
        event.getEntry(iEntry);
        cout<<" entry "<<iEntry<<endl;

        /**/ // event info not available in the test xaod

        const xAOD::EventInfo *info = nullptr;
        event.retrieve(info, "EventInfo");
        printf("\n\n>>> entry %zu, run %d, LB %d, event %llu, bc_id %d, lvl1_id %d, l1tt %d\n",
               iEntry,
               info->runNumber(),
               info->lumiBlock(),
               info->eventNumber(),
               info->bcid(),
               info->extendedLevel1ID(),
               info->level1TriggerType());
        /**/

        const xAOD::L1TopoRawDataContainer *l1toporawdatas = nullptr;
        event.retrieve(l1toporawdatas, "L1TopoRawData");
        int number_of_headers = 0;
        int number_of_l1topotob = 0;
        // initialise header: beware, this can make a valid-looking
        // header and be misinterpreted; set version 15, BCN -7, which
        // is unlikely:
        L1Topo::Header header(0xf,0,0,0,0,1,0x7);
        static const unsigned int nTopoCTPOutputs = 128; //! Number of CTP outputs
        std::bitset<nTopoCTPOutputs> triggerBits; //! trigger bits sent to CTP
        std::bitset<nTopoCTPOutputs> overflowBits; //! overflow bits corresponding to CTP output
        struct {
            bool operator()(const unsigned int v) {
                const std::vector<unsigned int> l1TopoDaqRobIds = {0x00910000, 0x00910010, 0x00910020};
                return std::find(l1TopoDaqRobIds.begin(), l1TopoDaqRobIds.end(), v)!= l1TopoDaqRobIds.end();
            }
        } isL1topoDaqRobId;
        // first loop: just count the words
        unsigned int nl1topowords = 0;
        for(auto &l1topo : *l1toporawdatas) {
            if(not isL1topoDaqRobId(l1topo->sourceID()))
                continue;
            for(auto &word : l1topo->dataWords()){
                if(L1Topo::BlockTypes::L1TOPO_TOB==L1Topo::blockType(word))
                    nl1topowords++;
            }
        }
        // second loop: just print out the words
        cout<<"L1Topo data words: "<<nl1topowords<<endl;
        for(auto &l1topo : *l1toporawdatas) {
            if(not isL1topoDaqRobId(l1topo->sourceID()))
                continue;
            for(auto &word : l1topo->dataWords()){
                if(L1Topo::BlockTypes::L1TOPO_TOB==L1Topo::blockType(word))
                    cout<<L1Topo::formatHex8(word)<<" (sourceId: "<<L1Topo::formatHex8(l1topo->sourceID())<<")"<<endl;
            }
        }
        int previousBunchNumber = 0;
        for(auto &l1topo : *l1toporawdatas) {
            if(not isL1topoDaqRobId(l1topo->sourceID()))
                continue;
            cout<<"l1topo.sourceID "<<std::hex<<l1topo->sourceID() <<std::dec<<endl;
            for(auto word : l1topo->dataWords()){
                const L1Topo::BlockTypes blockType = L1Topo::blockType(word);
                switch(blockType) {
                case L1Topo::BlockTypes::HEADER:     cout<<"HEADER"<<endl;      break;
                case L1Topo::BlockTypes::FIBRE:      cout<<"FIBRE"<<endl;       break;
                case L1Topo::BlockTypes::STATUS:     cout<<"STATUS"<<endl;      break;
                case L1Topo::BlockTypes::EM_TOB:     cout<<"EM_TOB"<<endl;      break;
                case L1Topo::BlockTypes::TAU_TOB:    cout<<"TAU_TOB"<<endl;     break;
                case L1Topo::BlockTypes::MUON_TOB:   cout<<"MUON_TOB"<<endl;    break;
                case L1Topo::BlockTypes::JET1_TOB:   cout<<"JET1_TOB"<<endl;    break;
                case L1Topo::BlockTypes::JET2_TOB:   cout<<"JET2_TOB"<<endl;    break;
                case L1Topo::BlockTypes::ENERGY_TOB: cout<<"ENERGY_TOB"<<endl;  break;
                case L1Topo::BlockTypes::L1TOPO_TOB: cout<<"L1TOPO_TOB"<<endl;  break;
                default:
                    cout<<"unknown blockType for word "<<word<<endl;
                }
                if (L1Topo::BlockTypes::HEADER==blockType){
                    header = L1Topo::Header(word);
                    cout<<"DG L1Topo::header: "<<header<<endl;
                    number_of_headers += 1;
                } else if(L1Topo::BlockTypes::FIBRE==blockType) {
                    auto fibreBlock = L1Topo::Fibre(word);
                    cout<<fibreBlock<<endl;
                } else if(L1Topo::BlockTypes::STATUS==blockType) {
                    auto status = L1Topo::Status(word);
                    cout<<status<<endl;
                } else if(L1Topo::BlockTypes::L1TOPO_TOB==blockType) {
                    if(header.bcn_offset()!=0) {
                        cout<<"skipping L1TopoTOB "<< L1Topo::formatHex8(word)
                            <<" because it comes after header with bcn : "<<header.bcn_offset()
                            <<endl;
                        continue;
                    }
                    auto tob = L1Topo::L1TopoTOB(word);
                    cout<<"dataword 0x"<<std::hex<<word<<std::dec<<endl;
                    cout<<tob<<endl;
                    cout<<"DG L1Topo::L1TopoTOB "<<std::hex<<word<<std::dec<<endl;
                    number_of_l1topotob +=1;
                    // todo
                    // unsigned int index = L1Topo::triggerBitIndexNew(rdo.getSourceID(),tob,i);
                    // check fibers errors (see fibers in monitoring) -- status flag
                    // run monitoring with debug on and get the output to compare against
                } else if(L1Topo::BlockTypes::MUON_TOB==L1Topo::blockType(word)) {
                    cout<<"got moun tob "<<std::hex<<word<<std::dec<<endl;
                }
            } // for(word)
        } // for(l1topo)
        printf("\nin this event we retrieved"
               " %d L1Topo::BlockTypes::HEADER words"
               " and %d L1Topo::BlockTypes::L1TOPO_TOB words\n",
               number_of_headers,
               number_of_l1topotob);

        // const xAOD::ElectronContainer* electrons;
        // CHECK( event.retrieve(electrons, "ElectronCollection") );
        // auto el_it      = electrons->begin();
        // auto el_it_last = electrons->end();
        // unsigned int i = 0;
        // for (; el_it != el_it_last; ++el_it, ++i) {
        //     const xAOD::Electron* el = *el_it;
        //     std::cout << "Electron " << i << std::endl;
        //     std::cout << "xAOD/raw pt = " << el->pt() << std::endl;
        //     Info (APP_NAME,"Electron #%d", i);

        //     const Root::TResult result= myEgCorrections.calculate(*el);

        //     Info( APP_NAME,"===>>> Result 0 position %f ",result.getResult(0));

        // }

    } // for iEntry
}
