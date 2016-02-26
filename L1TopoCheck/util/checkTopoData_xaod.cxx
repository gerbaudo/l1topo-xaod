#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/tools/ReturnCheck.h"

#include "xAODRootAccess/TEvent.h"
#include "xAODEventInfo/EventInfo.h"

#include "xAODTrigL1Calo/L1TopoRawData.h"
#include "xAODTrigL1Calo/L1TopoRawDataContainer.h"
#include "L1TopoRDO/BlockTypes.h"
#include "L1TopoRDO/Header.h"
#include "L1TopoRDO/Helpers.h"
#include "L1TopoRDO/L1TopoTOB.h"

#include "TChain.h"
#include "TError.h"
#include "TFile.h"
#include "TH1F.h"

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

        // auto l1rd_it = l1toporawdatas->begin();
        // auto l1rd_it_end = l1toporawdatas->end();

        for(auto &l1topo : *l1toporawdatas) {
            cout<<"l1topo.sourceID "<<std::hex<<l1topo->sourceID() <<std::dec<<endl;
            for(auto word : l1topo->dataWords()){
                if (L1Topo::BlockTypes::HEADER==L1Topo::blockType(word)){
                    auto header = L1Topo::Header(word);
                    cout<<header<<endl;
                } else if(L1Topo::BlockTypes::L1TOPO_TOB==L1Topo::blockType(word)) {
                    auto tob = L1Topo::L1TopoTOB(word);
                    cout<<tob<<endl;
                    cout<<"dataword "<<std::hex<<word<<std::dec<<endl;
                    // collect trigger and overflow bits in bitsets
                    for (unsigned int i=0; i<8; ++i){  // 8 bits/word?
                        unsigned int index = L1Topo::triggerBitIndexNew(l1topo->sourceID(), tob, i); // implement this
                        const bool trigger_bit = (tob.trigger_bits()>>i)&1;
                        const bool overflow_bit = (tob.overflow_bits()>>i)&1;
                    }
                    // check fibers errors (see fibers in monitoring) -- status flag
                    // run monitoring with debug on and get the output to compare against


                }
            } // for(word)
        } // for(l1topo)

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
