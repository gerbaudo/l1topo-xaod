#!/bin/env python
# basic check of L1Topo data in file

# Used in ATN tests:
topoFile='/afs/cern.ch/atlas/project/trigger/pesa-sw/validation/atn-test/data15_13TeV.00278748.physics_Main.daq.RAW._lb0103._SFO-1._0001_selected.data'

#topoFile='/afs/cern.ch/work/s/sgeorge/atlaspublic/L1TopoCnv/data_test.00253868.physics_Main/data_test.00253868.physics_Main.daq.RAW._lb0003._SFO-10._0001.data'
#topoFile='/afs/cern.ch/work/s/sgeorge/atlaspublic/L1TopoCnv/data_test.00253206.physics_CosmicCalo/data_test.00253206.physics_CosmicCalo.daq.RAW._lb0002._SFO-10._0001.data'

# root://eosatlas//eos/atlas/atlastier0/rucio/data15_13TeV/physics_EnhancedBias/00266904/data15_13TeV.00266904.physics_EnhancedBias.merge.RAW/data15_13TeV.00266904.physics_EnhancedBias.merge.RAW._lb0410._SFO-1._0001.1

import sys
import eformat
import libpyevent_storage as EventStorage

files=[]
#print len(sys.argv),sys.argv
if len(sys.argv)<=1:
    files=[topoFile]
else:
    files=sys.argv[1:]
maxevents=500

streamT = eformat.istream(files)
print 'The total number of events in files is %d' % (len(streamT))
if len(streamT)>maxevents:
    print 'Max events=',maxevents
print 'Files:', files
fragmentSizes={}
tt_mismatches=0
counter=0
sid_counts={}
for eventT in streamT:       
    counter+=1
    if counter>maxevents:
        break
    print "\n>>> counter, event global_id, LB, bc_id, lvl1_id, l1tt:", counter, eventT.global_id(), eventT.lumi_block(), eventT.bc_id(), eventT.lvl1_id(), eventT.lvl1_trigger_type()
    # CTP
    ctpSIDs= [r.source_id() for r in eventT.children() if r.source_id().subdetector_id()==eformat.helper.SubDetector.TDAQ_CTP]
    for sid in ctpSIDs:
        hsid=hex(sid.code())
        print "CTP sourceID", hsid        
        if hsid in sid_counts:
            sid_counts[hsid]+=1
        else:
            sid_counts[hsid]=1
        rob=eventT.find(sid)
        if len(rob.rod_data())==0:
            print hsid,"  --warning no data"
        try:
            rob.check_rob()
            rob.check_rod()
        except RuntimeError as e:
            print hsid,'  --warning check_ro[bd] Exception!',e, e.args
        print counter, eventT.global_id(), hsid,'fragment_size, rod_data size: ', rob.fragment_size_word(), len(rob.rod_data())
    #continue
    # Topo
    topoSIDs= [r.source_id() for r in eventT.children() if r.source_id().subdetector_id()==eformat.helper.SubDetector.TDAQ_CALO_TOPO_PROC]
    for sid in topoSIDs:
        hsid=hex(sid.code())
        print "L1Topo sourceID", hsid
        if hsid in sid_counts:
            sid_counts[hsid]+=1
        else:
            sid_counts[hsid]=1
        #continue
        rob=eventT.find(sid)
        if len(rob.rod_data())==0:
            print hsid,"  --warning no data"
        try:
            rob.check_rob()
            rob.check_rod()
        except RuntimeError as e:
            print hsid,'  --warning check_ro[bd] Exception!',e, e.args
        print hsid,'rob,rod sids:',hex(rob.rob_source_id().code()), hex(rob.rod_source_id().code())
        print hsid,'bcid, l1id, l1tt:', rob.rod_bc_id(), rob.rod_lvl1_id(), rob.rod_lvl1_trigger_type()
        if rob.rod_lvl1_trigger_type() != eventT.lvl1_trigger_type():
            print hsid, eventT.global_id(), "  --warning: trigger type mismatch: rod=%s event=%s" % (rob.rod_lvl1_trigger_type(),eventT.lvl1_trigger_type())
            tt_mismatches+=1
        if rob.rod_bc_id() != eventT.bc_id():
            print hsid,"  --warning: bcid mismatch: rod=%s event=%s" % (rob.rod_bc_id(),eventT.bc_id())
        print hsid,'check_rob, check_rod, checksum: ', rob.check_rob_noex(), rob.check_rod_noex(), rob.checksum()
        print hsid,'rob, rod status: ', rob.status(), rob.rod_status()
        #print hsid,'rob,rod header: ', rob.header(),rob.rod_header(), 'payload size, payload: ', rob.payload_size_word(),rob.payload()
        print hsid,'rod_data: ',rob.rod_data()
        print counter, eventT.global_id(), hsid,'fragment_size, rod_data size: ', rob.fragment_size_word(), len(rob.rod_data())
        if sid.code()==0x910000 and (rob.fragment_size_word()<150 or rob.fragment_size_word()>220):
            print hsid," --warning strange fragment size! event=%s size=%s" % ( eventT.global_id(), rob.fragment_size_word() )
        if not sid.code() in fragmentSizes:
            fragmentSizes[sid.code()]=[]
        fragmentSizes[sid.code()].append(rob.fragment_size_word())
        print hsid, 'rod_data: event# ', counter, ' trigger bits set ', len([x for x in rob.rod_data() if (x>>28 == 0x8) and (x & 0xff)]), ' overflow bits set ', len([x for x in rob.rod_data() if (x>>28 == 0x8) and (x>>8 & 0xff)])


print "tt mismatches=",tt_mismatches
#print "L1Topo fragment sizes: "
#for k,v in fragmentSizes.items():
#    print "  %s average %d min %d max %d " % (hex(k), sum(v)/len(v),min(v),max(v))

print "Seen these source IDs with counts given:\n ", sid_counts


def dump_events_with_long_roi_frag():
    import sys
    import eformat
    import libpyevent_storage as EventStorage
    fname='data15_cos.00262656.debug_DcmL1IdMismatchError.daq.RAW._lb0000._SFO-1._0001.data'
    streamT = eformat.istream(fname)
    sid=eformat.helper.SourceIdentifier(0x910081)
    events=[ (eventT.header(),eventT.payload()) for eventT in streamT if len((eventT.find(sid)).rod_data())==14 ]
    
    f=open('evdump.txt','w')
    f.write('header:\n')
    for w in events[0][0]:
        f.write('0x%08x\n' % w)
        
    f.write('payload:\n')
    for w in events[0][1]:
        f.write('0x%08x\n' % w)
            
    f.close()
