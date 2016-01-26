
```
setupATLAS
rcSetup Base,2.3.40
svn co svn+ssh://svn.cern.ch/reps/atlasgroups/Trigger/L1CTUpgrade/MuCTPIPhase0Upgrade/tags/MuCTPIPhase0Upgrade-00-01-01 MuCTPIPhase0Upgrade
rc find_packages
rc compile 
python skimL1.py `cat input_filelist.txt` 2>&1 | tee skimL1.log
```

The input file is from
```
rucio list-files data15_hi.00287924.physics_L1Calo.merge.AOD.f657_m1533/ | egrep 'lb022'
rucio download data15_hi:data15_hi.00287924.physics_L1Calo.merge.AOD.f657_m1533._lb0220-lb0222._0001.1
```

davide.gerbaudo@gmail.com

Jan 2016
