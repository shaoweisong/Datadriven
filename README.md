Datadriven scripts for HHWWgg background estimation

4 steps for datadriven
step1:
    run mk_pp_gjet.py to create input files (DiphotonJetBox, Data, GJets)
step2:
    run makeDataDrivenQCD.py to use GJets fakephoton pdf to reweight the data sideband events(QCD)
step3:
    run deriveSFs.C to get the Scale factors for pp and dd background
step4:
    apply sfs to pp and dd by runing addSFs.py