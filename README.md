### Datadriven scripts for HHWWgg background estimation

## step1: create input files (DiphotonJetBox, Data, GJets)
```
python mk_pp_gjet.py
```
## step2: use GJets fakephoton pdf to reweight the data sideband events(QCD)
```
python makeDataDrivenQCD.py
```
## step3: get the Scale factors for pp and dd background
```    
root -L deriveSFs.C
```
## step4: apply sfs to pp and dd 
```
python addSFs.py
```