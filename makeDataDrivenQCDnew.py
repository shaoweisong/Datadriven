#/cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc11-opt/setup.sh
import ROOT
from array import array
import awkward
import sys
sys.path.append("/eos/user/s/shsong/pkgs_condor/parquet_to_root-0.3.0/")
from parquet_to_root import parquet_to_root
import pandas
#run mk_pp_gjet.py firstly to get the Data_2017.root, GJet.root
GJets1=ROOT.TFile.Open("/eos/user/s/shsong/HiggsDNA/GJet/GJet.root")
# GJets1=ROOT.TFile.Open("/eos/user/s/shsong/HiggsDNA/UL16dd/GJet16.root")
GJets1_tree=GJets1.Get("GJet")
cuts="(Diphoton_mass<115 || Diphoton_mass>135)*(minID_genPartFlav!=1)"
# cuts="(Diphoton_mass<115 || Diphoton_mass>135)*(minID_genPartFlav!=1)"
h_minphotonID=ROOT.TH1F("h_minphotonID_gjet","h_minphotonID_gjet",19,-0.9,1)
GJets1_tree.Project("h_minphotonID_gjet","Diphoton_minID_modified","weight_central*"+cuts)
photonIDPDF_fake=ROOT.TF1("photonIDPDF_fake","pol10",-0.9,1.)
h_minphotonID.Fit(photonIDPDF_fake,"R")
fit_result = h_minphotonID.Fit(photonIDPDF_fake, "RS")  
chi2 = photonIDPDF_fake.GetChisquare()
ndf = photonIDPDF_fake.GetNDF()
p_value = ROOT.TMath.Prob(chi2, ndf)
print("chisquare:" ,chi2)
print("freedom degree:", ndf)
print("p_value:", p_value)
c1=ROOT.TCanvas("c1","c1",600,800)
h_minphotonID.Draw("E1")
c1.SaveAs("/eos/user/s/shsong/www/datadriven/fakephoton_pdf17.png")
Data=ROOT.TFile.Open("/eos/user/s/shsong/HiggsDNA/UL17data/Data_2017_cat3.root")
Data_tree=Data.Get("cat3")
nevents=Data_tree.GetEntries()
new_weight=-999
weights=[]

minID=[]
maxID=[]
hasMaxLead=[]
originalminID=[]
print(nevents)
for i in range(0,nevents):
    Data_tree.GetEntry(i)
    # weights.append(1)
    if(Data_tree.LeadPhoton_mvaID_modified < Data_tree.SubleadPhoton_mvaID_modified):
        hasleadIDmin=True
        original_Photon_MVA_min = Data_tree.LeadPhoton_mvaID_modified
        Photon_MVA_max = Data_tree.SubleadPhoton_mvaID_modified
    else:
        hasleadIDmin=False
        original_Photon_MVA_min = Data_tree.SubleadPhoton_mvaID_modified
        Photon_MVA_max = Data_tree.LeadPhoton_mvaID_modified
    originalminID.append(original_Photon_MVA_min)
    maxID.append(Photon_MVA_max)
    # weights.append(1)
    if(not (original_Photon_MVA_min<-0.7 and Photon_MVA_max>-0.7)):
        new_weight=-999
        minID.append(-999)
        hasMaxLead.append(-999)
    else:
        if(hasleadIDmin):
            hasMaxLead.append(0)
            LeadPhoton_mvaID_modified=photonIDPDF_fake.GetRandom(-0.7,Photon_MVA_max)
            PhotonID_min=LeadPhoton_mvaID_modified
        else:
            SubleadPhoton_mvaID_modified=photonIDPDF_fake.GetRandom(-0.7,Photon_MVA_max)
            PhotonID_min=SubleadPhoton_mvaID_modified
            hasMaxLead.append(1)
        minID.append(PhotonID_min)
        new_weight = photonIDPDF_fake.Integral(-0.7,Photon_MVA_max) / photonIDPDF_fake.Integral(-0.9,-0.7);
    weights.append(new_weight)
    if(i%100000==0):
        print("Read entry:",i,new_weight)

d={"new_weight":weights,"minID":minID,"maxID":maxID,"originalminID":originalminID,"hasMaxLead":hasMaxLead}
dataframe=pandas.DataFrame(d) 

DataFile=awkward.from_parquet("/eos/user/s/shsong/HiggsDNA/UL17data/merged_nominal.parquet")
DataFile=DataFile[DataFile.category==3]
# DataFile['year']=int(DataFile.year[0])
DataFile["Diphoton_maxID_modified"]=dataframe.maxID
DataFile["Diphoton_minID_modified"]=dataframe.minID
DataFile["originalminID"]=dataframe.originalminID
DataFile["weight_central"]=dataframe.new_weight
DataFile=DataFile[DataFile.weight_central!=-999]
#calculate the total weight
print(awkward.sum(DataFile.weight_central))


# awkward.to_parquet(DataFile,"/eos/user/s/shsong/HiggsDNA/datadriven/DatadrivenQCD_2017_cat3.parquet")
awkward.to_parquet(DataFile,"/eos/user/s/shsong/HiggsDNA/dd/DatadrivenQCD_2017_cat3.parquet")

# parquet_to_root("/eos/user/s/shsong/HiggsDNA/datadriven/DatadrivenQCD_2017_cat3.parquet","/eos/user/s/shsong/HiggsDNA/datadriven/DatadrivenQCD_2017_cat3.root",treename="cat3",verbose=False)
parquet_to_root("/eos/user/s/shsong/HiggsDNA/dd/DatadrivenQCD_2017_cat3.parquet","/eos/user/s/shsong/HiggsDNA/dd/DatadrivenQCD_2017_cat3.root",treename="cat3",verbose=False)

