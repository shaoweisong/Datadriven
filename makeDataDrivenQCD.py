import ROOT
from array import array
import awkward 
from parquet_to_root import parquet_to_root
import pandas
#run mk_pp_gjet.py firstly to get the Data_2017.root, GJet.root
# GJets1=ROOT.TFile.Open("/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/GJets.root")
GJets1=ROOT.TFile.Open("/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/GJets.root")
# GJets1_tree=GJets1.Get("GJet")
GJets1_tree=GJets1.Get("GJet")
cuts="(Diphoton_mass<115 || Diphoton_mass>135)*(minID_genPartFlav!=1)"
h_minphotonID=ROOT.TH1F("h_minphotonID_gjet","h_minphotonID_gjet",19,-0.9,1)
GJets1_tree.Project("h_minphotonID_gjet","Diphoton_minID","weight_central*"+cuts)
photonIDPDF_fake=ROOT.TF1("photonIDPDF_fake","pol7",-0.9,1.)
h_minphotonID.Fit(photonIDPDF_fake,"R")
c1=ROOT.TCanvas("c1","c1",600,800)
h_minphotonID.Draw("E1")
c1.SaveAs("fakephoton_pdf_cat1.png")
Data=ROOT.TFile.Open("/eos/user/s/shsong/HiggsDNA/UL17data/Data_2017_cat1.root")
Data_tree=Data.Get("cat1")
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

print(sum(weights))
d={"new_weight":weights,"minID":minID,"maxID":maxID,"originalminID":originalminID,"hasMaxLead":hasMaxLead}
dataframe=pandas.DataFrame(d) 

DataFile=awkward.from_parquet("/eos/user/s/shsong/HiggsDNA/UL17data/merged_nominal.parquet")
DataFile=DataFile[DataFile.category==1]
DataFile["Diphoton_maxID_modified"]=dataframe.maxID
DataFile["Diphoton_minID_modified"]=dataframe.minID
DataFile["originalminID"]=dataframe.originalminID
DataFile["weight_central"]=dataframe.new_weight
DataFile=DataFile[DataFile.weight_central!=-999]

awkward.to_parquet(DataFile,"/eos/user/s/shsong/HiggsDNA/datadriven/DatadrivenQCD_2017_cat1.parquet")

parquet_to_root("/eos/user/s/shsong/HiggsDNA/datadriven/DatadrivenQCD_2017_cat1.parquet","/eos/user/s/shsong/HiggsDNA/datadriven/DatadrivenQCD_2017_cat1.root",treename="cat1",verbose=False)

