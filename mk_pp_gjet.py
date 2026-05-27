# source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc11-opt/setup.sh

from email.utils import decode_rfc2231
import awkward as ak
import numpy as np
import vector
vector.register_awkward()
import pandas as pd
import sys
sys.path.append("/eos/user/s/shsong/pkgs_condor/parquet_to_root-0.3.0/")
from parquet_to_root import parquet_to_root
import sys
from random import choice
from math import *
import glob
import os

def load_parquet(fname):
    print("loading events from %s" % fname)
    events=ak.from_parquet(fname)
    return events
def get_minmaxID(event):
    pho_pt=ak.concatenate([ak.unflatten(event.LeadPhoton_pt,counts=1),ak.unflatten(event.SubleadPhoton_pt,counts=1)],axis=1)
    pho_eta=ak.concatenate([ak.unflatten(event.LeadPhoton_eta,counts=1),ak.unflatten(event.SubleadPhoton_eta,counts=1)],axis=1)
    pho_phi=ak.concatenate([ak.unflatten(event.LeadPhoton_phi,counts=1),ak.unflatten(event.SubleadPhoton_phi,counts=1)],axis=1)
    pho_mass=ak.concatenate([ak.unflatten(event.LeadPhoton_mass,counts=1),ak.unflatten(event.SubleadPhoton_mass,counts=1)],axis=1)
    pho_ID=ak.concatenate([ak.unflatten(event.LeadPhoton_mvaID_modified,counts=1),ak.unflatten(event.SubleadPhoton_mvaID_modified,counts=1)],axis=1)
    pho_genPartFlav=ak.concatenate([ak.unflatten(event.LeadPhoton_genPartFlav,counts=1),ak.unflatten(event.SubleadPhoton_genPartFlav,counts=1)],axis=1)
    pho_genPartIdx=ak.concatenate([ak.unflatten(event.LeadPhoton_genPartIdx,counts=1),ak.unflatten(event.SubleadPhoton_genPartIdx,counts=1)],axis=1)
    photon = ak.zip({"pt":pho_pt,"eta":pho_eta,"phi":pho_phi,"mass":pho_mass,"ID":pho_ID,"genPartFlav":pho_genPartFlav,"genPartIdx":pho_genPartIdx})
    photon=photon[ak.argsort(photon.ID,ascending=False,axis=-1)]
    event['maxIDpt']=photon.pt[:,0]
    event['maxIDeta']=photon.eta[:,0]
    event['maxIDphi']=photon.phi[:,0]
    event['maxIDmass']=photon.mass[:,0]
    event['maxID_genPartFlav']=photon.genPartFlav[:,0]
    event['maxID_genPartIdx']=photon.genPartIdx[:,0]
    event['minIDpt']=photon.pt[:,1]
    event['minIDeta']=photon.eta[:,1]
    event['minIDphi']=photon.phi[:,1]
    event['minIDmass']=photon.mass[:,1]
    event['minID_genPartFlav']=photon.genPartFlav[:,1]
    event['minID_genPartIdx']=photon.genPartIdx[:,1]
    return event

###step1, generate gjet
# inputgjet="/eos/user/s/shsong/HiggsDNA/GJet/GJet.parquet"
# inputgjet="/eos/user/s/shsong/HiggsDNA/GJet16post/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8_2016post/merged_nominal.parquet"
# if not os.path.exists(inputgjet.replace("parquet","root")):
#     gjet=load_parquet(inputgjet)
#     gjet=get_minmaxID(gjet)
#     ak.to_parquet(gjet,inputgjet)
#     parquet_to_root(inputgjet,inputgjet.replace("parquet","root"),treename="GJet",verbose=False)
#     # ak.to_parquet(gjet,inputgjet.split("GJet_Pt-40toInf_DoubleEMEnriched_")[0]+"GJet.parquet")  
#     # parquet_to_root(inputgjet.split("GJet_Pt-40toInf_DoubleEMEnriched_")[0]+"GJet.parquet" ,(inputgjet.split("GJet_Pt-40toInf_DoubleEMEnriched_")[0]+"GJet.parquet").replace("parquet","root"),treename="GJet",verbose=False)
# else:
#     print("GJet.root already exists")


###step2, generate pp  bkg
inputppfile = '/eos/cms/store/group/phys_b2g/shsong/HiggsDNA/PP16post/merged_nominal.parquet'
diphotonjetsbox=ak.from_parquet(inputppfile)
diphotonjetsbox_cat1=diphotonjetsbox[diphotonjetsbox.category==1]
diphotonjetsbox_cat2=diphotonjetsbox[diphotonjetsbox.category==2]
diphotonjetsbox_cat3=diphotonjetsbox[diphotonjetsbox.category==3]
diphotonjetsbox_cat4=diphotonjetsbox[diphotonjetsbox.category==4]
ak.to_parquet(diphotonjetsbox_cat1,"/eos/cms/store/group/phys_b2g/shsong/HiggsDNA/PP16post/DiPhotonJetsBox_cat1.parquet")
ak.to_parquet(diphotonjetsbox_cat2,"/eos/cms/store/group/phys_b2g/shsong/HiggsDNA/PP16post/DiPhotonJetsBox_cat2.parquet")
ak.to_parquet(diphotonjetsbox_cat3,"/eos/cms/store/group/phys_b2g/shsong/HiggsDNA/PP16post/DiPhotonJetsBox_cat3.parquet")
ak.to_parquet(diphotonjetsbox_cat4,"/eos/cms/store/group/phys_b2g/shsong/HiggsDNA/PP16post/DiPhotonJetsBox_cat4.parquet")
parquet_to_root("/eos/cms/store/group/phys_b2g/shsong/HiggsDNA/PP16post/DiPhotonJetsBox_cat1.parquet","/eos/cms/store/group/phys_b2g/shsong/HiggsDNA/PP16post/DiPhotonJetsBox_cat1.root",treename="cat1",verbose=False)
parquet_to_root("/eos/cms/store/group/phys_b2g/shsong/HiggsDNA/PP16post/DiPhotonJetsBox_cat2.parquet","/eos/cms/store/group/phys_b2g/shsong/HiggsDNA/PP16post/DiPhotonJetsBox_cat2.root",treename="cat2",verbose=False)
parquet_to_root("/eos/cms/store/group/phys_b2g/shsong/HiggsDNA/PP16post/DiPhotonJetsBox_cat3.parquet","/eos/cms/store/group/phys_b2g/shsong/HiggsDNA/PP16post/DiPhotonJetsBox_cat3.root",treename="cat3",verbose=False)
parquet_to_root("/eos/cms/store/group/phys_b2g/shsong/HiggsDNA/PP16post/DiPhotonJetsBox_cat4.parquet","/eos/cms/store/group/phys_b2g/shsong/HiggsDNA/PP16post/DiPhotonJetsBox_cat4.root",treename="cat4",verbose=False)
print("Successfully get the gjet and pp bkg")


data=ak.from_parquet("/eos/user/s/shsong/HiggsDNA/UL16postdata/merged_nominal.parquet")
print("loading events:",len(data))
data_cat1=data[data.category==1]
data_cat2=data[data.category==2]
data_cat3=data[data.category==3]
data_cat4=data[data.category==4]
ak.to_parquet(data_cat1,"/eos/user/s/shsong/HiggsDNA/UL16postdata/Data_cat1.parquet")
ak.to_parquet(data_cat2,"/eos/user/s/shsong/HiggsDNA/UL16postdata/Data_cat2.parquet")
ak.to_parquet(data_cat3,"/eos/user/s/shsong/HiggsDNA/UL16postdata/Data_cat3.parquet")
ak.to_parquet(data_cat4,"/eos/user/s/shsong/HiggsDNA/UL16postdata/Data_cat4.parquet")
parquet_to_root("/eos/user/s/shsong/HiggsDNA/UL16postdata/Data_cat1.parquet","/eos/user/s/shsong/HiggsDNA/UL16postdata/Data_2016post_cat1.root",treename="cat1",verbose=False)
parquet_to_root("/eos/user/s/shsong/HiggsDNA/UL16postdata/Data_cat2.parquet","/eos/user/s/shsong/HiggsDNA/UL16postdata/Data_2016post_cat2.root",treename="cat2",verbose=False)
parquet_to_root("/eos/user/s/shsong/HiggsDNA/UL16postdata/Data_cat3.parquet","/eos/user/s/shsong/HiggsDNA/UL16postdata/Data_2016post_cat3.root",treename="cat3",verbose=False)
parquet_to_root("/eos/user/s/shsong/HiggsDNA/UL16postdata/Data_cat4.parquet","/eos/user/s/shsong/HiggsDNA/UL16postdata/Data_2016post_cat4.root",treename="cat4",verbose=False)

