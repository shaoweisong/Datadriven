from email.utils import decode_rfc2231
import awkward as ak
import numpy as np
import vector
vector.register_awkward()
import pandas as pd
from parquet_to_root import parquet_to_root
import sys
from random import choice
from math import *
import glob
import os
# inputfile=['/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/DiPhotonJetsBox_M40_80-sherpa_2017/merged_nominal.parquet',
#            '/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/DiPhotonJetsBox_MGG-80toInf_13TeV-sherpa_2017/merged_nominal.parquet',
#            '/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8_2017/merged_nominal.parquet',
#            '/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8_2017/merged_nominal.parquet',
#            '/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8_2017/merged_nominal.parquet',
#            '/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8_2017/merged_nominal.parquet',
#            '/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8_2017/merged_nominal.parquet'
#            ]
inputfile=['/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/DiPhotonJetsBox_M40_80-sherpa_2017/merged_nominal.parquet',
           '/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/DiPhotonJetsBox_MGG-80toInf_13TeV-sherpa_2017/merged_nominal.parquet'
           ]
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
    event['Diphoton_minID_modified']=photon.ID[:,1]
    event['Diphoton_maxID_modified']=photon.ID[:,0]
    return event

for file in inputfile:
    event =load_parquet(file)
    print("Got event")
    if "GJets" in file:
        event=get_minmaxID(event)
        print('Got get_minmaxID')
    if "DiPhoton" in file:
        photonID=ak.concatenate([ak.unflatten(event.LeadPhoton_mvaID_modified,counts=1),ak.unflatten(event.SubleadPhoton_mvaID_modified,counts=1)],axis=1)
        event['Diphoton_minID_modified']=photonID[ak.argsort(photonID,ascending=True)][:,0]
        event['Diphoton_maxID_modified']=photonID[ak.argsort(photonID,ascending=True)][:,1]
    outputname=file.split("/merged_nominal.parquet")[0].replace("-","_")+".parquet"
    ak.to_parquet(event, outputname)
print("start to get the gjet and pp bkg")
gjets=glob.glob("/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/GJets*.parquet")
diphotonjetsbox=glob.glob("/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/DiPhotonJetsBox*.parquet")
gjets=ak.concatenate([ak.from_parquet(file) for file in gjets])
diphotonjetsbox=ak.concatenate([ak.from_parquet(file) for file in diphotonjetsbox])
print("Got all the events")
ak.to_parquet(gjets,"/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/GJets.parquet")
ak.to_parquet(diphotonjetsbox,"/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/DiPhotonJetsBox.parquet")
gjet_cat1=gjets[gjets.category==1]
gjet_cat2=gjets[gjets.category==2]
gjet_cat3=gjets[gjets.category==3]
gjet_cat4=gjets[gjets.category==4]
ak.to_parquet(gjet_cat1,"/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/GJets_cat1.parquet")
ak.to_parquet(gjet_cat2,"/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/GJets_cat2.parquet")
ak.to_parquet(gjet_cat3,"/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/GJets_cat3.parquet")
ak.to_parquet(gjet_cat4,"/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/GJets_cat4.parquet")
parquet_to_root("/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/GJets.parquet","/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/GJets.root",treename="GJet",verbose=False)
parquet_to_root("/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/GJets_cat1.parquet","/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/GJets_cat1.root",treename="cat1",verbose=False)
parquet_to_root("/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/GJets_cat2.parquet","/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/GJets_cat2.root",treename="cat2",verbose=False)
parquet_to_root("/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/GJets_cat3.parquet","/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/GJets_cat3.root",treename="cat3",verbose=False)
parquet_to_root("/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/GJets_cat4.parquet","/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/GJets_cat4.root",treename="cat4",verbose=False)
diphotonjetsbox_cat1=diphotonjetsbox[diphotonjetsbox.category==1]
diphotonjetsbox_cat2=diphotonjetsbox[diphotonjetsbox.category==2]
diphotonjetsbox_cat3=diphotonjetsbox[diphotonjetsbox.category==3]
diphotonjetsbox_cat4=diphotonjetsbox[diphotonjetsbox.category==4]
ak.to_parquet(diphotonjetsbox_cat1,"/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/DiPhotonJetsBox_cat1.parquet")
ak.to_parquet(diphotonjetsbox_cat2,"/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/DiPhotonJetsBox_cat2.parquet")
ak.to_parquet(diphotonjetsbox_cat3,"/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/DiPhotonJetsBox_cat3.parquet")
ak.to_parquet(diphotonjetsbox_cat4,"/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/DiPhotonJetsBox_cat4.parquet")
parquet_to_root("/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/DiPhotonJetsBox.parquet","/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/DiPhotonJetsBox.root",treename="DiPhotonJetsBox",verbose=False)
parquet_to_root("/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/DiPhotonJetsBox_cat1.parquet","/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/DiPhotonJetsBox_cat1.root",treename="cat1",verbose=False)
parquet_to_root("/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/DiPhotonJetsBox_cat2.parquet","/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/DiPhotonJetsBox_cat2.root",treename="cat2",verbose=False)
parquet_to_root("/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/DiPhotonJetsBox_cat3.parquet","/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/DiPhotonJetsBox_cat3.root",treename="cat3",verbose=False)
parquet_to_root("/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/DiPhotonJetsBox_cat4.parquet","/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/DiPhotonJetsBox_cat4.root",treename="cat4",verbose=False)
print("Successfully get the gjet and pp bkg")


data=ak.from_parquet("/eos/user/s/shsong/HiggsDNA/UL17data/merged_nominal.parquet")
print("loading events:",len(data))
photonID=ak.concatenate([ak.unflatten(data.LeadPhoton_mvaID_modified,counts=1),ak.unflatten(data.SubleadPhoton_mvaID_modified,counts=1)],axis=1)
print(len(photonID))
data['Diphoton_minID_modified']=photonID[ak.argsort(photonID,ascending=True)][:,0]
data['Diphoton_maxID_modified']=photonID[ak.argsort(photonID,ascending=True)][:,1]
data_cat1=data[data.category==1]
data_cat2=data[data.category==2]
data_cat3=data[data.category==3]
data_cat4=data[data.category==4]
ak.to_parquet(data_cat1,"/eos/user/s/shsong/HiggsDNA/UL17data/Data_cat1.parquet")
ak.to_parquet(data_cat2,"/eos/user/s/shsong/HiggsDNA/UL17data/Data_cat2.parquet")
ak.to_parquet(data_cat3,"/eos/user/s/shsong/HiggsDNA/UL17data/Data_cat3.parquet")
ak.to_parquet(data_cat4,"/eos/user/s/shsong/HiggsDNA/UL17data/Data_cat4.parquet")
parquet_to_root("/eos/user/s/shsong/HiggsDNA/UL17data/Data_cat1.parquet","/eos/user/s/shsong/HiggsDNA/UL17data/Data_2017_cat1.root",treename="cat1",verbose=False)
parquet_to_root("/eos/user/s/shsong/HiggsDNA/UL17data/Data_cat2.parquet","/eos/user/s/shsong/HiggsDNA/UL17data/Data_2017_cat2.root",treename="cat2",verbose=False)
parquet_to_root("/eos/user/s/shsong/HiggsDNA/UL17data/Data_cat3.parquet","/eos/user/s/shsong/HiggsDNA/UL17data/Data_2017_cat3.root",treename="cat3",verbose=False)
parquet_to_root("/eos/user/s/shsong/HiggsDNA/UL17data/Data_cat4.parquet","/eos/user/s/shsong/HiggsDNA/UL17data/Data_2017_cat4.root",treename="cat4",verbose=False)

