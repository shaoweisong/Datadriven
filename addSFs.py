import awkward as ak
import numpy as np
import pandas as pd
import sys
from parquet_to_root import parquet_to_root
import sys
from random import choice
from math import *
import concurrent.futures

def load_parquet(fname):
    print("loading events from %s" % fname)
    events=ak.from_parquet(fname)
    events=events[events.Diphoton_minID_modified>-0.7]
    return events
def add_sale_factor(event,sclae_factor):
    event['weight_central']=sclae_factor*event.weight_central
    return event
ppsf=[1.39774,1.17699,6.94627,1.47747]
ddsf=[0.810103,1.05269,0.862671,1.30649]
for i in range(1,5):
    ppfile="/eos/user/s/shsong/HiggsDNA/FakePhotonbkg17/DiPhotonJetsBox_cat"+str(i)+".parquet"
    ddfile="/eos/user/s/shsong/HiggsDNA/datadriven/DatadrivenQCD_2017_cat"+str(i)+".parquet"
    pp=load_parquet(ppfile)
    dd=load_parquet(ddfile)
    pp_reweighted=add_sale_factor(pp,ppsf[i-1])
    dd_reweighted=add_sale_factor(dd,ddsf[i-1])
    ak.to_parquet(pp_reweighted,"/eos/user/s/shsong/HiggsDNA/datadriven/DiPhotonJetsBox_2017_cat"+str(i)+"_reweighted.parquet")
    ak.to_parquet(dd_reweighted,"/eos/user/s/shsong/HiggsDNA/datadriven/DatadrivenQCD_2017_cat"+str(i)+"_reweighted.parquet")
    parquet_to_root("/eos/user/s/shsong/HiggsDNA/datadriven/DiPhotonJetsBox_2017_cat"+str(i)+"_reweighted.parquet","/eos/user/s/shsong/HiggsDNA/datadriven/DiPhotonJetsBox_2017_cat"+str(i)+"_reweighted.root",treename="cat"+str(i),verbose=False)
    parquet_to_root("/eos/user/s/shsong/HiggsDNA/datadriven/DatadrivenQCD_2017_cat"+str(i)+"_reweighted.parquet","/eos/user/s/shsong/HiggsDNA/datadriven/DatadrivenQCD_2017_cat"+str(i)+"_reweighted.root",treename="cat"+str(i),verbose=False)
